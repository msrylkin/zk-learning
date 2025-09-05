use std::ops::{Div, Mul, Sub};
use ark_ec::{AffineRepr, CurveGroup, PrimeGroup};
use ark_ec::pairing::Pairing;
use ark_std::{UniformRand, Zero};
use ark_test_curves::bls12_381::{Bls12_381, Fr, G1Affine as G1, G2Affine as G2};
use ark_ff::{FftField, One, PrimeField};
use ark_poly::univariate::{DensePolynomial};
use ark_ec::short_weierstrass::{Affine};
use ark_ff::Field;
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_std::iterable::Iterable;
use ark_test_curves::bls12_381::g1::Config as G1Config;
use ark_test_curves::bls12_381::g2::Config as G2Config;

struct Config<P: Pairing> {
    srs: Vec<P::G1>,
    verifier_key: P::G2,
}

struct KZG<P: Pairing> {
    config: Config<P>,
    G1: P::G1,
    G2: P::G2,
}

fn setup_2<P: Pairing>(n: usize, toxic_waste: P::ScalarField) -> Config<P> {
    // let toxic_waste = P::ScalarField::from(333444555);

    Config {
        srs: (0..n).map(|i| P::G1::generator() * toxic_waste.pow([i as u64])).collect(),
        verifier_key: P::G2::generator() * toxic_waste,
    }
}

type Opening<P> = (<P as Pairing>::ScalarField, <P as Pairing>::G1);

// struct Opening<P: Pairing>(P::ScalarField, P::G1);

// struct Opening<P: Pairing> {
//     opening: P::ScalarField,
//     evaluation_proof: P::G1,
// }
//
// struct BatchOpening<P: Pairing> {
//     openings: Vec<Opening<P>>,
//     evaluation_proof: P::G1,
// }

type BatchOpening<P> = (Vec<<P as Pairing>::ScalarField>, <P as Pairing>::G1);

impl<P: Pairing> KZG<P> {
    pub fn new(config: Config<P>) -> Self {
        Self { config, G1: P::G1::generator(), G2: P::G2::generator() }
    }

    pub fn commit(
        &self,
        poly: &DensePolynomial<P::ScalarField>,
    ) -> P::G1 {
        println!("\nsrs.len() {}", self.config.srs.len());
        println!("poly degree {}", poly.degree());
        println!("poly coefs {}", poly.coeffs.len());
        if self.config.srs.len() < poly.degree() + 1 {
            panic!("srs lower than poly degree");
        }

        poly.coeffs
            .iter()
            .zip(&self.config.srs)
            .map(|(x_i, h)| *h * x_i)
            .fold(P::G1::zero(), |a, b| a + b)
    }

    pub fn linearize_elements(elements: &[P::ScalarField], linearization_scalar: &P::ScalarField) -> P::ScalarField {
        elements
            .iter()
            .enumerate()
            .fold(P::ScalarField::zero(), |acc, (i, e)| {
                acc + *e * linearization_scalar.pow([i as u64])
            })
    }

    pub fn linearize_commitments(elements: &[P::G1], linearization_scalar: &P::ScalarField) -> P::G1 {
        elements
            .iter()
            .enumerate()
            .fold(P::G1::zero(), |acc, (i, e)| {
                println!("i {i} {}", (*e * linearization_scalar.pow([i as u64])).into_affine());
                println!("i {i} comm {}", (*e).into_affine());
                println!("lc {i} {}", linearization_scalar.pow([i as u64]));

                acc + *e * linearization_scalar.pow([i as u64])
            })
    }

    pub fn batch_open(
        &self,
        polys: &[&DensePolynomial<P::ScalarField>],
        point: &P::ScalarField,
        linearization_scalar: &P::ScalarField,
    ) -> BatchOpening<P> {
        let mut linearized_poly = DensePolynomial::zero();
        let mut evaluations = vec![];

        for (i, poly) in polys.iter().enumerate() {
            let evaluation = poly.evaluate(point);
            evaluations.push(evaluation);
            let quotient_poly = (*poly - DensePolynomial::from_coefficients_slice(&[evaluation]))
                / DensePolynomial::from_coefficients_slice(&[-*point, P::ScalarField::one()]);

            linearized_poly = linearized_poly + quotient_poly * linearization_scalar.pow([i as u64]);
        }

        let eval_proof = self.commit(&linearized_poly);

        (evaluations, eval_proof)
    }

    pub fn open(
        &self,
        poly: &DensePolynomial<P::ScalarField>,
        point: &P::ScalarField,
    ) -> Opening<P> {
        let (evaluations, eval_proof) = self.batch_open(
            &[poly],
            point,
            &P::ScalarField::one(),
        );

        (evaluations[0], eval_proof)
    }

    pub fn check(
        &self,
        point: &P::ScalarField,
        commitment: &P::G1,
        eval_proof: &Opening<P>,
    ) -> bool {
        self.check_batched(
            point,
            &[*commitment],
            &(vec![eval_proof.0], eval_proof.1),
            &P::ScalarField::one(),
        )
    }

    pub fn check_batched(
        &self,
        point: &P::ScalarField,
        commitments: &[P::G1],
        eval_proofs: &BatchOpening<P>,
        linearization_scalar: &P::ScalarField,
    ) -> bool {
        let linearized_commitments = Self::linearize_commitments(commitments, linearization_scalar);
        let linearized_evaluations = Self::linearize_elements(&eval_proofs.0, linearization_scalar);

        println!("real linearized_commitments: {:?}", linearized_commitments.into_affine());
        println!("combined evals {:?}", (P::G1::generator() * linearized_evaluations).into_affine());
        let real_p1_left = linearized_commitments - P::G1::generator() * linearized_evaluations;
        println!("\nreal p1_left {:?}", real_p1_left.into_affine());
        println!("real p1_right {:?}", P::G2::generator().into_affine());

        println!("real p2_left {:?}", eval_proofs.1.into_affine());
        println!("real p2_right {:?}", (self.config.verifier_key - P::G2::generator() * point).into_affine());

        P::pairing(
            linearized_commitments - P::G1::generator() * linearized_evaluations,
            P::G2::generator(),
        ) == P::pairing(
            eval_proofs.1,
            self.config.verifier_key - P::G2::generator() * point,
        )
    }
}

pub fn run_kzg() {
    let polynomial = DensePolynomial::from_coefficients_vec(
        vec![Fr::from(11), Fr::from(12), Fr::from(345)]
    );

    let (srs, tau_g2) = setup();
    let commitment = commit_to_poly(&polynomial, &srs);

    let random_a = Fr::from(777);

    let (f_a, evaluation_proof) = open_commitment(&srs, &random_a, &polynomial);

    let is_valid = verify_commitment(&commitment, &tau_g2, &random_a, &f_a, &evaluation_proof);

    println!("is_valid = {}", is_valid);
}

fn setup() -> (Vec<Affine<G1Config>>, Affine<G2Config>) {
    let max_degree = 8;
    let rng = &mut ark_std::test_rng();

    let tau = Fr::rand(rng);

    let tau_g2 = (G2::generator() * tau).into_affine();
    let srs = (0..max_degree).map(|d| {
        let tau_times_d  = tau.pow(&[d]);
        G1::generator().mul(tau_times_d).into_affine()
    }).collect::<Vec<_>>();

    (srs, tau_g2)
}

fn commit_to_poly(
    poly: &DensePolynomial<Fr>,
    srs: &Vec<Affine<G1Config>>
) -> Affine<G1Config> {
    let comm_f = eval_poly_on_srs(poly, srs);

    comm_f
}

fn open_commitment(
    srs: &Vec<Affine<G1Config>>,
    a: &Fr,
    poly: &DensePolynomial<Fr>,
) -> (Fr, Affine<G1Config>) {
    let f_a = poly.evaluate(a);

    let f_a_poly = DensePolynomial::from_coefficients_vec(vec![f_a]);
    let x_minus_a_poly = DensePolynomial::from_coefficients_vec(vec![-*a, Fr::one()]);
    let q_x = poly.sub(&f_a_poly).div(&x_minus_a_poly);

    let evaluation_proof = eval_poly_on_srs(&q_x, srs);

    (f_a, evaluation_proof)
}

fn verify_commitment(
    commitment: &Affine<G1Config>,
    tau_g2: &Affine<G2Config>,
    a: &Fr,
    f_a: &Fr,
    evaluation_proof: &Affine<G1Config>,
) -> bool {
    let l_pairing = Bls12_381::pairing(
        evaluation_proof,
        *tau_g2 - G2::generator() * a,
    );
    let r_pairing = Bls12_381::pairing(
        *commitment - G1::generator() * f_a,
        G2::generator(),
    );

    l_pairing == r_pairing
}

fn eval_poly_on_srs(
    poly: &DensePolynomial<Fr>,
    srs: &[Affine<G1Config>]
) -> Affine<G1Config> {
    if srs.len() < poly.degree() {
        panic!("srs lower than poly degree");
    }

    poly.coeffs
        .iter()
        .zip(srs)
        .map(|(x_i, h)| h.mul(x_i))
        .reduce(|a, b| a + b)
        .unwrap()
        .into_affine()
}

#[cfg(test)]
mod tests {
    use ark_ec::{CurveGroup, PrimeGroup};
    use ark_ec::pairing::Pairing;
    use ark_ff::Field;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use ark_poly::univariate::DensePolynomial;
    use ark_std::{One, Zero};
    use ark_test_curves::bls12_381::{Bls12_381, Fr, G1Projective as G1, G2Projective as G2};
    use crate::kzg::lib::{commit_to_poly, setup_2, BatchOpening, KZG};
    use crate::poly_utils::to_f;

    struct TestData {
        kzg: KZG<Bls12_381>,
        tau: Fr,
        linearization_scalar: Fr,
        opening_point: Fr,
    }

    fn test_setup() -> TestData {
        let linearization_scalar = Fr::from(133);
        let tau = Fr::from(333444555);
        let config = setup_2::<Bls12_381>(11, tau);
        let kzg = KZG::<Bls12_381>::new(config);
        let opening_point = Fr::from(123);

        TestData {
            kzg,
            tau,
            linearization_scalar,
            opening_point
        }
    }

    fn test_polys() -> (DensePolynomial<Fr>, DensePolynomial<Fr>, DensePolynomial<Fr>) {
        let poly_1 = DensePolynomial::from_coefficients_vec(
            to_f::<Fr>(vec![2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
        );
        let poly_2 = poly_1.clone() * Fr::from(123) * DensePolynomial::from_coefficients_vec(vec![
            Fr::from(999), Fr::from(8919)
        ]);
        let poly_3 = DensePolynomial::from_coefficients_vec(
            to_f::<Fr>(vec![22, 33, 44, 55, 66, 77, 88, 99, 1000, 1231]),
        );

        (poly_1, poly_2, poly_3)
    }

    #[test]
    pub fn test_commitment() {
        let TestData { kzg, tau, .. } = test_setup();
        let (poly, _, _) = test_polys();

        let comm = kzg.commit(&poly);

        assert_eq!(comm, G1::generator() * poly.evaluate(&tau));
    }

    #[test]
    pub fn test_batch_commitment() {
        let TestData { kzg, tau, linearization_scalar, .. } = test_setup();
        let (poly_1, poly_2, poly_3) = test_polys();

        let comm_1 = kzg.commit(&poly_1);
        let comm_2 = kzg.commit(&poly_2);
        let comm_3 = kzg.commit(&poly_3);
        let comm = KZG::<Bls12_381>::linearize_commitments(
            &[comm_1, comm_2, comm_3],
            &linearization_scalar,
        );

        assert_eq!(
            comm,
            G1::generator() * (poly_1.evaluate(&tau) + poly_2.evaluate(&tau) * linearization_scalar + poly_3.evaluate(&tau) * linearization_scalar.square())
        );
    }

    #[test]
    pub fn test_open() {
        let TestData { kzg, tau, opening_point, .. } = test_setup();
        let (poly, ..) = test_polys();

        let opening = kzg.open(&poly, &opening_point);

        assert_eq!(
            opening.0,
            poly.evaluate(&opening_point),
        );
        assert_eq!(
            opening.1,
            G1::generator() * ((poly.evaluate(&tau) - poly.evaluate(&opening_point)) / (tau - opening_point))
        )
    }

    #[test]
    pub fn test_batch_open() {
        let TestData { kzg, tau, linearization_scalar, opening_point } = test_setup();
        let (poly_1, poly_2, poly_3) = test_polys();

        let opening = kzg.batch_open(
            &[&poly_1, &poly_2, &poly_3],
            &opening_point,
            &linearization_scalar,
        );

        let openings = &[poly_1.evaluate(&opening_point), poly_2.evaluate(&opening_point), poly_3.evaluate(&opening_point)];
        let linearized_commitments_raw = poly_1.evaluate(&tau)
            + poly_2.evaluate(&tau) * linearization_scalar
            + poly_3.evaluate(&tau) * linearization_scalar.square();
        let linearized_openings_raw = poly_1.evaluate(&opening_point) +
            poly_2.evaluate(&opening_point) * linearization_scalar
            + poly_3.evaluate(&opening_point) * linearization_scalar.square();

        assert_eq!(
            opening.0,
            openings,
        );
        assert_eq!(
            opening.1,
            G1::generator() * ((linearized_commitments_raw - linearized_openings_raw) / (tau - opening_point))
        );
    }

    #[test]
    pub fn test_kzg_single_open() {
        let TestData { kzg, .. } = test_setup();
        let (poly, _, _) = test_polys();
        let commitment = kzg.commit(&poly);
        let point = Fr::from(123);
        let opening = kzg.open(&poly, &point);

        let is_valid = kzg.check(
            &point,
            &commitment,
            &opening,
        );

        assert!(is_valid);
    }

    #[test]
    pub fn test_kzg_batch_open() {
        let TestData { kzg, linearization_scalar, opening_point, .. } = test_setup();
        let (poly_1, poly_2, poly_3) = test_polys();
        let commitment_1 = kzg.commit(&poly_1);
        let commitment_2 = kzg.commit(&poly_2);
        let commitment_3 = kzg.commit(&poly_3);

        let openings = kzg.batch_open(
            &[&poly_1, &poly_2, &poly_3],
            &opening_point,
            &linearization_scalar,
        );

        let is_valid = kzg.check_batched(
            &opening_point,
            &[commitment_1, commitment_2, commitment_3],
            &openings,
            &linearization_scalar
        );

        assert!(is_valid);
    }
}