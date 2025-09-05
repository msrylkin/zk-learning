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
}

fn setup<P: Pairing>(n: usize, toxic_waste: P::ScalarField) -> Config<P> {
    Config {
        srs: (0..n).map(|i| P::G1::generator() * toxic_waste.pow([i as u64])).collect(),
        verifier_key: P::G2::generator() * toxic_waste,
    }
}

struct Opening<P: Pairing> {
    evaluation: P::ScalarField,
    evaluation_proof: P::G1,
}

impl<P: Pairing> Opening<P> {
    pub fn new(evaluation: P::ScalarField, evaluation_proof: P::G1) -> Self {
        Self {
            evaluation,
            evaluation_proof,
        }
    }
}

struct BatchOpening<P: Pairing> {
    evaluations: Vec<P::ScalarField>,
    evaluation_proof: P::G1,
}

impl<P: Pairing> BatchOpening<P> {
    pub fn new(evaluations: Vec<P::ScalarField>, evaluation_proof: P::G1) -> Self {
        Self {
            evaluations,
            evaluation_proof,
        }
    }
}

impl<P: Pairing> From<Opening<P>> for BatchOpening<P> {
    fn from(opening: Opening<P>) -> Self {
        Self {
            evaluations: vec![opening.evaluation],
            evaluation_proof: opening.evaluation_proof,
        }
    }
}

impl<P: Pairing> KZG<P> {
    pub fn new(config: Config<P>) -> Self {
        Self { config }
    }

    pub fn commit(
        &self,
        poly: &DensePolynomial<P::ScalarField>,
    ) -> P::G1 {
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

        BatchOpening::new(evaluations, eval_proof)
    }

    pub fn open(
        &self,
        poly: &DensePolynomial<P::ScalarField>,
        point: &P::ScalarField,
    ) -> Opening<P> {
        let BatchOpening { evaluations, evaluation_proof } = self.batch_open(
            &[poly],
            point,
            &P::ScalarField::one(),
        );

        Opening::new(evaluations[0], evaluation_proof)
    }

    pub fn check_single(
        &self,
        point: &P::ScalarField,
        commitment: &P::G1,
        opening: &Opening<P>,
    ) -> bool {
        self.check(
            point,
            &[*commitment],
            &[opening.evaluation],
            &opening.evaluation_proof,
            &P::ScalarField::one(),
        )
    }

    pub fn check_batched(
        &self,
        point: &P::ScalarField,
        commitments: &[P::G1],
        openings: &BatchOpening<P>,
        linearization_scalar: &P::ScalarField,
    ) -> bool {
        self.check(
            point,
            commitments,
            &openings.evaluations,
            &openings.evaluation_proof,
            linearization_scalar,
        )
    }

    fn check(
        &self,
        point: &P::ScalarField,
        commitments: &[P::G1],
        evaluations: &[P::ScalarField],
        evaluation_proof: &P::G1,
        linearization_scalar: &P::ScalarField,
    ) -> bool {
        let linearized_commitments = Self::linearize_commitments(commitments, linearization_scalar);
        let linearized_evaluations = Self::linearize_elements(evaluations, linearization_scalar);

        P::pairing(
            linearized_commitments - P::G1::generator() * linearized_evaluations,
            P::G2::generator(),
        ) == P::pairing(
            evaluation_proof,
            self.config.verifier_key - P::G2::generator() * point,
        )
    }
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
    use crate::kzg::{setup, BatchOpening, KZG};
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
        let config = setup::<Bls12_381>(11, tau);
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
            opening.evaluation,
            poly.evaluate(&opening_point),
        );
        assert_eq!(
            opening.evaluation_proof,
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
            opening.evaluations,
            openings,
        );
        assert_eq!(
            opening.evaluation_proof,
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

        let is_valid = kzg.check_single(
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