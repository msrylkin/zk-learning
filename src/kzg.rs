use std::ops::{Div, Mul, Sub};
use ark_ec::{AffineRepr, CurveGroup, PrimeGroup};
use ark_ec::pairing::Pairing;
use ark_std::{UniformRand, Zero};
use ark_test_curves::bls12_381::{Bls12_381, Fr, G1Affine as G1, G2Affine as G2};
use ark_ff::{FftField, One, PrimeField};
use ark_poly::univariate::{DensePolynomial};
use ark_ff::Field;
use ark_poly::{DenseUVPolynomial, Polynomial};

pub struct Config<P: Pairing> {
    srs: Vec<P::G1>,
    pub verifier_key: P::G2,
}

pub struct KZG<P: Pairing> {
    pub config: Config<P>,
}

pub struct MultipointOpening<'a, P: Pairing> {
    pub batch_opening: &'a BatchOpening<P>,
    pub linearization_scalar: &'a P::ScalarField,
    pub opening_point: &'a P::ScalarField,
    pub commitments: &'a [P::G1],
}

pub fn setup<P: Pairing>(n: usize, toxic_waste: P::ScalarField) -> Config<P> {
    Config {
        srs: (0..n).map(|i| P::G1::generator() * toxic_waste.pow([i as u64])).collect(),
        verifier_key: P::G2::generator() * toxic_waste,
    }
}

#[derive(Debug)]
pub struct Opening<P: Pairing> {
    pub evaluation: P::ScalarField,
    pub evaluation_proof: P::G1,
}

impl<P: Pairing> Opening<P> {
    pub fn new(evaluation: P::ScalarField, evaluation_proof: P::G1) -> Self {
        Self {
            evaluation,
            evaluation_proof,
        }
    }
}

pub struct BatchOpening<P: Pairing> {
    evaluations: Vec<P::ScalarField>,
    evaluation_proof: P::G1,
}

impl<P: Pairing> BatchOpening<P> {
    pub fn proof(&self) -> P::G1 {
        self.evaluation_proof.clone()
    }
}

impl<P: Pairing> Opening<P> {
    pub fn proof(&self) -> P::G1 {
        self.evaluation_proof.clone()
    }
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

impl<P: Pairing> From<&Opening<P>> for BatchOpening<P> {
    fn from(opening: &Opening<P>) -> Self {
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
        self.check_multipoint(
            &[
                MultipointOpening {
                    batch_opening: &BatchOpening::from(opening),
                    linearization_scalar: &P::ScalarField::one(),
                    opening_point: point,
                    commitments: &[*commitment],
                }
            ],
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
        self.check_multipoint(
            &[
                MultipointOpening {
                    batch_opening: openings,
                    linearization_scalar,
                    opening_point: point,
                    commitments,
                }
            ],
            &P::ScalarField::one(),
        )
    }

    fn linearized_comm_evals_diff(
        commitments: &[P::G1],
        evaluations: &[P::ScalarField],
        linearization_scalar: &P::ScalarField,
    ) -> P::G1 {
        let linearized_commitments = Self::linearize_commitments(commitments, linearization_scalar);
        let linearized_evaluations = Self::linearize_elements(evaluations, linearization_scalar);

        linearized_commitments - P::G1::generator() * linearized_evaluations
    }

    pub fn check_multipoint(
        &self,
        // commitments: &[P::G1],
        multipoint_openings: &[MultipointOpening<P>],
        multipoint_linearization_scalar: &P::ScalarField,
    ) -> bool {
        let (l1_pair, l2_pair) = multipoint_openings
            .iter()
            .enumerate()
            .fold((P::G1::zero(), P::G1::zero()), |(l1_pair, l2_pair), (i, opening) | {
                let diff = Self::linearized_comm_evals_diff(
                    opening.commitments,
                    &opening.batch_opening.evaluations,
                    &opening.linearization_scalar,
                );

                let scale = multipoint_linearization_scalar.pow([i as u64]);

                (
                    l1_pair + diff * scale + opening.batch_opening.evaluation_proof * opening.opening_point * scale,
                    l2_pair + opening.batch_opening.evaluation_proof * scale,
                )
            });

        P::pairing(
            l1_pair,
            P::G2::generator(),
        ) == P::pairing(
            l2_pair,
            self.config.verifier_key,
        )
    }
}

#[cfg(test)]
mod tests {
    use std::ops::Mul;
    use ark_ec::{CurveGroup, PrimeGroup};
    use ark_ec::pairing::Pairing;
    use ark_ff::Field;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use ark_poly::univariate::DensePolynomial;
    use ark_test_curves::bls12_381::{Bls12_381, Fr, G1Projective as G1,};
    use crate::kzg::{setup, BatchOpening, MultipointOpening, KZG};
    use crate::poly_utils::to_f;

    struct TestData {
        kzg: KZG<Bls12_381>,
        tau: Fr,
    }

    fn test_setup() -> TestData {
        let tau = Fr::from(333444555);
        let config = setup::<Bls12_381>(30, tau);
        let kzg = KZG::<Bls12_381>::new(config);

        TestData {
            kzg,
            tau,
        }
    }

    fn scalars(i: usize) -> (Fr, Fr) {
        let ii = Fr::from(i as u64);
        let linearization_scalar = (Fr::from(133) * ii + ii).pow([i as u64]);
        let opening_point = (Fr::from(123) * ii + ii).pow([i as u64]);

        (linearization_scalar, opening_point)
    }

    fn test_polys() -> (DensePolynomial<Fr>, DensePolynomial<Fr>, DensePolynomial<Fr>, DensePolynomial<Fr>, DensePolynomial<Fr>) {
        let poly_1 = DensePolynomial::from_coefficients_vec(
            to_f::<Fr>(vec![2, 3, 4, 5, 6, 7, 8, 9, 10, 11]),
        );
        let poly_2 = poly_1.clone() * Fr::from(123) * DensePolynomial::from_coefficients_vec(vec![
            Fr::from(999), Fr::from(8919)
        ]);
        let poly_3 = DensePolynomial::from_coefficients_vec(
            to_f::<Fr>(vec![22, 33, 44, 55, 66, 77, 88, 99, 1000, 1231]),
        );
        let poly_4 = DensePolynomial::from_coefficients_vec(
            to_f::<Fr>(vec![222, 333, 444, 5555, 6676, 717, 88, 9129, 101100, 1231231]),
        );
        let poly_5 = poly_4.clone() * poly_3.clone() * Fr::from(33) + DensePolynomial::from_coefficients_slice(&[Fr::from(4455)]);

        (poly_1, poly_2, poly_3, poly_4, poly_5)
    }

    #[test]
    pub fn test_commitment() {
        let TestData { kzg, tau, .. } = test_setup();
        let (poly, ..) = test_polys();

        let comm = kzg.commit(&poly);

        assert_eq!(comm, G1::generator() * poly.evaluate(&tau));
    }

    #[test]
    pub fn test_batch_commitment() {
        let TestData { kzg, tau, .. } = test_setup();
        let (linearization_scalar, .. ) = scalars(1);
        let (poly_1, poly_2, poly_3, ..) = test_polys();

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
        let TestData { kzg, tau, .. } = test_setup();
        let (_, opening_point ) = scalars(1);

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
        let TestData { kzg, tau, .. } = test_setup();
        let (linearization_scalar, opening_point ) = scalars(1);
        let (poly_1, poly_2, poly_3, ..) = test_polys();

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
        let (poly, ..) = test_polys();
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
        let TestData { kzg, .. } = test_setup();
        let (linearization_scalar, opening_point ) = scalars(1);

        let (poly_1, poly_2, poly_3, ..) = test_polys();
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

    #[test]
    pub fn test_kzg_batch_open_multipoint() {
        let TestData { kzg, .. } = test_setup();
        let (linearization_scalar_1, opening_point_1 ) = scalars(1);
        let (linearization_scalar_2, opening_point_2 ) = scalars(2);

        let (poly_1, poly_2, poly_3, poly_4, poly_5) = test_polys();

        let commitment_1 = kzg.commit(&poly_1);
        let commitment_2 = kzg.commit(&poly_2);
        let commitment_3 = kzg.commit(&poly_3);
        let commitment_4 = kzg.commit(&poly_4);
        let commitment_5 = kzg.commit(&poly_5);

        let openings_1 = kzg.batch_open(
            &[&poly_1, &poly_2, &poly_3],
            &opening_point_1,
            &linearization_scalar_1,
        );
        let openings_2 = kzg.batch_open(
            &[&poly_4, &poly_5,],
            &opening_point_2,
            &linearization_scalar_2,
        );

        let is_valid = kzg.check_multipoint(
            &[
                MultipointOpening {
                    batch_opening: &openings_1,
                    linearization_scalar: &linearization_scalar_1,
                    opening_point: &opening_point_1,
                    commitments: &[commitment_1, commitment_2, commitment_3]
                },
                MultipointOpening {
                    batch_opening: &openings_2,
                    linearization_scalar: &linearization_scalar_2,
                    opening_point: &opening_point_2,
                    commitments: &[commitment_4, commitment_5]
                }
            ],
            &Fr::from(777),
        );

        assert!(is_valid);
    }
}