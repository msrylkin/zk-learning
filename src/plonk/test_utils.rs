use ark_ec::pairing::Pairing;
use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;
use ark_test_curves::bls12_381::{Bls12_381, Fr};
use crate::kzg::{setup, KZG};
use crate::poly_utils::const_poly;

#[cfg(test)]
pub fn get_test_kzg<P: Pairing>(n: usize) -> KZG<P> {
    let tau = P::ScalarField::from(777);
    let config = setup::<P>(n * 2, tau);
    let kzg = KZG::new(config);

    kzg
}

#[cfg(test)]
pub fn hash_permutation_poly<F: Field>(
    f: &DensePolynomial<F>,
    perm_poly: &DensePolynomial<F>,
    beta: F,
    gamma: F,
) -> DensePolynomial<F> {
    f + perm_poly * beta + const_poly(gamma)
}