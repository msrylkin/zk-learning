use ark_ec::pairing::Pairing;
use ark_test_curves::bls12_381::{Bls12_381, Fr};
use crate::kzg::{setup, KZG};

#[cfg(test)]
pub fn get_test_kzg<P: Pairing>(n: usize) -> KZG<P> {
    let tau = P::ScalarField::from(777);
    let config = setup::<P>(n * 2, tau);
    let kzg = KZG::new(config);

    kzg
}