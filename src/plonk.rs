mod prover;
mod circuit;
mod permutation;
mod evaluation_domain;
mod blinder;
mod proof;
mod srs;
mod verifier;

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::{Bls12_381, Fr};
    use crate::kzg::{setup, KZG};
    use crate::plonk::circuit::get_test_circuit;
    use crate::plonk::prover::prove;
    use crate::plonk::verifier::verify;
    use crate::poly_utils::generate_multiplicative_subgroup;

    #[test]
    pub fn protocol_test() {
        let circuit = get_test_circuit();
        let tau = Fr::from(777);
        let domain = generate_multiplicative_subgroup::<{ 1 << 3 }, Fr>();
        let config = setup::<Bls12_381>(domain.len() * 2, tau);
        let kzg = KZG::new(config);
        let proof = prove(
            &circuit,
            &kzg,
            &domain,
        );

        let is_valid = verify(
            &circuit,
            &domain,
            &kzg,
            proof,
        );
    }
}