mod prover;
mod circuit;
mod permutation;
mod blinder;
mod proof;
mod srs;
mod verifier;

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::{Bls12_381, Fr};
    use crate::kzg::{setup, KZG};
    use crate::plonk::circuit::{get_test_circuit, get_test_solution};
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::plonk::prover::prove;
    use crate::plonk::verifier::verify;

    #[test]
    pub fn protocol_test() {
        let domain = generate_multiplicative_subgroup::<{ 1 << 3 }, Fr>();
        let circuit = get_test_circuit(&domain);
        let solution = get_test_solution(&domain);
        let tau = Fr::from(777);
        let config = setup::<Bls12_381>(domain.len() * 2, tau);
        let kzg = KZG::new(config);
        let pi = solution.public_input.pi.clone();
        let proof = prove(
            &circuit,
            solution,
            &kzg,
            &domain,
        );

        let is_valid = verify(
            &circuit,
            &pi,
            &domain,
            &kzg,
            proof,
        );
    }
}