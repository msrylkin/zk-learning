mod prover;
mod circuit;
mod permutation;
mod blinder;
mod proof;
mod srs;
mod verifier;
mod transcript_protocol;
mod protocol;
mod domain;
mod test_utils;
mod big_quotient_poly;

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::{Bls12_381, Fr};
    use crate::kzg::{setup, KZG};
    use crate::plonk::circuit::{get_test_circuit, get_test_solution};
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::plonk::domain::PlonkDomain;
    use crate::plonk::prover::PlonkProver;
    use crate::plonk::verifier::PlonkVerifier;
    // use crate::plonk::verifier::verify;

    #[test]
    pub fn protocol_test() {
        let domain = generate_multiplicative_subgroup::<{ 1 << 3 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);
        let circuit = get_test_circuit(&domain);
        let solution = get_test_solution(&domain);
        let tau = Fr::from(777);
        let config = setup::<Bls12_381>(domain.len() * 2, tau);
        let kzg = KZG::new(config);
        let public_input = solution.public_input.clone();

        let prover = PlonkProver::new(
            &kzg,
            &domain,
            &circuit,
        );
        
        let proof = prover.prove(solution);
        
        let verifier = PlonkVerifier::new(
            &kzg,
            &domain,
            &circuit,
        );

        verifier.verify(
            &public_input,
            &proof
        );
    }
}