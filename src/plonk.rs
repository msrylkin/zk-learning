mod prover;
mod circuit;
mod permutation;
mod blinder;
mod proof;
mod verifier;
mod transcript_protocol;
mod protocol;
mod domain;
mod test_utils;
mod big_quotient_poly;

pub use circuit::*;
pub use domain::*;
pub use protocol::*;
pub use proof::*;

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::{Bls12_381, Fr};
    use crate::kzg::{setup, KZG};
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::plonk::domain::PlonkDomain;
    use crate::plonk::protocol::{PlonkProtocol};
    use crate::plonk::test_utils::{get_test_circuit, get_test_solution};

    #[test]
    pub fn protocol_test() {
        let domain = generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);
        let circuit = get_test_circuit(&domain);

        let solution = get_test_solution(&domain);
        let tau = Fr::from(777);
        let config = setup::<Bls12_381>(domain.len() * 2, tau);
        let kzg = KZG::new(config);
        let public_input = solution.public_witness.clone();

        let protocol = PlonkProtocol::new(kzg, domain);
        let instance = protocol.create_instance(&circuit);

        let proof = instance.prove(solution);
        instance.verify(&proof, &public_input);
    }
}