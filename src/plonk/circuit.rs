pub mod circuit_description;
pub mod solution;
pub mod compiled_circuit;
mod preprocess;

use ark_ff::{FftField, PrimeField};
use ark_poly::univariate::{DensePolynomial};
pub use crate::plonk::circuit::compiled_circuit::CompiledCircuit;
pub use preprocess::*;

#[derive(Clone, Debug)]
pub struct PublicWitness<F: FftField + PrimeField> {
    pub pi_combined: DensePolynomial<F>,
    pub pi_vector: Vec<F>,
    pub output_vector: Vec<F>,
}

pub struct PublicInput<F: FftField + PrimeField> {
    pub pi: DensePolynomial<F>,
    pub pi_vector: Vec<F>,
}

#[derive(Debug, Clone)]
pub struct GateSolution<F: FftField + PrimeField> {
    left: F,
    right: F,
    out: F,
}

#[cfg(test)]
mod tests {
    use ark_poly::Polynomial;
    use ark_std::iterable::Iterable;
    use ark_test_curves::bls12_381::Fr;
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::plonk::domain::PlonkDomain;
    use crate::plonk::test_utils::get_test_circuit;

    #[test]
    pub fn test_s_id() {
        let domain = generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);
        let k1_coset = domain.iter().map(|e| domain.k1() * *e).collect::<Vec<_>>();
        let k2_coset = domain.iter().map(|e| domain.k2() * *e).collect::<Vec<_>>();
        let test_circuit = get_test_circuit(&domain);

        for e in &domain {
            assert_eq!(*e, test_circuit.sid_1.evaluate(e));
        }

        for (ke, e) in k1_coset.iter().zip(domain.iter()) {
            assert_eq!(*ke, test_circuit.sid_2.evaluate(e));
        }

        for (ke, e) in k2_coset.iter().zip(domain.iter()) {
            assert_eq!(*ke, test_circuit.sid_3.evaluate(e));
        }
    }

    #[test]
    pub fn test_sigma_poly() {
        let domain = generate_multiplicative_subgroup::<{ 1 << 3 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);
        let k1_coset = domain.iter().map(|e| domain.k1() * *e).collect::<Vec<_>>();
        let k2_coset = domain.iter().map(|e| domain.k2() * *e).collect::<Vec<_>>();
        let test_circuit = get_test_circuit(&domain);

        let sigma_1_permutations = [
            domain[3],
            k2_coset[5],
            k1_coset[3],
            domain[0],
            k1_coset[4],
            k2_coset[4],
            domain[6],
            domain[7]
        ];
        let sigma_2_permutations = [
            k1_coset[1],
            k1_coset[2],
            k2_coset[0],
            k1_coset[5],
            k2_coset[3],
            domain[2],
            k1_coset[6],
            k1_coset[7],
        ];
        let sigma_3_permutations = [
            k2_coset[1],
            k2_coset[2],
            k1_coset[0],
            domain[4],
            domain[5],
            domain[1],
            k2_coset[6],
            k2_coset[7],
        ];

        for (x, y) in domain.iter().zip(sigma_1_permutations.iter()) {
            assert_eq!(test_circuit.s_sigma_1.evaluate(x), y);
        }

        for (x, y) in domain.iter().zip(sigma_2_permutations.iter()) {
            assert_eq!(test_circuit.s_sigma_2.evaluate(x), y);
        }

        for (x, y) in domain.iter().zip(sigma_3_permutations.iter()) {
            assert_eq!(test_circuit.s_sigma_3.evaluate(x), y);
        }
    }
}