mod circuit_description;
mod solution;
mod compiled_circuit;
mod preprocess;

use ark_ff::{FftField, PrimeField};
use ark_poly::univariate::{DensePolynomial};
pub use compiled_circuit::CompiledCircuit;
pub use circuit_description::CircuitDescription;
pub use solution::PlonkSolution;
pub use preprocess::{preprocess_circuit, PreprocessedCircuit};

/// Representation of the public witness for a circuit.
///
/// Contains:
/// - `inputs_vector`: the public input values,
/// - `output_vector`: the circuit output values,
/// - `pi_combined`: a polynomial obtained by interpolating inputs and outputs together.
#[derive(Clone, Debug)]
pub struct PublicWitness<F: FftField + PrimeField> {
    /// Polynomial interpolating both inputs and outputs.
    pub pi_combined: DensePolynomial<F>,
    /// Vector of circuit inputs.
    pub inputs_vector: Vec<F>,
    /// Vector of circuit outputs.
    pub output_vector: Vec<F>,
}

/// Witness values for a single gate, including its left input, right input, and output.
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
    use crate::plonk::test_utils::test_utils::get_test_circuit;

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
        let domain = generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);
        let k1_coset = domain.iter().map(|e| domain.k1() * *e).collect::<Vec<_>>();
        let k2_coset = domain.iter().map(|e| domain.k2() * *e).collect::<Vec<_>>();
        let test_circuit = get_test_circuit(&domain);

        let sigma_1_permutations = [
            domain[5],
            domain[7],
            k2_coset[9],
            k1_coset[5],
            k1_coset[8],
            domain[0],
            k1_coset[6],
            domain[9],
            k2_coset[7],
            k2_coset[6],
            domain[10],
            domain[11],
            domain[12],
            domain[13],
            domain[14],
            domain[15],
        ];
        let sigma_2_permutations = [
            k1_coset[1],
            k1_coset[2],
            k1_coset[3],
            k1_coset[4],
            k2_coset[0],
            domain[3],
            k2_coset[5],
            k1_coset[7],
            domain[4],
            k2_coset[8],
            k1_coset[10],
            k1_coset[11],
            k1_coset[12],
            k1_coset[13],
            k1_coset[14],
            k1_coset[15],
        ];
        let sigma_3_permutations = [
            k2_coset[1],
            k2_coset[2],
            k2_coset[3],
            k2_coset[4],
            k1_coset[0],
            domain[6],
            domain[1],
            domain[8],
            k1_coset[9],
            domain[2],
            k2_coset[10],
            k2_coset[11],
            k2_coset[12],
            k2_coset[13],
            k2_coset[14],
            k2_coset[15],
        ];

        assert_eq!(sigma_1_permutations.len(), domain.len());
        assert_eq!(sigma_2_permutations.len(), domain.len());

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