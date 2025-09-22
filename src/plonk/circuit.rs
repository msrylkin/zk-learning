pub mod circuit_description;
pub mod solution;
pub mod compiled_circuit;

use ark_ff::{FftField, PrimeField};
use ark_poly::univariate::{DensePolynomial};
use crate::evaluation_domain::MultiplicativeSubgroup;
use crate::plonk::circuit::circuit_description::CircuitDescription;
pub(crate) use crate::plonk::circuit::compiled_circuit::CompiledCircuit;
use crate::plonk::circuit::solution::Solution;
use crate::plonk::domain::PlonkDomain;

enum CircuitOperation {
    Mul,
    Add,
}

#[derive(Debug, Clone)]
enum Operation {
    Multiplication,
    Addition,
}

fn build_test_circuit<F: FftField + PrimeField>() -> CircuitDescription<F> {
    let mut circuit_builder = CircuitDescription::new();

    let a = circuit_builder.add_variable();
    let b = circuit_builder.constant_var(F::from(82));
    circuit_builder.make_public(a);

    let mul_result_1 = circuit_builder.multiplication_gate(a, b);
    let mul_result_2 = circuit_builder.multiplication_gate(mul_result_1, mul_result_1);
    let add_result = circuit_builder.addition_gate(mul_result_2, b);

    circuit_builder.make_public(add_result);

    circuit_builder
}

#[derive(Clone)]
pub struct PublicInput<F: FftField + PrimeField> {
    pub pi: DensePolynomial<F>,
    pub pi_vector: Vec<F>,
}

#[derive(Debug)]
pub struct GateSolution<F: FftField + PrimeField> {
    left: F,
    right: F,
    out: F,
}

pub fn get_test_circuit<'a, F: FftField + PrimeField>(domain: &'a PlonkDomain<F>) -> CompiledCircuit<'a, F> {
    let circuit_description = build_test_circuit();
    circuit_description.compile(&domain)
}

pub fn get_test_solution<F: FftField + PrimeField>(domain: &PlonkDomain<F>) -> Solution<F> {
    let compiled_circuit = get_test_circuit(domain);
    compiled_circuit.solve(&[F::from(9), F::from(544726)], &[])
}

#[cfg(test)]
mod tests {
    use ark_poly::Polynomial;
    use ark_std::iterable::Iterable;
    use ark_std::Zero;
    use ark_test_curves::bls12_381::Fr;
    use crate::plonk::circuit::{get_test_circuit};
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::plonk::domain::PlonkDomain;

    #[test]
    pub fn test_s_id() {
        let domain = generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>();
        let domain = PlonkDomain::new(&domain);
        // let (k1, k2) = pick_coset_shifters(&domain);
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
        let domain = PlonkDomain::new(&domain);
        // let (k1, k2) = pick_coset_shifters(&domain);
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