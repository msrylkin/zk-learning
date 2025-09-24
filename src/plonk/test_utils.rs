use ark_ec::pairing::Pairing;
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::univariate::DensePolynomial;
use ark_test_curves::bls12_381::{Bls12_381, Fr};
use crate::kzg::{setup, KZG};
use crate::plonk::circuit::circuit_description::CircuitDescription;
use crate::plonk::circuit::CompiledCircuit;
use crate::plonk::circuit::solution::Solution;
use crate::plonk::domain::PlonkDomain;
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

#[cfg(test)]
pub fn build_test_circuit<F: FftField + PrimeField>() -> CircuitDescription<F> {
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

#[cfg(test)]
pub fn get_test_circuit<F: FftField + PrimeField>(domain: &PlonkDomain<F>) -> CompiledCircuit<F> {
    let circuit_description = build_test_circuit();
    circuit_description.compile(&domain)
}

#[cfg(test)]
pub fn get_test_solution<F: FftField + PrimeField>(domain: &PlonkDomain<F>) -> Solution<F> {
    let compiled_circuit = get_test_circuit(domain);
    compiled_circuit.solve(&[F::from(9), F::from(544726)], &[], &domain)
}