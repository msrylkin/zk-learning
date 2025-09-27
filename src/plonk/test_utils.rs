#[cfg(test)]
pub mod test_utils {
    use ark_ec::pairing::Pairing;
    use ark_ff::{FftField, Field, PrimeField};
    use ark_poly::univariate::DensePolynomial;
    use crate::kzg::{setup, KZG};
    use crate::plonk::circuit::circuit_description::{CircuitDescription};
    use crate::plonk::circuit::CompiledCircuit;
    use crate::plonk::circuit::solution::Solution;
    use crate::plonk::domain::PlonkDomain;
    use crate::poly_utils::const_poly;

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
        let min1 = circuit_builder.constant_var(F::from(-1));
        circuit_builder.make_public(a);
        let c = circuit_builder.add_variable();

        let mul_result_1 = circuit_builder.multiplication_gate(a, b);
        let mul_result_2 = circuit_builder.multiplication_gate(mul_result_1, mul_result_1);
        let add_result = circuit_builder.addition_gate(mul_result_2, c);
        let mul_result_3 = circuit_builder.multiplication_gate(add_result, min1);
        let mul_result_4 = circuit_builder.multiplication_gate(mul_result_2, mul_result_3);

        circuit_builder.make_output(mul_result_2);
        circuit_builder.make_output(mul_result_4);

        circuit_builder
    }

    #[cfg(test)]
    pub fn get_test_circuit<F: FftField + PrimeField>(domain: &PlonkDomain<F>) -> CompiledCircuit<F> {
        build_test_circuit().compile(domain)
    }

    #[cfg(test)]
    pub fn get_test_solution<F: FftField + PrimeField>(domain: &PlonkDomain<F>) -> Solution<F> {
        build_test_circuit().solve(&[F::from(9)], &[F::from(-544000)], domain)
    }

    #[cfg(test)]
    pub fn get_test_kzg<P: Pairing>(n: usize) -> KZG<P> {
        let tau = P::ScalarField::from(777);
        let config = setup::<P>(n * 2, tau);
        let kzg = KZG::new(config);

        kzg
    }
}


