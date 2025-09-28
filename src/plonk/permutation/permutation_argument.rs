use ark_ff::{FftField, PrimeField};
use ark_poly::univariate::DensePolynomial;
use crate::evaluation_domain::MultiplicativeSubgroup;
use crate::plonk::circuit::CompiledCircuit;
use crate::plonk::circuit::PlonkSolution;
use crate::plonk::permutation::hash_permutation::HashPermutationTerm;
use crate::plonk::permutation::permutation_product::PermutationProduct;

#[derive(Clone)]
pub struct PermutationArgument<'a, F: PrimeField + FftField> {
    domain: &'a MultiplicativeSubgroup<F>,
    beta: F,
    gamma: F,
    solution: &'a PlonkSolution<F>,
    circuit: &'a CompiledCircuit<F>,
}

impl<'a, F: PrimeField + FftField> PermutationArgument<'a, F> {
    pub fn new(
        domain: &'a MultiplicativeSubgroup<F>,
        beta: F,
        gamma: F,
        circuit: &'a CompiledCircuit<F>,
        solution: &'a PlonkSolution<F>,
    ) -> Self {
        Self {
            domain,
            beta,
            gamma,
            solution,
            circuit,
        }
    }

    pub fn numerator(&self) -> PermutationProduct<F> {
        PermutationProduct::new(
            HashPermutationTerm::new(&self.solution.a, &self.circuit.sid_1, self.beta, self.gamma),
            HashPermutationTerm::new(&self.solution.b, &self.circuit.sid_2, self.beta, self.gamma),
            HashPermutationTerm::new(&self.solution.c, &self.circuit.sid_3, self.beta, self.gamma),
        )
    }

    pub fn denominator(&self) -> PermutationProduct<F> {
        PermutationProduct::new(
            HashPermutationTerm::new(&self.solution.a, &self.circuit.s_sigma_1, self.beta, self.gamma),
            HashPermutationTerm::new(&self.solution.b, &self.circuit.s_sigma_2, self.beta, self.gamma),
            HashPermutationTerm::new(&self.solution.c, &self.circuit.s_sigma_3, self.beta, self.gamma),
        )
    }

    pub fn z_poly(&self) -> DensePolynomial<F> {
        let numerator = self.numerator();
        let denominator = self.denominator();
        let mut prev_value = F::one();
        let mut values = vec![prev_value];

        for w in self.domain.iter().take(self.domain.len() - 1) {
            prev_value *= numerator.evaluate(w) / denominator.evaluate(w);
            values.push(prev_value);
        }

        self.domain.interpolate_univariate(&values)
    }
}

#[cfg(test)]
mod tests {
    use ark_poly::Polynomial;
    use ark_std::One;
    use ark_test_curves::bls12_381::Fr;
    use crate::evaluation_domain::{generate_multiplicative_subgroup};
    use crate::plonk::domain::PlonkDomain;
    use crate::plonk::permutation::permutation_argument::PermutationArgument;
    use crate::plonk::test_utils::{get_test_circuit, get_test_solution, hash_permutation_poly};

    #[test]
    fn test_hash_permutation_poly() {
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 4 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);

        let test_circuit = get_test_circuit(&domain);

        let beta = Fr::from(43);
        let gamma = Fr::from(35);

        let solution = get_test_solution(&domain);
        let permutation = PermutationArgument::new(&domain, beta, gamma, &test_circuit, &solution);

        let num_poly = permutation.numerator().combined();
        let denom_poly = permutation.denominator().combined();

        let custom_num_poly = hash_permutation_poly(&solution.a, &test_circuit.sid_1, beta, gamma)
            * hash_permutation_poly(&solution.b, &test_circuit.sid_2, beta, gamma)
            * hash_permutation_poly(&solution.c, &test_circuit.sid_3, beta, gamma);

        let mut reduced_num_values = vec![];
        for w in &domain {
            reduced_num_values.push(custom_num_poly.evaluate(&w));
        }
        let reduced_num_poly = domain.interpolate_univariate(&reduced_num_values);

        let custom_denom_poly = hash_permutation_poly(&solution.a, &test_circuit.s_sigma_1, beta, gamma)
            * hash_permutation_poly(&solution.b, &test_circuit.s_sigma_2, beta, gamma)
            * hash_permutation_poly(&solution.c, &test_circuit.s_sigma_3, beta, gamma);

        let mut reduced_denom_values = vec![];
        for w in &domain {
            reduced_denom_values.push(custom_denom_poly.evaluate(&w));
        }
        let reduced_denom_poly = domain.interpolate_univariate(&reduced_denom_values);

        for w in &domain {
            assert_eq!(num_poly.evaluate(&w), reduced_num_poly.evaluate(&w));
            assert_eq!(num_poly.evaluate(&w), custom_num_poly.evaluate(&w));
            assert_eq!(denom_poly.evaluate(&w), reduced_denom_poly.evaluate(&w));
            assert_eq!(denom_poly.evaluate(&w), custom_denom_poly.evaluate(&w));
        }

        assert_eq!(
            num_poly,
            custom_num_poly,
        );
        assert_ne!(
            num_poly,
            reduced_num_poly,
        );

        assert_eq!(
            denom_poly,
            custom_denom_poly,
        );
        assert_ne!(
            denom_poly,
            reduced_denom_poly,
        );
    }

    #[test]
    fn test_z_poly() {
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 4 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);

        let test_circuit = get_test_circuit(&domain);
        let solution = get_test_solution(&domain);

        let beta = Fr::from(43);
        let gamma = Fr::from(35);

        let permutation = PermutationArgument::new(&domain, beta, gamma, &test_circuit, &solution);

        let z_poly = permutation.z_poly();

        assert_eq!(z_poly.evaluate(domain.last().unwrap()), Fr::one());
        assert_eq!(z_poly.evaluate(&Fr::one()), Fr::one());

        let omega= domain.generator();
        let num_poly = permutation.numerator().combined();
        let denom_poly = permutation.denominator().combined();

        for x in &domain {
            assert_eq!(num_poly.evaluate(&x), permutation.numerator().evaluate(x));
            assert_eq!(denom_poly.evaluate(&x), permutation.denominator().evaluate(x));
            assert_eq!(
                z_poly.evaluate(&(omega * x)) * denom_poly.evaluate(x),
                z_poly.evaluate(x) * num_poly.evaluate(x),
            );
        }
    }
}
