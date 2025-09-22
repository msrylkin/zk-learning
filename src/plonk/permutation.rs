use ark_ff::{FftField, Field, PrimeField};
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;
use crate::evaluation_domain::MultiplicativeSubgroup;
use crate::plonk::circuit::{CompiledCircuit};
use crate::plonk::circuit::solution::Solution;
use crate::poly_utils::{const_poly};

#[derive(Clone)]
pub struct PermutationArgument<'a, F: PrimeField + FftField> {
    domain: &'a MultiplicativeSubgroup<F>,
    beta: &'a F,
    gamma: &'a F,
    solution: &'a Solution<F>,
    circuit: &'a CompiledCircuit<'a, F>,
}

impl<'a, F: PrimeField + FftField> PermutationArgument<'a, F> {
    pub fn new(
        domain: &'a MultiplicativeSubgroup<F>,
        beta: &'a F,
        gamma: &'a F,
        circuit: &'a CompiledCircuit<F>,
        solution: &'a Solution<F>,
    ) -> Self {
        Self {
            domain,
            beta,
            gamma,
            solution,
            circuit,
        }
    }

    pub fn hash_permutation(
        &self,
        f: &DensePolynomial<F>,
        perm_poly: &DensePolynomial<F>,
        point: &F,
    ) -> F {
        f.evaluate(&point) + *self.beta * perm_poly.evaluate(&point) + *self.gamma
    }

    pub fn hash_permutation_poly_old(
        &self,
        f: &DensePolynomial<F>,
        perm_poly: &DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        let mut values = vec![];

        for w in self.domain {
            values.push(self.hash_permutation(f, perm_poly, w));
        }

        self.domain.interpolate_univariate(&values)
    }

    pub fn hash_permutation_poly(
        &self,
        f: &DensePolynomial<F>,
        perm_poly: &DensePolynomial<F>,
    ) -> DensePolynomial<F> {
        f + perm_poly * const_poly(*self.beta) + const_poly(*self.gamma)
    }

    fn numerator_acc(
        &self,
        max_i: usize,
    ) -> F {
        self.permutation_product_acc(
            &self.circuit.sid_1,
            &self.circuit.sid_2,
            &self.circuit.sid_3,
            max_i,
        )
    }

    pub fn denominator_acc(
        &self,
        max_i: usize,
    ) -> F {
        self.permutation_product_acc(
            &self.circuit.s_sigma_1,
            &self.circuit.s_sigma_2,
            &self.circuit.s_sigma_3,
            max_i,
        )
    }

    fn permutation_product_acc(
        &self,
        perm_poly_1: &DensePolynomial<F>,
        perm_poly_2: &DensePolynomial<F>,
        perm_poly_3: &DensePolynomial<F>,
        max_i: usize,
    ) -> F {
        self.domain
            .iter()
            .take(max_i)
            .map(|w| self.permutation_product(perm_poly_1, perm_poly_2, perm_poly_3, &w))
            .product()
    }

    pub fn denominator_poly_old(
        &self,
    ) -> DensePolynomial<F> {
        let values = self.domain.iter().map(|w| {
            self.permutation_product(
                &self.circuit.s_sigma_1,
                &self.circuit.s_sigma_2,
                &self.circuit.s_sigma_3,
                &w
            )
        }).collect::<Vec<_>>();

        self.domain.interpolate_univariate(&values)
    }

    pub fn denominator_poly(&self) -> DensePolynomial<F> {
        self.hash_permutation_poly(&self.solution.a, &self.circuit.s_sigma_1)
            * self.hash_permutation_poly(&self.solution.b, &self.circuit.s_sigma_2)
            * self.hash_permutation_poly(&self.solution.c, &self.circuit.s_sigma_3)
    }

    pub fn numerator_poly_old(&self) -> DensePolynomial<F> {
        let values = self.domain.iter().map(|w| {
            self.permutation_product(
                &self.circuit.sid_1,
                &self.circuit.sid_2,
                &self.circuit.sid_3,
                &w
            )
        }).collect::<Vec<_>>();

        self.domain.interpolate_univariate(&values)
    }

    pub fn numerator_poly(&self) -> DensePolynomial<F> {
        self.hash_permutation_poly(&self.solution.a, &self.circuit.sid_1)
            * self.hash_permutation_poly(&self.solution.b, &self.circuit.sid_2)
            * self.hash_permutation_poly(&self.solution.c, &self.circuit.sid_3)
    }

    pub fn permutation_product(
        &self,
        perm_poly_1: &DensePolynomial<F>,
        perm_poly_2: &DensePolynomial<F>,
        perm_poly_3: &DensePolynomial<F>,
        point: &F,
    ) -> F {
        self.hash_permutation(&self.solution.a, &perm_poly_1, point)
            * self.hash_permutation(&self.solution.b, &perm_poly_2, point)
            * self.hash_permutation(&self.solution.c, &perm_poly_3, point)
    }

    pub fn z_poly(&self) -> DensePolynomial<F> {
        let values = (0..self.domain.len())
            .map(|i| {
                let num = self.numerator_acc(i);
                let denom = self.denominator_acc(i);

                num / denom
            })
            .collect::<Vec<_>>();

        self.domain.interpolate_univariate(&values)
    }

    pub fn beta(&self) -> F {
        self.beta.clone()
    }

    pub fn gamma(&self) -> F {
        self.gamma.clone()
    }
}

#[cfg(test)]
mod tests {
    use ark_poly::Polynomial;
    use ark_poly::univariate::DensePolynomial;
    use ark_std::One;
    use ark_test_curves::bls12_381::Fr;
    use crate::plonk::circuit::{get_test_circuit, get_test_solution};
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::plonk::domain::PlonkDomain;
    use crate::plonk::permutation::PermutationArgument;
    // use crate::plonk::prover::{pick_coset_shifters};

    #[test]
    fn test_permutation_poly_acc() {
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 3 }, Fr>();
        let domain = PlonkDomain::new(&domain);
        let test_circuit = get_test_circuit(&domain);
        // let (k1, k2) = pick_coset_shifters(&domain);

        let beta = Fr::from(43);
        let gamma = Fr::from(35);

        let solution= get_test_solution(&domain);
        let permutation = PermutationArgument::new(&domain, &beta, &gamma, &test_circuit, &solution);

        let num = permutation.numerator_acc(
            domain.len(),
        );
        let denom = permutation.denominator_acc(
            domain.len(),
        );

        assert_eq!(num / denom, Fr::one());
    }

    #[test]
    fn test_hash_permutation_poly() {
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 3 }, Fr>();
        let domain = PlonkDomain::new(&domain);

        let test_circuit = get_test_circuit(&domain);
        let Zh = domain.get_vanishing_polynomial();
        let Zh = DensePolynomial::from(Zh);
        // let (k1, k2) = pick_coset_shifters(&domain);

        let beta = Fr::from(43);
        let gamma = Fr::from(35);

        let solution = get_test_solution(&domain);
        let permutation = PermutationArgument::new(&domain, &beta, &gamma, &test_circuit, &solution);

        let num_poly = permutation.numerator_poly();
        let denom_poly = permutation.denominator_poly();

        let custom_num_poly = permutation.hash_permutation_poly(&solution.a, &test_circuit.sid_1)
            * permutation.hash_permutation_poly(&solution.b, &test_circuit.sid_2)
            * permutation.hash_permutation_poly(&solution.c, &test_circuit.sid_3);
        
        let mut reduced_num_values = vec![];
        for w in &domain {
            reduced_num_values.push(custom_num_poly.evaluate(&w));
        }
        let reduced_num_poly = domain.interpolate_univariate(&reduced_num_values);
        
        let custom_denom_poly = permutation.hash_permutation_poly(&solution.a, &test_circuit.s_sigma_1)
            * permutation.hash_permutation_poly(&solution.b, &test_circuit.s_sigma_2)
            * permutation.hash_permutation_poly(&solution.c, &test_circuit.s_sigma_3);

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
}