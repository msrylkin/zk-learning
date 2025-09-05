use ark_ff::{FftField, Field, PrimeField};
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;
use crate::plonk::circuit::{CompiledCircuit, Solution};
use crate::poly_utils::interpolate_univariate;

pub struct PermutationArgument<'a, F: PrimeField + FftField> {
    domain: &'a [F],
    beta: &'a F,
    gamma: &'a F,
    solution: &'a Solution<F>,
}

impl<'a, F: PrimeField + FftField> PermutationArgument<'a, F> {
    pub fn new(domain: &'a [F], beta: &'a F, gamma: &'a F, solution: &'a Solution<F>) -> Self {
        Self {
            domain,
            beta,
            gamma,
            solution
        }
    }

    fn hash_permutation(
        &self,
        f: &DensePolynomial<F>,
        perm_poly: &DensePolynomial<F>,
        point: &F,
    ) -> F {
        f.evaluate(&point) + *self.beta * perm_poly.evaluate(&point) + *self.gamma
    }

    pub fn numerator_acc(
        &self,
        max_i: usize,
    ) -> F {
        self.permutation_product_acc(
            &self.solution.sid_1,
            &self.solution.sid_2,
            &self.solution.sid_3,
            max_i,
        )
    }

    pub fn denominator_acc(
        &self,
        max_i: usize,
    ) -> F {
        self.permutation_product_acc(
            &self.solution.s_sigma_1,
            &self.solution.s_sigma_2,
            &self.solution.s_sigma_3,
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

    pub fn denominator_poly(
        &self,
    ) -> DensePolynomial<F> {
        let values = self.domain.iter().map(|w| {
            self.permutation_product(
                &self.solution.s_sigma_1,
                &self.solution.s_sigma_2,
                &self.solution.s_sigma_3,
                &w
            )
        }).collect::<Vec<_>>();

        interpolate_univariate(self.domain, &values)
    }

    pub fn numerator_poly(&self) -> DensePolynomial<F> {
        let values = self.domain.iter().map(|w| {
            self.permutation_product(
                &self.solution.sid_1,
                &self.solution.sid_2,
                &self.solution.sid_3,
                &w
            )
        }).collect::<Vec<_>>();

        interpolate_univariate(self.domain, &values)
    }

    fn permutation_product(
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

        interpolate_univariate(self.domain, &values)
    }
}

#[cfg(test)]
mod tests {
    use ark_std::One;
    use ark_test_curves::bls12_381::Fr;
    use crate::plonk::circuit::get_test_circuit;
    use crate::plonk::permutation::PermutationArgument;
    use crate::plonk::prover::{pick_coset_shifters};
    use crate::poly_utils::generate_multiplicative_subgroup;

    #[test]
    fn test_permutation_poly_acc() {
        let test_circuit = get_test_circuit();
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 3 }, Fr>();
        let (k1, k2) = pick_coset_shifters(&domain);

        let beta = Fr::from(43);
        let gamma = Fr::from(35);

        let solution = test_circuit.get_solution(&domain, k1, k2);
        let permutation = PermutationArgument::new(&domain, &beta, &gamma, &solution);

        let num = permutation.numerator_acc(
            domain.len(),
        );
        let denom = permutation.denominator_acc(
            domain.len(),
        );

        assert_eq!(num / denom, Fr::one());
    }
}