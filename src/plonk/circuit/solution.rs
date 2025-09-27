use ark_ff::{FftField, PrimeField};
use ark_poly::univariate::DensePolynomial;
use crate::evaluation_domain::MultiplicativeSubgroup;
use crate::plonk::circuit::{GateSolution, PublicWitness};

pub struct Solution<F: FftField + PrimeField> {
    pub solution_gates: Vec<GateSolution<F>>,
    pub a: DensePolynomial<F>,
    pub b: DensePolynomial<F>,
    pub c: DensePolynomial<F>,
    pub public_witness: PublicWitness<F>,
    // pub output: Vec<F>,
    // pub output: Vec<F>,
}

impl<F: FftField + PrimeField> Solution<F> {
    pub fn new(
        solution_gates: Vec<GateSolution<F>>,
        domain: &MultiplicativeSubgroup<F>,
        public_witness: PublicWitness<F>,
        // output: Vec<F>,
        // output: Vec<F>,
    ) -> Self {
        let (a, b, c) = Self::get_abc_polys(&solution_gates, domain);

        Self {
            a,
            b,
            c,
            solution_gates,
            public_witness,
            // output,
        }
    }

    pub fn get_abc_vectors(solution_gates: &[GateSolution<F>], domain: &MultiplicativeSubgroup<F>) -> (Vec<F>, Vec<F>, Vec<F>) {
        let mut a = vec![];
        let mut b = vec![];
        let mut c = vec![];

        for gate in solution_gates {
            a.push(gate.left);
            b.push(gate.right);
            c.push(gate.out);
        }

        for _ in solution_gates.len()..domain.len() {
            a.push(F::zero());
            b.push(F::zero());
            c.push(F::zero());
        }

        (a, b, c)
    }

    pub fn get_abc_polys(solution_gates: &[GateSolution<F>], domain: &MultiplicativeSubgroup<F>) -> (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>) {
        let (a, b, c) = Self::get_abc_vectors(solution_gates, domain);

        let a = domain.interpolate_univariate(&a);
        let b = domain.interpolate_univariate(&b);
        let c = domain.interpolate_univariate(&c);

        (a, b, c)
    }
}
