use ark_ff::{FftField, PrimeField};
use ark_poly::univariate::DensePolynomial;
use crate::evaluation_domain::MultiplicativeSubgroup;
use crate::plonk::circuit::{GateSolution, PublicWitness};

/// Represents the complete solution for a Plonk circuit.
///
/// Contains all gates including constraint gates and synthetic gates for inputs,
/// as well as the column polynomials `a`, `b`, `c` and a public witness representation.
pub struct PlonkSolution<F: FftField + PrimeField> {
    /// Solved gate values for left, right, and output columns.
    pub solution_gates: Vec<GateSolution<F>>,
    /// `a` polynomial representing the left column of gates.
    pub a: DensePolynomial<F>,
    /// `b` polynomial representing the right column of gates.
    pub b: DensePolynomial<F>,
    /// `c` polynomial representing the output column of gates.
    pub c: DensePolynomial<F>,
    /// Representation of public inputs, outputs, and combined `pi` polynomial.
    pub public_witness: PublicWitness<F>,
}

impl<F: FftField + PrimeField> PlonkSolution<F> {
    /// Constructs a new `PlonkSolution` from solved gates and a public witness.
    ///
    /// Interpolates the gate values into univariate polynomials over the given `domain`.
    pub fn new(
        solution_gates: Vec<GateSolution<F>>,
        domain: &MultiplicativeSubgroup<F>,
        public_witness: PublicWitness<F>,
    ) -> Self {
        let (a, b, c) = Self::get_abc_polys(&solution_gates, domain);

        Self {
            a,
            b,
            c,
            solution_gates,
            public_witness,
        }
    }

    /// Returns vectors of gate values for the `a`, `b`, and `c` columns,
    /// padding with zeros to match the domain size.
    fn get_abc_vectors(solution_gates: &[GateSolution<F>], domain: &MultiplicativeSubgroup<F>) -> (Vec<F>, Vec<F>, Vec<F>) {
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

    /// Interpolates vectors of `a`, `b`, and `c` gate values into univariate polynomials.
    fn get_abc_polys(solution_gates: &[GateSolution<F>], domain: &MultiplicativeSubgroup<F>) -> (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>) {
        let (a, b, c) = Self::get_abc_vectors(solution_gates, domain);

        let a = domain.interpolate_univariate(&a);
        let b = domain.interpolate_univariate(&b);
        let c = domain.interpolate_univariate(&c);

        (a, b, c)
    }
}
