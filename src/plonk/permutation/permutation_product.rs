use ark_ff::{FftField, PrimeField};
use ark_poly::univariate::DensePolynomial;
use crate::plonk::permutation::hash_permutation::HashPermutationTerm;

pub struct PermutationProduct<'a, F: PrimeField + FftField> {
    left: HashPermutationTerm<'a, F>,
    right: HashPermutationTerm<'a, F>,
    out: HashPermutationTerm<'a, F>,
}

impl<'a, F: PrimeField + FftField> PermutationProduct<'a, F> {
    pub fn new(left: HashPermutationTerm<'a, F>, right: HashPermutationTerm<'a, F>, out: HashPermutationTerm<'a, F>) -> Self {
        Self {
            left,
            right,
            out,
        }
    }
    pub fn combined(&self) -> DensePolynomial<F> {
        self.left.as_poly() * self.right.as_poly() * self.out.as_poly()
    }

    pub fn evaluate(&self, point: &F) -> F {
        self.left.evaluate(point) * self.right.evaluate(point) * self.out.evaluate(point)
    }

    pub fn evaluate_left(&self, point: &F) -> F {
        self.left.evaluate(point)
    }

    pub fn evaluate_right(&self, point: &F) -> F {
        self.right.evaluate(point)
    }

    pub fn evaluate_out(&self, point: &F) -> F {
        self.out.evaluate(point)
    }

    pub fn linearize_left_on_witness(&self, point: &F) -> DensePolynomial<F> {
        self.left.linearize_on_witness(point)
    }

    pub fn linearize_right_on_witness(&self, point: &F) -> DensePolynomial<F> {
        self.right.linearize_on_witness(point)
    }
    pub fn linearize_out_on_witness(&self, point: &F) -> DensePolynomial<F> {
        self.out.linearize_on_witness(point)
    }
}
