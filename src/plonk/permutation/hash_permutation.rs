use ark_ff::{FftField, PrimeField};
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;
use crate::poly_utils::const_poly;

#[derive(Clone)]
pub struct HashPermutationTerm<'a, F: PrimeField + FftField> {
    f: &'a DensePolynomial<F>,
    shift_poly: &'a DensePolynomial<F>,
    beta: F,
    gamma: F,
}

impl<'a, F: PrimeField + FftField> HashPermutationTerm<'a, F> {
    pub fn new(f: &'a DensePolynomial<F>, shift_poly: &'a DensePolynomial<F>, beta: F, gamma: F) -> Self {
        Self {
            f,
            shift_poly,
            beta,
            gamma,
        }
    }
    pub fn as_poly(&self) -> DensePolynomial<F> {
        self.f + self.shift_poly * self.beta + const_poly(self.gamma)
    }

    pub fn evaluate(&self, point: &F) -> F {
        self.f.evaluate(point) + self.shift_poly.evaluate(point) * self.beta + self.gamma
    }

    pub fn linearize_on_witness(&self, point: &F) -> DensePolynomial<F> {
        const_poly(self.f.evaluate(point)) + self.shift_poly * self.beta + const_poly(self.gamma)
    }
}
