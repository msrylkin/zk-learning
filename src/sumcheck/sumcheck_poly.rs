use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

/// Trait representing the final evaluation of an oracle
/// in the sumcheck protocol.
///
/// The oracle produces a single field element after all
/// randomness values from the protocol have been applied.
pub trait OracleEvaluation<F: Field> {
    /// Computes the final evaluation given all randomness `r` used in the protocol.
    fn final_eval(&self, r: &[F]) -> F;
}

/// Trait for multivariate polynomials used in the sumcheck protocol.
///
/// Provides an interface for accessing evaluations, partially summing
/// over variables, fixing variables, and evaluating the polynomial at a point.
pub trait SumCheckPoly<F: Field> {
    // TODO: remove?
    /// Returns the polynomial in its full evaluation form.
    fn get_evaluations(&self) -> Vec<F>;

    /// Returns the number of variables in the polynomial.
    fn num_vars(&self) -> usize;

    /// Computes the partial sum polynomial, summing over all variables
    /// except for the first one.
    fn get_partial_sum_poly(&self) -> DensePolynomial<F>;

    /// Returns a new polynomial with the first variable fixed to the value `e`.
    fn fix_variable(&self, e: F) -> Self;

    /// Evaluates the polynomial at the given `point`.
    fn evaluate(&self, point: &[F]) -> F;
}
