use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

pub trait OracleEvaluation<F: Field> {
    fn final_eval(&self, r: &[F]) -> F;
}

pub trait SumCheckPoly<F: Field> {
    fn get_evaluations(&self) -> Vec<F>;

    fn num_vars(&self) -> usize;

    fn get_partial_sum_poly(&self) -> DensePolynomial<F>;

    fn fix_variable(&self, e: F) -> Self;

    fn evaluate(&self, point: &[F]) -> F;
}
