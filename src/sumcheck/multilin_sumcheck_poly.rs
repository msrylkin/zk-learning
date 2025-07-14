use std::ops::Deref;
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use crate::poly_utils::interpolate_univariate_on_evals;
use crate::sumcheck::OracleEvaluation;
use crate::sumcheck::sumcheck_poly::SumCheckPoly;

// 000 |
// 001 --> 00[1 + 0]
// 010 |
// 011 --> 01[1 + 0]
fn sum_over_last_variable<F: Field>(poly: &DenseMultilinSumcheckPoly<F>) -> DenseMultilinSumcheckPoly<F> {
    let sums = poly
        .evaluations
        .iter()
        .enumerate()
        .fold(vec![F::zero(); poly.evaluations.len() / 2], |mut sums, (i, e)| {
            if i % 2 == 0 {
                sums[(i / 4) * 2] += e;
            } else {
                sums[(i / 4) * 2 + 1] += e;
            }
            sums
        });

    let res = DenseMultilinSumcheckPoly::from_evaluations_vec(poly.num_vars - 1, sums);

    res
}

#[derive(Clone)]
pub struct DenseMultilinSumcheckPoly<F: Field>(DenseMultilinearExtension<F>);

impl<F: Field> DenseMultilinSumcheckPoly<F> {
    pub fn new(poly: DenseMultilinearExtension<F>) -> Self {
        Self(poly)
    }
    pub fn from_evaluations_vec(num_vars: usize, evals: Vec<F>) -> Self {
        Self(DenseMultilinearExtension::from_evaluations_vec(num_vars, evals))
    }
}

impl<F: Field> Deref for DenseMultilinSumcheckPoly<F> {
    type Target = DenseMultilinearExtension<F>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<F: Field> SumCheckPoly<F> for DenseMultilinSumcheckPoly<F> {
    fn get_evaluations(&self) -> Vec<F> {
        self.0.to_evaluations()
    }

    fn num_vars(&self) -> usize {
        self.0.num_vars()
    }

    fn get_partial_sum_poly(&self) -> DensePolynomial<F> {
        let sum_vars = self.num_vars() - 1;
        let mut current_poly = self.clone();

        for _ in 0..sum_vars {
            let summed_poly = sum_over_last_variable(&current_poly);
            current_poly = summed_poly;
        }

        interpolate_univariate_on_evals(&current_poly.get_evaluations().try_into().unwrap())
    }

    fn fix_variables(&self, partial_point: &[F]) -> Self {
        DenseMultilinSumcheckPoly(self.0.fix_variables(partial_point))
    }

    fn evaluate(&self, point: &[F]) -> F {
        self.0.evaluate(&point.to_vec())
    }
}

pub struct DenseMultilinFinalEvaluationOracle<F: Field>(DenseMultilinSumcheckPoly<F>);

impl<F: Field> DenseMultilinFinalEvaluationOracle<F> {
    pub fn new(poly: DenseMultilinSumcheckPoly<F>) -> Self {
        Self(poly)
    }
}

impl<F: Field> OracleEvaluation<F> for DenseMultilinFinalEvaluationOracle<F> {
    fn final_eval(&self, r: &[F]) -> F {
        self.0.evaluate(&r.to_vec())
    }
}
