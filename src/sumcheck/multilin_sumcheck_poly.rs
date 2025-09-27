use std::ops::Deref;
use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use crate::poly_utils::{interpolate_univariate_on_evals};
use crate::sumcheck::OracleEvaluation;
use crate::sumcheck::sumcheck_poly::SumCheckPoly;

#[derive(Clone, Debug)]
pub struct DenseMultilinSumcheckPoly<F: Field>(DenseMultilinearExtension<F>);

impl<F: Field> DenseMultilinSumcheckPoly<F> {
    pub fn new(poly: DenseMultilinearExtension<F>) -> Self {
        Self(poly)
    }
    pub fn from_evaluations_vec(num_vars: usize, evals: Vec<F>) -> Self {
        Self(DenseMultilinearExtension::from_evaluations_vec(num_vars, evals))
    }
    
    pub fn sum_over_last_variable(&self) -> DenseMultilinSumcheckPoly<F> {
        let new_vars_num = self.num_vars() - 1;

        let res = (0..1 << new_vars_num).map(|i| {
            let index = i;
            let partner_index = i | (1 << new_vars_num);

            self.evaluations[index] + self.evaluations[partner_index]
        }).collect::<Vec<_>>();
        
        DenseMultilinSumcheckPoly::from_evaluations_vec(new_vars_num, res)
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
        let (zeroes, ones) = self
            .evaluations
            .iter()
            .enumerate()
            .fold((F::zero(), F::zero()), |(zeroes, ones), (i, e)| {
                if i & 1 == 0 {
                    (zeroes + e, ones)
                } else {
                    (zeroes, ones + e)
                }
            });

        interpolate_univariate_on_evals(&[zeroes, ones])
    }

    fn fix_variable(&self, e: F) -> Self {
        DenseMultilinSumcheckPoly(self.0.fix_variables(&[e]))
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
        self.0.evaluate(r)
    }
}


#[cfg(test)]
mod tests {
    use ark_ff::Field;
    use ark_poly::{MultilinearExtension};
    use ark_test_curves::bls12_381::Fr;
    use crate::poly_utils::{to_f};
    use crate::sumcheck::multilin_sumcheck_poly::{DenseMultilinSumcheckPoly};
    use crate::sumcheck::SumCheckPoly;
    
    fn test_poly<const VARS: usize, F: Field>() -> DenseMultilinSumcheckPoly<F> {
        let values = match VARS {
            1 => vec![10, 20],
            2 => vec![10, 300, 20, 400],
            3 => vec![10, 55, 300, 77, 20, 66, 400, 88],
            _ => panic!(),
        };
        
        DenseMultilinSumcheckPoly::from_evaluations_vec(VARS, to_f(values))
    }
    
    #[test]
    fn sum_over_last_variable_test_3_vars() {
        let poly = test_poly::<3, Fr>();
        
        let poly = poly.sum_over_last_variable();
        
        assert_eq!(poly.num_vars(), 2);
        assert_eq!(poly.to_evaluations(), to_f::<Fr>(vec![30, 121, 700, 165]));

        let poly = poly.sum_over_last_variable();

        assert_eq!(poly.num_vars(), 1);
        assert_eq!(poly.to_evaluations(), to_f::<Fr>(vec![730, 286]));
    }

    #[test]
    fn test_partial_sum_poly_3var() {
        let poly = test_poly::<3, Fr>();

        let partial_sum_poly = poly.get_partial_sum_poly();
        assert_eq!(partial_sum_poly.coeffs, to_f::<Fr>(vec![730, -444]));
    }

    #[test]
    fn test_partial_sum_poly_2var() {
        let poly = test_poly::<2, Fr>();

        let partial_sum_poly = poly.get_partial_sum_poly();
        assert_eq!(partial_sum_poly.coeffs, to_f::<Fr>(vec![30, 670]));
    }
}