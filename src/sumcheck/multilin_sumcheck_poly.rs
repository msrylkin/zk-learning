use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_test_curves::bls12_381::Fr;
use crate::poly_utils::interpolate_univariate_on_evals;
// use crate::sumcheck::{interpolate_univariate_on_evals, sum_over_last_variable};
use crate::sumcheck::sumcheck_poly::SumCheckPoly;

fn get_partial_sum_poly(
    poly: &DenseMultilinearExtension<Fr>,
) -> DensePolynomial<Fr> {
    let sum_vars = SumCheckPoly::num_vars(poly) - 1;
    let mut current_poly = poly.clone();

    for i in 0..sum_vars {
        let summed_poly = sum_over_last_variable(&current_poly);
        current_poly = summed_poly;
    }

    interpolate_univariate_on_evals(&current_poly.evaluations.try_into().unwrap())
}

// 000 |
// 001 --> 00[1 + 0]
// 010 |
// 011 --> 01[1 + 0]
fn sum_over_last_variable<F: Field>(poly: &DenseMultilinearExtension<F>) -> DenseMultilinearExtension<F> {
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

    let res = DenseMultilinearExtension::from_evaluations_vec(poly.num_vars - 1, sums);

    res
}

impl<F: Field> SumCheckPoly<F> for DenseMultilinearExtension<F> {
    fn get_evaluations(&self) -> Vec<F> {
        self.to_evaluations()
    }

    fn num_vars(&self) -> usize {
        self.num_vars
    }

    fn get_partial_sum_poly(&self) -> DensePolynomial<F> {
        let sum_vars = SumCheckPoly::num_vars(self) - 1;
        let mut current_poly = self.clone();

        for i in 0..sum_vars {
            let summed_poly = sum_over_last_variable(&current_poly);
            current_poly = summed_poly;
        }

        interpolate_univariate_on_evals(&current_poly.evaluations.try_into().unwrap())
    }

    // fn get_partial_sum_poly(&self) -> DensePolynomial<F> {
    //     let mut current_sc_poly = self.clone();
    //     let num_vars = MultilinearExtension::num_vars(&current_sc_poly);
    // 
    //     let evals = self.get_evaluations();
    //     println!("evals {:?}", evals);
    // 
    //     let (zeroes, ones) = evals.iter().enumerate().fold((F::zero(), F::zero()), |(zeroes, ones), (i, eval)| {
    //         if i % 2 == 0 {
    //             (zeroes + eval, ones)
    //         } else {
    //             (zeroes, ones + eval)
    //         }
    //     });
    // 
    //     interpolate_univariate_on_evals(&[zeroes, ones])
    // }

    fn fix_variables(&self, partial_point: &[F]) -> Self {
        MultilinearExtension::fix_variables(self, partial_point)
    }

    fn evaluate(&self, point: &[F]) -> F {
        Polynomial::evaluate(self, &point.to_vec())
    }
}