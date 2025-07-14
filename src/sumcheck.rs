mod protocol;
mod prover;
mod verifier;
mod sumcheck_poly;
mod multilin_sumcheck_poly;
mod proof;

pub use proof::*;
pub use sumcheck_poly::*;
pub use protocol::*;

use std::fmt::Debug;
use ark_ff::{Field, One};
use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_std::Zero;
use ark_test_curves::bls12_381::Fr;
use crate::random_oracle::FixedRandomOracle;
use crate::sumcheck::multilin_sumcheck_poly::{DenseMultilinFinalEvaluationOracle, DenseMultilinSumcheckPoly};

pub fn test_sc() {
    let mut poly = DenseMultilinearExtension::from_evaluations_vec(
        4,
        vec![91, 62, 13, 431, 98, 123, 2871, 7, 512, 12, 63 ,982, 474, 2847, 912, 744].into_iter().map(Fr::from).collect()
    );
    
    let sum = poly.iter().sum::<Fr>();
    let dense_sc_poly = DenseMultilinSumcheckPoly::new(poly);
    let random_oracle = FixedRandomOracle::new(vec![10, 20, 30, 40, 50].into_iter().map(Fr::from).collect());
    let sc_protocol = SumCheckProtocol::new(
        1,
        &random_oracle
    );

    let sc_proof = sc_protocol.prove(&dense_sc_poly);
    
    println!("\nnew res {:?}", sc_proof);
    
    let final_evaluation_oracle = DenseMultilinFinalEvaluationOracle::new(dense_sc_poly);
    
    sc_protocol.verify(&final_evaluation_oracle, &sc_proof, sum);
    
    println!("poly proved");
}

// impl<F: Field> SumCheckPoly<F> for DenseMultilinearExtension<F> {
//     fn get_evaluations(&self) -> Vec<F> {
//         self.to_evaluations()
//     }
// 
//     fn num_vars(&self) -> usize {
//         self.num_vars
//     }
// 
//     fn get_partial_sum_poly(&self) -> DensePolynomial<F> {
//         let sum_vars = SumCheckPoly::num_vars(self) - 1;
//         let mut current_poly = self.clone();
//     
//         for i in 0..sum_vars {
//             let summed_poly = sum_over_last_variable(&current_poly);
//             current_poly = summed_poly;
//         }
//     
//         interpolate_univariate_on_evals(&current_poly.evaluations.try_into().unwrap())
//     }
// 
//     // fn get_partial_sum_poly(&self) -> DensePolynomial<F> {
//     //     let mut current_sc_poly = self.clone();
//     //     let num_vars = MultilinearExtension::num_vars(&current_sc_poly);
//     // 
//     //     let evals = self.get_evaluations();
//     //     println!("evals {:?}", evals);
//     // 
//     //     let (zeroes, ones) = evals.iter().enumerate().fold((F::zero(), F::zero()), |(zeroes, ones), (i, eval)| {
//     //         if i % 2 == 0 {
//     //             (zeroes + eval, ones)
//     //         } else {
//     //             (zeroes, ones + eval)
//     //         }
//     //     });
//     // 
//     //     interpolate_univariate_on_evals(&[zeroes, ones])
//     // }
// 
//     fn fix_variables(&self, partial_point: &[F]) -> Self {
//         MultilinearExtension::fix_variables(self, partial_point)
//     }
// 
//     fn evaluate(&self, point: &[F]) -> F {
//         Polynomial::evaluate(self, &point.to_vec())
//     }
// }

// pub trait OracleEvaluation<F: Field> {
//     fn final_eval(&self, r: &[F]) -> F;
// }
// 
// pub trait SumCheckPoly<F: Field> {
//     fn get_evaluations(&self) -> Vec<F>;
// 
//     fn num_vars(&self) -> usize;
// 
//     fn get_partial_sum_poly(&self) -> DensePolynomial<F>;
// 
//     fn fix_variables(&self, partial_point: &[F]) -> Self;
// 
//     fn evaluate(&self, point: &[F]) -> F;
// }

// fn prove_for_val<F: Field, S: SumCheckPoly<F> + Clone>(poly: S) {
//     let mut current_poly = poly.clone();
//     // let mut current_poly = poly;
//     // TODO: Fiat-Shamir
//     let random_points = vec![21,23,35,48].into_iter()
//         .map(F::from)
//         .collect::<Vec<_>>();
//     let mut current_sum = current_poly.get_evaluations().iter().sum::<F>();
//     let mut total_rounds = 0;
// 
//     for round in 0..poly.num_vars() {
//         // let partial_sum_poly = get_partial_sum_poly(&current_poly);
//         let partial_sum_poly = current_poly.get_partial_sum_poly();
// 
//         if partial_sum_poly.degree() != 1 {
//             panic!("degree does not match");
//         }
// 
//         let eval_0 = partial_sum_poly.evaluate(&F::zero());
//         let eval_1 = partial_sum_poly.evaluate(&F::one());
// 
//         let eval_sum = eval_0 + eval_1;
// 
//         if eval_sum != current_sum {
//             panic!("sum does not match");
//         }
// 
//         current_poly = current_poly.fix_variables(&[random_points[round]]);
//         current_sum = current_poly.get_evaluations().iter().sum::<F>();
//         total_rounds += 1;
//     }
// 
//     let last_evaluation = current_poly.get_evaluations().pop().unwrap();
//     let points = random_points.into_iter().take(total_rounds).collect::<Vec<_>>();
//     let original_poly_evaluation = poly.evaluate(&points);
// 
//     if last_evaluation != original_poly_evaluation {
//         panic!("last evaluation not match");
//     }
//     // println!("")
// }

// #[derive(Debug)]
// pub struct SumCheckStep<F: Field> {
//     pub poly: DensePolynomial<F>,
//     pub r: F,
// }
// 
// impl<F: Field> SumCheckStep<F> {
//     fn new(poly: DensePolynomial<F>, r: F) -> Self {
//         SumCheckStep {
//             poly,
//             r
//         }
//     }
// }
// 
// #[derive(Debug)]
// pub struct SumCheckProof<F: Field> {
//     pub first_round_poly: DensePolynomial<F>,
//     pub last_round_r: F,
//     pub steps: Vec<SumCheckStep<F>>
// }
// 
// impl<F: Field> SumCheckProof<F> {
//     pub fn get_used_randomness(&self) -> Vec<F> {
//         let mut used_r = self.steps
//             .iter()
//             .map(|step| step.r)
//             .collect::<Vec<_>>();
//         used_r.push(self.last_round_r);
//         
//         used_r
//     }
// } 
// 
// impl<F: Field> SumCheckProof<F> {
//     fn new(first_round_poly: DensePolynomial<F>, last_round_r: F, steps: Vec<SumCheckStep<F>>) -> Self {
//         SumCheckProof {
//             first_round_poly,
//             last_round_r,
//             steps,
//         }
//     }
// }
// 
// pub fn prove<F: Field, S: SumCheckPoly<F> + Clone + Debug>(poly: &S) -> SumCheckProof<F> {
//     let random_points = vec![3,4,35,100].into_iter()
//         .map(F::from)
//         .collect::<Vec<_>>();
//     let num_vars = poly.num_vars();
// 
//     // round 1
//     let first_round_poly = poly.get_partial_sum_poly();
//     // let g1_0 = first_round_poly.evaluate(&F::zero());
//     // let g1_1 = first_round_poly.evaluate(&F::one());
// 
//     // assert_eq!(H, g1_0 + g1_1);
//     
//     let mut sc_steps = vec![];
//     
//     let mut current_poly = poly.clone();
//     let mut current_r = random_points[0];
// 
//     for i in 1..num_vars {
//         // let previous_partial_sum_poly = current_poly.get_partial_sum_poly();
//         
//         current_poly = current_poly.fix_variables(&[current_r]);
//         let partial_sum_poly = current_poly.get_partial_sum_poly();
//         
//         // let gi_0 = partial_sum_poly.evaluate(&F::zero());
//         // let gi_1 = partial_sum_poly.evaluate(&F::one());
//         // let previous_poly_at_r = previous_partial_sum_poly.evaluate(&current_r);
// 
//         // assert_eq!(gi_1 + gi_0, previous_poly_at_r);
//         
//         sc_steps.push(SumCheckStep::new(partial_sum_poly, current_r));
//         
//         current_r = random_points[i];
//     }
//     
//     // let mut r_vals = sc_steps
//     //     .iter()
//     //     .map(|step| step.r)
//     //     .collect::<Vec<_>>();
//     // r_vals.push(current_r);
// 
//     // let last_eval = poly.evaluate(&r_vals);
// 
//     // assert_eq!(sc_steps.last().unwrap().poly.evaluate(&current_r), last_eval);
// 
//     SumCheckProof::new(first_round_poly, current_r, sc_steps)
// }

// pub fn verify<F: Field, O: OracleEvaluation<F>>(oracle: &O, proof: &SumCheckProof<F>, H: F) {
//     let g1_0 = proof.first_round_poly.evaluate(&F::zero());
//     let g1_1 = proof.first_round_poly.evaluate(&F::one());
//     
//     assert_eq!(g1_0 + g1_1, H);
//     // check degree
//     
//     let mut previous_poly = &proof.first_round_poly;
//     
//     for step in &proof.steps {
//         let previous_poly_at_r = previous_poly.evaluate(&step.r);
//         let gi_0 = step.poly.evaluate(&F::zero());
//         let gi_1 = step.poly.evaluate(&F::one());
//         
//         assert_eq!(gi_0 + gi_1, previous_poly_at_r);
//         // check degree
//         
//         previous_poly = &step.poly;
//     }
// 
//     let mut r_vals = proof.steps
//         .iter()
//         .map(|step| step.r)
//         .collect::<Vec<_>>();
//     r_vals.push(proof.last_round_r);
//     
//     // let last_eval = poly.evaluate(&r_vals);
//     let last_eval = oracle.final_eval(&r_vals);
//     
//     assert_eq!(last_eval, proof.steps.last().unwrap().poly.evaluate(&proof.last_round_r));
// }

// pub fn get_partial_sum_poly(
//     poly: &DenseMultilinearExtension<Fr>,
// ) -> DensePolynomial<Fr> {
//     let sum_vars = SumCheckPoly::num_vars(poly) - 1;
//     let mut current_poly = poly.clone();
// 
//     for i in 0..sum_vars {
//         let summed_poly = sum_over_last_variable(&current_poly);
//         current_poly = summed_poly;
//     }
// 
//     interpolate_univariate_on_evals(&current_poly.evaluations.try_into().unwrap())
// }

// // 000 |
// // 001 --> 00[1 + 0]
// // 010 |
// // 011 --> 01[1 + 0]
// pub fn sum_over_last_variable<F: Field>(poly: &DenseMultilinearExtension<F>) -> DenseMultilinearExtension<F> {
//     let sums = poly
//         .evaluations
//         .iter()
//         .enumerate()
//         .fold(vec![F::zero(); poly.evaluations.len() / 2], |mut sums, (i, e)| {
//             if i % 2 == 0 {
//                 sums[(i / 4) * 2] += e;
//             } else {
//                 sums[(i / 4) * 2 + 1] += e;
//             }
//             sums
//         });
// 
//     let res = DenseMultilinearExtension::from_evaluations_vec(poly.num_vars - 1, sums);
// 
//     res
// }

// // f(x) = (1 - x) * a + x * b ==> a - x * a + x * b ==> a + x * (b - a)
// pub fn interpolate_univariate_on_evals<F: Field>(
//     evals: &[F; 2]
// ) -> DensePolynomial<F> {
//     let a = evals[0];
//     let b = evals[1];
// 
//     DensePolynomial::from_coefficients_vec(vec![a, b - a])
// }