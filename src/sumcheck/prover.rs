use std::fmt::Debug;
use ark_ff::Field;
use crate::random_oracle::RandomOracle;
use crate::sumcheck::proof::{SumCheckProof, SumCheckStep};
use crate::sumcheck::sumcheck_poly::SumCheckPoly;

pub fn prove<F: Field, S: SumCheckPoly<F> + Clone, R: RandomOracle>(
    poly: &S,
    random_oracle: &R,
) -> SumCheckProof<F> {
    let random_points = vec![3,4,35,100].into_iter()
        .map(F::from)
        .collect::<Vec<_>>();
    let num_vars = poly.num_vars();

    // round 1
    let first_round_poly = poly.get_partial_sum_poly();
    // let g1_0 = first_round_poly.evaluate(&F::zero());
    // let g1_1 = first_round_poly.evaluate(&F::one());

    // assert_eq!(H, g1_0 + g1_1);

    let mut sc_steps = vec![];

    let mut current_poly = poly.clone();
    let mut current_r = random_points[0];

    for i in 1..num_vars {
        // let previous_partial_sum_poly = current_poly.get_partial_sum_poly();

        current_poly = current_poly.fix_variables(&[current_r]);
        let partial_sum_poly = current_poly.get_partial_sum_poly();

        // let gi_0 = partial_sum_poly.evaluate(&F::zero());
        // let gi_1 = partial_sum_poly.evaluate(&F::one());
        // let previous_poly_at_r = previous_partial_sum_poly.evaluate(&current_r);

        // assert_eq!(gi_1 + gi_0, previous_poly_at_r);

        sc_steps.push(SumCheckStep::new(partial_sum_poly, current_r));

        current_r = random_points[i];
    }

    // let mut r_vals = sc_steps
    //     .iter()
    //     .map(|step| step.r)
    //     .collect::<Vec<_>>();
    // r_vals.push(current_r);

    // let last_eval = poly.evaluate(&r_vals);

    // assert_eq!(sc_steps.last().unwrap().poly.evaluate(&current_r), last_eval);

    SumCheckProof::new(first_round_poly, current_r, sc_steps)
}