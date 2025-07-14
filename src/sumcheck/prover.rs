use ark_ff::Field;
use crate::random_oracle::RandomOracle;
use crate::sumcheck::proof::{SumCheckProof, SumCheckStep};
use crate::sumcheck::sumcheck_poly::SumCheckPoly;

pub fn prove<F: Field, S: SumCheckPoly<F> + Clone, R: RandomOracle<Item = F>>(
    poly: &S,
    random_oracle: &R,
) -> SumCheckProof<F> {
    let num_vars = poly.num_vars();

    let first_round_poly = poly.get_partial_sum_poly();

    let mut sc_steps = vec![];

    let mut current_poly = poly.clone();
    let mut current_r = random_oracle.get_randomness(1)[0];

    for _ in 1..num_vars {
        current_poly = current_poly.fix_variables(&[current_r]);
        let partial_sum_poly = current_poly.get_partial_sum_poly();

        sc_steps.push(SumCheckStep::new(partial_sum_poly, current_r));

        current_r = random_oracle.get_randomness(1)[0];
    }

    SumCheckProof::new(first_round_poly, current_r, sc_steps)
}