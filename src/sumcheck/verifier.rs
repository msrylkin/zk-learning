use ark_ff::Field;
use ark_poly::Polynomial;
use crate::sumcheck::proof::SumCheckProof;
use crate::sumcheck::sumcheck_poly::OracleEvaluation;

pub fn verify<F: Field, O: OracleEvaluation<F>>(
    oracle: &O,
    proof: &SumCheckProof<F>,
    H: F,
    max_step_partial_poly_degree: usize,
) {
    let g1_0 = proof.first_round_poly.evaluate(&F::zero());
    let g1_1 = proof.first_round_poly.evaluate(&F::one());

    assert_eq!(g1_0 + g1_1, H);
    assert!(proof.first_round_poly.degree() <= max_step_partial_poly_degree);

    let mut previous_poly = &proof.first_round_poly;

    for step in &proof.steps {
        let previous_poly_at_r = previous_poly.evaluate(&step.r);
        let gi_0 = step.poly.evaluate(&F::zero());
        let gi_1 = step.poly.evaluate(&F::one());

        assert_eq!(gi_0 + gi_1, previous_poly_at_r);
        assert!(step.poly.degree() <= max_step_partial_poly_degree);

        previous_poly = &step.poly;
    }

    let r_vals = proof.get_used_randomness();
    let last_eval = oracle.final_eval(&r_vals);

    assert_eq!(last_eval, proof.steps.last().unwrap().poly.evaluate(&proof.last_round_r));
}