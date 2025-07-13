use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use crate::gkr::circuit::Circuit;
use crate::gkr::common::{GKRProof, GKRProofLayer};
use crate::poly_utils::{interpolate, line};
use crate::sumcheck;
use crate::sumcheck::{OracleEvaluation, SumCheckProtocol};

#[derive(Debug, Clone)]
struct GKRFinalOracle<F: Field> {
    mul_i: DenseMultilinearExtension<F>,
    add_i: DenseMultilinearExtension<F>,
    q: DensePolynomial<F>
}

impl<F: Field> GKRFinalOracle<F> {
    fn new(add_i: DenseMultilinearExtension<F>, mul_i: DenseMultilinearExtension<F>, q: DensePolynomial<F>) -> Self {
        Self { mul_i, add_i, q }
    }
}

impl<F: Field> OracleEvaluation<F> for GKRFinalOracle<F> {
    fn final_eval(&self, r: &[F]) -> F {
        let q_0 = self.q.evaluate(&F::zero());
        let q_1 = self.q.evaluate(&F::one());
        let add_eval = ark_poly::Polynomial::evaluate(&self.add_i, &r.to_vec()) * (q_0 + q_1);
        let mul_eval = ark_poly::Polynomial::evaluate(&self.mul_i, &r.to_vec()) * (q_0 * q_1);

        add_eval + mul_eval
    }
}

pub fn verify<F: Field>(
    circuit: &Circuit<F>,
    proof: &GKRProof<F>,
    sum_check_protocol: &SumCheckProtocol,
) {
    let mut mi = Polynomial::evaluate(&proof.W0, &proof.r0);
    let mut ri = proof.r0.clone();

    for (i, GKRProofLayer { q, r_star, sumcheck_proof }) in proof.layers.iter().enumerate().rev() {
        let add_poly = circuit.add_i(i);
        let mul_poly = circuit.mul_i(i);
        let add_fixed = MultilinearExtension::fix_variables(&add_poly, &ri);
        let mul_fixed = MultilinearExtension::fix_variables(&mul_poly, &ri);
        let used_r = sumcheck_proof.get_used_randomness();
        let (b, c) = used_r.split_at(used_r.len() / 2);
        let l = line(b, c);

        let final_oracle = GKRFinalOracle::new(add_fixed, mul_fixed, q.clone());

        // sumcheck::verify(&final_oracle, &sumcheck_proof, mi);
        sum_check_protocol.verify(&final_oracle, &sumcheck_proof, mi);

        ri = l.iter().map(|li| li.evaluate(r_star)).collect::<Vec<_>>();
        mi = q.evaluate(r_star);
    }

    let last_w = interpolate(&proof.inputs);
    let last_w_r = Polynomial::evaluate(&last_w, &ri.to_vec());

    assert_eq!(last_w_r, mi);
}