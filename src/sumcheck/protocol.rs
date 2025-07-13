use ark_ff::Field;
use crate::random_oracle::RandomOracle;
use crate::sumcheck::prover::prove;
use crate::sumcheck::{OracleEvaluation, SumCheckPoly, SumCheckProof};
use crate::sumcheck::verifier::verify;

pub struct SumCheckProtocol<R> {
    max_step_partial_poly_degree: usize,
    random_oracle: R,
}

impl<R: RandomOracle> SumCheckProtocol<R> {
    pub fn new(
        max_step_partial_poly_degree: usize,
        random_oracle: R,
    ) -> Self {
        Self { max_step_partial_poly_degree, random_oracle }
    }
    
    pub fn prove<S: SumCheckPoly<R::Item> + Clone>(
        &self,
        poly: &S,
    ) -> SumCheckProof<R::Item> {
        prove(poly, &self.random_oracle)
    }
    
    pub fn verify<E: OracleEvaluation<R::Item>>(
        &self,
        final_evaluation_oracle: &E,
        proof: &SumCheckProof<R::Item>,
        claimed_sum: R::Item,
    ) {
        verify(final_evaluation_oracle, proof, claimed_sum, self.max_step_partial_poly_degree)
    }
}