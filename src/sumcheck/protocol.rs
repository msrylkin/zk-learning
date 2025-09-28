use crate::random_oracle::RandomOracle;
use crate::sumcheck::prover::prove;
use crate::sumcheck::{OracleEvaluation, SumCheckPoly, SumCheckProof};
use crate::sumcheck::verifier::verify;

/// SumCheck protocol object.
///
/// Handles both proving and verification for multivariate
/// sumcheck polynomials. Uses a provided `RandomOracle`
/// for deriving verifier randomness during the protocol.
pub struct SumCheckProtocol<'a, R: RandomOracle> {
    /// Maximum degree of partial sum polynomials at any step.
    max_step_partial_poly_degree: usize,
    /// Reference to the randomness oracle used in the protocol.
    random_oracle: &'a R,
}

impl<'a, R: RandomOracle> SumCheckProtocol<'a, R> {
    /// Creates a new `SumCheckProtocol` instance with a specified
    /// maximum partial polynomial degree and a `RandomOracle`.
    pub fn new(
        max_step_partial_poly_degree: usize,
        random_oracle: &'a R,
    ) -> Self {
        Self { max_step_partial_poly_degree, random_oracle }
    }

    /// Executes the sumcheck prover for a given sumncheck polynomial `poly`.
    ///
    /// Returns the generated `SumCheckProof`.
    pub fn prove<S: SumCheckPoly<R::Item> + Clone>(
        &self,
        poly: &S,
    ) -> SumCheckProof<R::Item> {
        prove(poly, self.random_oracle)
    }

    /// Verifies a sumcheck proof against a claimed sum using
    /// the final evaluation oracle `final_evaluation_oracle`.
    ///
    /// # Panics
    /// if the proof is invalid.
    pub fn verify<E: OracleEvaluation<R::Item>>(
        &self,
        final_evaluation_oracle: &E,
        proof: &SumCheckProof<R::Item>,
        claimed_sum: R::Item,
    ) {
        verify(final_evaluation_oracle, proof, claimed_sum, self.max_step_partial_poly_degree)
    }
}