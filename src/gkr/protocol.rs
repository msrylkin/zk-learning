use ark_ff::Field;
use crate::gkr::circuit::{Circuit, GkrSolution};
use crate::gkr::proof::GKRProof;
use crate::gkr::prover::prove;
use crate::gkr::verifier::verify;
use crate::random_oracle::{RandomOracle};
use crate::sumcheck::SumCheckProtocol;

/// Main GKR protocol object. Encapsulates the sumcheck protocol
/// and randomness oracle used for Fiat-Shamir challenges.
pub struct GKRProtocol<'a, R: RandomOracle> {
    random_oracle: &'a R,
    sum_check_protocol: SumCheckProtocol<'a, R>,
}

impl<'a, F: Field, R: RandomOracle<Item = F>> GKRProtocol<'a, R> {
    /// Creates a new GKR protocol instance with the given `random_oracle`.
    pub fn new(random_oracle: &'a R) -> Self {
        Self {
            random_oracle,
            sum_check_protocol: SumCheckProtocol::new(2, random_oracle),
        }
    }

    /// Generates a `GKRProof` for the given `circuit` and its `solution`.
    pub fn prove(
        &self,
        circuit: &Circuit<R::Item>,
        solution: &GkrSolution<R::Item>,
    ) -> GKRProof<R::Item> {
        prove(circuit, solution, self.random_oracle, &self.sum_check_protocol)
    }

    /// Verifies a GKR proof against the given circuit.
    ///
    /// # Panics
    /// if verification fails
    pub fn verify(
        &self,
        circuit: &Circuit<R::Item>,
        proof: &GKRProof<R::Item>,
    ) {
        verify(circuit, proof, &self.sum_check_protocol)
    }
}
