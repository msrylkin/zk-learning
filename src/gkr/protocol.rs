use ark_ff::Field;
use crate::gkr::circuit::{Circuit, Solution};
use crate::gkr::common::GKRProof;
use crate::gkr::prover::prove;
use crate::gkr::verifier::verify;
use crate::gkr::random_oracle::RandomOracle;

struct GKRProtocol<O> {
    random_oracle: O,
}

impl<F: Field, O: RandomOracle<Item = F>> GKRProtocol<O> {
    pub fn new(random_oracle: O) -> Self {
        Self {
            random_oracle,
        }
    }

    pub fn prove(
        &self,
        circuit: &Circuit<F>,
        solution: &Solution<F>,
    ) -> GKRProof<F> {
        prove(circuit, solution, &self.random_oracle)
    }
    
    pub fn verify(
        circuit: &Circuit<F>,
        proof: &GKRProof<F>,
    ) {
        verify(circuit, proof)
    }
}
