use ark_ff::Field;
use crate::gkr::circuit::{Circuit, Solution};
use crate::gkr::common::GKRProof;
use crate::gkr::prover::prove;
use crate::gkr::verifier::verify;
use crate::random_oracle::{RandomOracle};
use crate::sumcheck::SumCheckProtocol;

pub struct GKRProtocol<'a, R: RandomOracle> {
    random_oracle: &'a R,
    sum_check_protocol: SumCheckProtocol<'a, R>,
}

impl<'a, F: Field, R: RandomOracle<Item = F>> GKRProtocol<'a, R> {
    pub fn new(random_oracle: &'a R) -> Self {
        Self {
            random_oracle,
            sum_check_protocol: SumCheckProtocol::new(2, random_oracle),
        }
    }

    pub fn prove(
        &self,
        circuit: &Circuit<R::Item>,
        solution: &Solution<R::Item>,
    ) -> GKRProof<R::Item> {
        prove(circuit, solution, self.random_oracle, &self.sum_check_protocol)
    }
    
    pub fn verify(
        &self,
        circuit: &Circuit<R::Item>,
        proof: &GKRProof<R::Item>,
    ) {
        verify(circuit, proof, &self.sum_check_protocol)
    }
}
