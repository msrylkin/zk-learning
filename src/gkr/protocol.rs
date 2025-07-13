use std::pin::{pin, Pin};
use std::rc::Rc;
use ark_ff::Field;
use crate::gkr::circuit::{Circuit, Solution};
use crate::gkr::common::GKRProof;
use crate::gkr::prover::prove;
use crate::gkr::verifier::verify;
use crate::random_oracle::{FixedRandomOracle, ProxyRandomOracle, RandomOracle};
use crate::sumcheck::SumCheckProtocol;

struct GKRProtocol<'a, R: RandomOracle> {
    random_oracle: R,
    sum_check_protocol: SumCheckProtocol<ProxyRandomOracle<'a, R>>,
}

impl<F: Field, R: RandomOracle<Item = F>> GKRProtocol<'_, R> {
    pub fn new(random_oracle: &R) -> Self {
        let random_oracle = pin!(random_oracle);
        let proxy_oracle = ProxyRandomOracle::new(random_oracle);

        Self {
            random_oracle,
            sum_check_protocol: SumCheckProtocol::new(2, proxy_oracle),
        }
    }

    pub fn prove(
        &self,
        circuit: &Circuit<R::Item>,
        solution: &Solution<R::Item>,
    ) -> GKRProof<R::Item> {
        prove(circuit, solution, &self.random_oracle)
    }
    
    pub fn verify(
        &self,
        circuit: &Circuit<R::Item>,
        proof: &GKRProof<R::Item>,
    ) {
        verify(circuit, proof, &self.sum_check_protocol)
    }
}
