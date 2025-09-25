use ark_ec::pairing::Pairing;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use crate::kzg::KZG;
use crate::plonk::circuit::{preprocess_circuit, CompiledCircuit, PreprocessedCircuit, PublicWitness};
use crate::plonk::circuit::solution::Solution;
use crate::plonk::domain::PlonkDomain;
use crate::plonk::proof::Proof;
use crate::plonk::prover::PlonkProver;
use crate::plonk::verifier::PlonkVerifier;

pub struct Party<'a, P: Pairing> {
    pub kzg: &'a KZG<P>,
    pub domain: &'a PlonkDomain<P::ScalarField>,
    pub Zh: SparsePolynomial<P::ScalarField>,
    pub circuit: &'a PreprocessedCircuit<'a, P>,
    pub lagrange_1: DensePolynomial<P::ScalarField>,
}

impl<'a, P: Pairing> Party<'a, P> {
    pub fn new(
        kzg: &'a KZG<P>,
        domain: &'a PlonkDomain<P::ScalarField>,
        circuit: &'a PreprocessedCircuit<'a, P>,
    ) -> Self {
        Self {
            kzg,
            Zh: domain.get_vanishing_polynomial(),
            lagrange_1: domain.lagrange_polys().first().unwrap().clone(),
            domain,
            circuit,
        }
    }
}

pub struct PlonkProtocol<P: Pairing> {
    kzg: KZG<P>,
    domain: PlonkDomain<P::ScalarField>,
    Zh: SparsePolynomial<P::ScalarField>,
    lagrange_1: DensePolynomial<P::ScalarField>,
}

impl<P: Pairing> PlonkProtocol<P> {
    pub fn new(
        kzg: KZG<P>,
        domain: PlonkDomain<P::ScalarField>,
    ) -> Self {
        Self {
            kzg,
            Zh: domain.get_vanishing_polynomial(),
            lagrange_1: domain.lagrange_polys().first().unwrap().clone(),
            domain,
        }
    }

    pub fn create_instance<'a>(
        &'a self,
        circuit: &'a CompiledCircuit<P::ScalarField>,
    ) -> PlonkInstance<'a, P> {
        PlonkInstance {
            protocol: self,
            circuit: preprocess_circuit(circuit, &self.kzg),
        }
    }
}

pub struct PlonkInstance<'a, P: Pairing> {
    protocol: &'a PlonkProtocol<P>,
    circuit: PreprocessedCircuit<'a, P>,
}

impl<'a, P: Pairing> PlonkInstance<'a, P> {
    pub fn prove(&self, solution: Solution<P::ScalarField>) -> Proof<P> {
        PlonkProver::new(Party::new(
            &self.protocol.kzg,
            &self.protocol.domain,
            &self.circuit
        )).prove(solution)
    }

    pub fn verify(&self, proof: &Proof<P>, public_input: &PublicWitness<P::ScalarField>) {
        PlonkVerifier::new(Party::new(
            &self.protocol.kzg,
            &self.protocol.domain,
            &self.circuit
        )).verify(public_input, proof)
    }
}