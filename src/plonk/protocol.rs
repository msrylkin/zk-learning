use ark_ec::pairing::Pairing;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use crate::kzg::KZG;
use crate::plonk::circuit::{preprocess_circuit, CompiledCircuit, PreprocessedCircuit, PublicWitness};
use crate::plonk::circuit::PlonkSolution;
use crate::plonk::domain::PlonkDomain;
use crate::plonk::proof::Proof;
use crate::plonk::prover::PlonkProver;
use crate::plonk::verifier::PlonkVerifier;

/// Common data shared by a protocol participant (Prover or Verifier).
///
/// Encapsulates the KZG commitment scheme, evaluation domain, vanishing polynomial,
/// preprocessed circuit with polynomial commitments, and the first Lagrange polynomial.
pub struct Party<'a, P: Pairing> {
    /// KZG commitment scheme with public parameters
    pub kzg: &'a KZG<P>,
    /// Plonk evaluation domain with precomputed cosets
    pub domain: &'a PlonkDomain<P::ScalarField>,
    /// Vanishing polynomial over the `self.domain`
    pub Zh: SparsePolynomial<P::ScalarField>,
    /// Preprocessed circuit containing commitments to circuit polynomials
    pub circuit: &'a PreprocessedCircuit<'a, P>,
    /// First Lagrange polynomial for the domain (`L1(X)`)
    pub lagrange_1: DensePolynomial<P::ScalarField>,
}

impl<'a, P: Pairing> Party<'a, P> {
    /// Creates a new `Party` object.
    /// The `circuit` must be preprocessed.
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

/// Main Plonk protocol object.
pub struct PlonkProtocol<P: Pairing> {
    kzg: KZG<P>,
    domain: PlonkDomain<P::ScalarField>,
    Zh: SparsePolynomial<P::ScalarField>,
    lagrange_1: DensePolynomial<P::ScalarField>,
}

impl<P: Pairing> PlonkProtocol<P> {
    /// Creates a new Plonk protocol with the provided `kzg` and Plonk `domain`
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

    /// Creates Plonk Instance for proving and verification.
    /// Preprocesses the `circuit` by commiting to all `circuit` polynomials
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

/// A prepared PLONK instance ready for proving and verifying.
///
/// Encapsulates a reference to the protocol and the preprocessed circuit.
pub struct PlonkInstance<'a, P: Pairing> {
    protocol: &'a PlonkProtocol<P>,
    circuit: PreprocessedCircuit<'a, P>,
}

impl<'a, P: Pairing> PlonkInstance<'a, P> {
    /// Generates a `Proof` for the provided `Solution`.
    pub fn prove(&self, solution: PlonkSolution<P::ScalarField>) -> Proof<P> {
        PlonkProver::new(Party::new(
            &self.protocol.kzg,
            &self.protocol.domain,
            &self.circuit
        )).prove(solution)
    }

    /// Verifies a Plonk `Proof` against the given `public_input`.
    ///
    /// # Panics
    /// Panics if verification fails.
    pub fn verify(&self, proof: &Proof<P>, public_input: &PublicWitness<P::ScalarField>) {
        PlonkVerifier::new(Party::new(
            &self.protocol.kzg,
            &self.protocol.domain,
            &self.circuit
        )).verify(public_input, proof)
    }
}