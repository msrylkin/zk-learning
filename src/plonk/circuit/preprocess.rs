use std::ops::Deref;
use ark_ec::pairing::Pairing;
use crate::kzg::KZG;
use crate::plonk::circuit::CompiledCircuit;

/// Preprocessed circuit with commitments to circuit polynomials.
pub struct PreprocessedCircuit<'a, P: Pairing> {
    /// Underlying compiled circuit.
    pub circuit: &'a CompiledCircuit<P::ScalarField>,
    /// Commitment to the `qm` polynomial.
    pub qm_comm: P::G1,
    /// Commitment to the `ql` polynomial.
    pub ql_comm: P::G1,
    /// Commitment to the `qr` polynomial.
    pub qr_comm: P::G1,
    /// Commitment to the `qo` polynomial.
    pub qo_comm: P::G1,
    /// Commitment to the `qc` polynomial.
    pub qc_comm: P::G1,
    /// Commitment to the `sigma_1` permutation polynomial.
    pub sigma_1_comm: P::G1,
    /// Commitment to the `sigma_2` permutation polynomial.
    pub sigma_2_comm: P::G1,
    /// Commitment to the `sigma_3` permutation polynomial.
    pub sigma_3_comm: P::G1,
}

impl<'a, P: Pairing> Deref for PreprocessedCircuit<'a, P> {
    type Target = CompiledCircuit<P::ScalarField>;

    fn deref(&self) -> &Self::Target {
        self.circuit
    }
}

/// Preprocesses a compiled circuit by committing to its polynomials.
///
/// All commitments are created using `kzg.commit(circuit.*)`.
pub fn preprocess_circuit<'a, P: Pairing>(
    circuit: &'a CompiledCircuit<P::ScalarField>,
    kzg: &KZG<P>,
) -> PreprocessedCircuit<'a, P> {
    PreprocessedCircuit {
        circuit,
        ql_comm: kzg.commit(&circuit.ql),
        qm_comm: kzg.commit(&circuit.qm),
        qr_comm: kzg.commit(&circuit.qr),
        qo_comm: kzg.commit(&circuit.qo),
        qc_comm: kzg.commit(&circuit.qc),
        sigma_1_comm: kzg.commit(&circuit.s_sigma_1),
        sigma_2_comm: kzg.commit(&circuit.s_sigma_2),
        sigma_3_comm: kzg.commit(&circuit.s_sigma_3),
    }
}