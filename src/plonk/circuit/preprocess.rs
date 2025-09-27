use std::ops::Deref;
use ark_ec::pairing::Pairing;
use crate::kzg::KZG;
use crate::plonk::circuit::CompiledCircuit;

pub struct PreprocessedCircuit<'a, P: Pairing> {
    pub circuit: &'a CompiledCircuit<P::ScalarField>,
    pub qm_comm: P::G1,
    pub ql_comm: P::G1,
    pub qr_comm: P::G1,
    pub qo_comm: P::G1,
    pub qc_comm: P::G1,
    pub sigma_1_comm: P::G1,
    pub sigma_2_comm: P::G1,
    pub sigma_3_comm: P::G1,
}

impl<'a, P: Pairing> Deref for PreprocessedCircuit<'a, P> {
    type Target = CompiledCircuit<P::ScalarField>;

    fn deref(&self) -> &Self::Target {
        self.circuit
    }
}

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