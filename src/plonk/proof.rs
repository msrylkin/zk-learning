use ark_ec::pairing::Pairing;
use ark_ff::{FftField, PrimeField};

#[derive(Clone, Debug)]
pub struct Openings<F: PrimeField + FftField> {
    pub a: F,
    pub b: F,
    pub c: F,
    pub s_sigma_1: F,
    pub s_sigma_2: F,
    pub z_shifted: F,
}

#[derive(Clone, Debug)]
pub struct Commitments<P: Pairing> {
    pub a: P::G1,
    pub b: P::G1,
    pub c: P::G1,
    pub z: P::G1,
    pub t_lo: P::G1,
    pub t_mid: P::G1,
    pub t_hi: P::G1,
}

pub struct OpeningProofs<P: Pairing> {
    pub w_zeta: P::G1,
    pub w_zeta_omega: P::G1,
}

pub struct Proof<P: Pairing> {
    pub openings: Openings<P::ScalarField>,
    pub commitments: Commitments<P>,
    pub opening_proofs: OpeningProofs<P>,
}

impl<P: Pairing> Proof<P> {
    pub fn new(
        openings: Openings<P::ScalarField>,
        commitments: Commitments<P>,
        opening_proofs: OpeningProofs<P>,
    ) -> Self {
        Self {
            openings,
            commitments,
            opening_proofs,
        }
    }
}
