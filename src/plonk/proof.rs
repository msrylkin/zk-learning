use ark_ec::pairing::Pairing;
use ark_ff::{FftField, PrimeField};

/// Evaluations of polynomials at a specific point during the protocol.
///
/// Each field represents the evaluation of a specific column polynomial or permutation polynomial.
#[derive(Clone, Debug)]
pub struct Openings<F: PrimeField + FftField> {
    /// `a` poly evaluation (left column)
    pub a: F,
    /// `b` poly evaluation (right column)
    pub b: F,
    /// `c` poly evaluation (output column)
    pub c: F,
    /// `sigma_1` poly evaluation (left column permutation poly)
    pub s_sigma_1: F,
    /// `sigma_2` poly evaluation (right column permutation poly)
    pub s_sigma_2: F,
    /// `z_shifted` poly evaluation (grand product poly shifted by 1 evaluation)
    pub z_shifted: F,
}


/// Polynomials commitments.
/// Each commitment is a value in G1
#[derive(Clone, Debug)]
pub struct Commitments<P: Pairing> {
    /// `a` poly commitment (left column)
    pub a: P::G1,
    /// `b` poly commitment (left column)
    pub b: P::G1,
    /// `c` poly commitment (left column)
    pub c: P::G1,
    /// `z` poly commitment (grand permutation product poly)
    pub z: P::G1,
    /// `t_lo` poly commitment (first part of big quotient poly)
    pub t_lo: P::G1,
    /// `t_mid` poly commitment (second part of big quotient poly)
    pub t_mid: P::G1,
    /// `t_hi` poly commitment (third part of big quotient poly)
    pub t_hi: P::G1,
}

/// 2 Opening proofs for zeta and zeta * omega
pub struct OpeningProofs<P: Pairing> {
    /// `zeta` opening proof
    pub w_zeta: P::G1,
    /// `zeta * omega` opening proof
    pub w_zeta_omega: P::G1,
}

/// Complete PLONK proof structure.
///
/// Combines polynomial openings, commitments, and the corresponding opening proofs.
pub struct Proof<P: Pairing> {
    /// Opened polynomial values at `zeta` challenge point
    pub openings: Openings<P::ScalarField>,
    /// Commitments to polynomials used in the protocol
    pub commitments: Commitments<P>,
    /// Opening proofs for `zeta` and `zeta * omega`
    pub opening_proofs: OpeningProofs<P>,
}

impl<P: Pairing> Proof<P> {
    /// Constructs a new PLONK proof from provided openings, commitments, and opening proofs.
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
