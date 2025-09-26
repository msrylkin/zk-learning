use std::ops::Deref;
use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, PrimeGroup};
use ark_ff::{Field, PrimeField};
use ark_poly::Polynomial;
use ark_poly::univariate::{DensePolynomial};
use ark_std::{One, Zero};
use crate::kzg::{BatchOpening, MultipointOpening};
use crate::plonk::circuit::{PublicWitness};
use crate::plonk::proof::Proof;
use crate::plonk::protocol::Party;
use crate::plonk::transcript_protocol::TranscriptProtocol;
use crate::poly_utils::{const_poly};

pub struct PlonkVerifier<'a, P: Pairing> {
    party_data: Party<'a, P>,
}

impl<'a, P: Pairing> Deref for PlonkVerifier<'a, P> {
    type Target = Party<'a, P>;

    fn deref(&self) -> &Self::Target {
        &self.party_data
    }
}

struct DerivedChallenges<F: PrimeField> {
    beta: F,
    gamma: F,
    alpha: F,
    zeta: F,
    vi: F,
    u: F,
}

impl<'a, P: Pairing> PlonkVerifier<'a, P> {
    pub fn new(
        party_data: Party<'a, P>,
    ) -> Self {
        Self {
            party_data
        }
    }

    pub fn verify(
        &self,
        public_witness: &PublicWitness<P::ScalarField>,
        proof: &Proof<P>,
    ) {
        let omega = self.domain.generator();

        let challenges = derive_challenges(
            omega,
            &proof,
            &public_witness.pi_vector,
            &public_witness.output_vector,
        );

        let D = self.compute_linearized_commitment(
            &proof,
            &public_witness.pi_combined,
            &challenges
        );

        let is_valid = self.kzg.check_multipoint(
            &[
                MultipointOpening {
                    opening_point: &challenges.zeta,
                    commitments: &[D, proof.commitments.a, proof.commitments.b, proof.commitments.c, self.circuit.sigma_1_comm, self.circuit.sigma_2_comm],
                    batch_opening: &BatchOpening::new(
                        vec![P::ScalarField::zero(), proof.openings.a, proof.openings.b, proof.openings.c, proof.openings.s_sigma_1, proof.openings.s_sigma_2],
                        proof.opening_proofs.w_zeta,
                    ),
                    linearization_scalar: &challenges.vi
                },
                MultipointOpening {
                    opening_point: &(challenges.zeta * omega),
                    commitments: &[proof.commitments.z],
                    batch_opening: &BatchOpening::new(
                        vec![proof.openings.z_shifted],
                        proof.opening_proofs.w_zeta_omega,
                    ),
                    linearization_scalar: &P::ScalarField::one(),
                },
            ],
            &challenges.u,
        );

        assert!(is_valid);
    }

    fn compute_linearized_commitment(
        &self,
        proof: &Proof<P>,
        public_input_poly: &DensePolynomial<P::ScalarField>,
        challenges: &DerivedChallenges<P::ScalarField>,
    ) -> P::G1 {
        let (k1, k2) = (self.domain.k1(), self.domain.k2());
        let beta = &challenges.beta;
        let gamma = &challenges.gamma;
        let alpha = &challenges.alpha;
        let zeta = &challenges.zeta;

        // gates
        let mut D = self.circuit.qm_comm * proof.openings.a * proof.openings.b
            + self.circuit.ql_comm * proof.openings.a
            + self.circuit.qr_comm * proof.openings.b
            + self.circuit.qo_comm * proof.openings.c
            + P::G1::generator() * public_input_poly.evaluate(zeta)
            + self.circuit.qc_comm;

        // numerator
        D += proof.commitments.z * alpha * (proof.openings.a + *beta * *zeta + gamma)
            * (proof.openings.b + *beta * *zeta * k1 + gamma)
            * (proof.openings.c + *beta * *zeta * k2 + gamma);

        // denominator
        D -= (P::G1::generator() * proof.openings.c + self.circuit.sigma_3_comm * beta + P::G1::generator() * gamma) * alpha
            * (proof.openings.a + *beta * proof.openings.s_sigma_1 + gamma)
            * (proof.openings.b + *beta * proof.openings.s_sigma_2 + gamma) * proof.openings.z_shifted;

        // perm start
        D += (proof.commitments.z - P::G1::generator() * P::ScalarField::one()) * self.lagrange_1.evaluate(zeta) * alpha.square();

        // vanishing * big_t
        D -= (
            proof.commitments.t_lo
                + proof.commitments.t_mid * zeta.pow([self.domain.len() as u64])
                + proof.commitments.t_hi * zeta.pow([self.domain.len() as u64 * 2])
        ) * self.Zh.evaluate(zeta);

        D
    }
}

fn derive_challenges<P: Pairing>(
    omega: P::ScalarField,
    proof: &Proof<P>,
    public_input_vector: &[P::ScalarField],
    output_vector: &[P::ScalarField],
) -> DerivedChallenges<P::ScalarField> {
    let mut transcript = TranscriptProtocol::<P>::new(
        omega,
        [public_input_vector.clone(), output_vector.clone()].concat().as_slice(),
    );

    transcript.append_abc_commitments(proof.commitments.a.into_affine(), proof.commitments.b.into_affine(), proof.commitments.c.into_affine());
    transcript.append_z_commitment(proof.commitments.z.into_affine());
    transcript.append_t(proof.commitments.t_lo.into_affine(), proof.commitments.t_mid.into_affine(), proof.commitments.t_hi.into_affine());
    transcript.append_openings(&proof.openings);
    transcript.append_opening_proofs(proof.opening_proofs.w_zeta.into_affine(), proof.opening_proofs.w_zeta_omega.into_affine());

    let (beta, gamma) = transcript.get_beta_gamma();
    let alpha = transcript.get_alpha();
    let zeta = transcript.get_zeta();
    let vi = transcript.get_vi();
    let u = transcript.get_u();

    DerivedChallenges {
        beta,
        gamma,
        alpha,
        zeta,
        vi,
        u
    }
}
