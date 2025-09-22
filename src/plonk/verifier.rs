use ark_ec::pairing::Pairing;
use ark_ec::{CurveGroup, PrimeGroup};
use ark_ff::Field;
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;
use ark_std::{One, Zero};
use crate::kzg::{BatchOpening, MultipointOpening, KZG};
use crate::plonk::circuit::{CompiledCircuit, PublicInput};
use crate::evaluation_domain::MultiplicativeSubgroup;
use crate::plonk::proof::Proof;
use crate::plonk::prover::{pick_coset_shifters};
use crate::plonk::transcript_protocol::TranscriptProtocol;
use crate::poly_utils::{const_poly};

pub fn verify<P: Pairing>(
    circuit: &CompiledCircuit<P::ScalarField>,
    public_input: &PublicInput<P::ScalarField>,
    domain: &MultiplicativeSubgroup<P::ScalarField>,
    kzg: &KZG<P>,
    proof: Proof<P>,
) {
    let omega = domain.generator();

    let (beta, gamma, alpha, zeta, vi, u) = derive_challenges(
        omega,
        &proof,
        &public_input.pi_vector,
    );
    // let (beta, gamma) = generate_beta_gamma();
    // let alpha = generate_alpha::<P::ScalarField>();
    // let zeta = generate_zeta();
    // let vi = generate_vi::<P::ScalarField>();
    // let u = generate_u::<P::ScalarField>();

    let vanishing_poly = domain.get_vanishing_polynomial();

    let lagrange_1 = domain.lagrange_polys().first().unwrap();
    let lagrange_eval = lagrange_1.evaluate(&zeta);

    let (k1, k2) = pick_coset_shifters(&domain);

    // let solution = circuit.get_solution(&domain, k1, k2);
    // let pi = solution.pi.clone();

    let pi_zeta = public_input.pi.evaluate(&zeta);

    let _r0 = pi_zeta
        + alpha.square() * lagrange_eval
        - alpha
            * (proof.openings.a + beta * proof.openings.s_sigma_1 + gamma)
            * (proof.openings.b + beta * proof.openings.s_sigma_2 + gamma)
            * (proof.openings.c + gamma)
            * proof.openings.z_shifted;

    let D = compute_linearized_commitment(
        &kzg,
        &proof,
        circuit,
        &public_input.pi,
        beta,
        gamma,
        zeta,
        alpha,
        k1,
        k2,
        DensePolynomial::from(vanishing_poly),
        domain.len(),
        lagrange_1,
        u,
    );

    let is_valid = kzg.check_multipoint(
        &[
            MultipointOpening {
                opening_point: &zeta,
                commitments: &[D, proof.commitments.a, proof.commitments.b, proof.commitments.c, kzg.commit(&circuit.s_sigma_1), kzg.commit(&circuit.s_sigma_2)],
                batch_opening: &BatchOpening::new(
                    vec![P::ScalarField::zero(), proof.openings.a, proof.openings.b, proof.openings.c, proof.openings.s_sigma_1, proof.openings.s_sigma_2],
                    proof.opening_proofs.w_zeta,
                ),
                linearization_scalar: &vi
            },
            MultipointOpening {
                opening_point: &(zeta * omega),
                commitments: &[proof.commitments.z],
                batch_opening: &BatchOpening::new(
                    vec![proof.openings.z_shifted],
                    proof.opening_proofs.w_zeta_omega,
                ),
                linearization_scalar: &P::ScalarField::one(),
            },
        ],
        &u,
    );

    assert!(is_valid);
}

pub fn compute_linearized_commitment<P: Pairing>(
    kzg: &KZG<P>,
    proof: &Proof<P>,
    circuit: &CompiledCircuit<P::ScalarField>,
    public_input_poly: &DensePolynomial<P::ScalarField>,
    beta: P::ScalarField,
    gamma: P::ScalarField,
    zeta: P::ScalarField,
    alpha: P::ScalarField,
    k1: P::ScalarField,
    k2: P::ScalarField,
    Zh: DensePolynomial<P::ScalarField>,
    n: usize,
    lagrange_1: &DensePolynomial<P::ScalarField>,
    u: P::ScalarField,
) -> P::G1 {
    let mut D = kzg.commit(&circuit.qm) * proof.openings.a * proof.openings.b
        + kzg.commit(&circuit.ql) * proof.openings.a
        + kzg.commit(&circuit.qr) * proof.openings.b
        + kzg.commit(&circuit.qo) * proof.openings.c
        + kzg.commit(&const_poly(public_input_poly.evaluate(&zeta)))
        + kzg.commit(&circuit.qc);

    D += proof.commitments.z * alpha * (proof.openings.a + beta * zeta + gamma)
        * (proof.openings.b + beta * zeta * k1 + gamma)
        * (proof.openings.c + beta * zeta * k2 + gamma);

    D -= (P::G1::generator() * proof.openings.c + kzg.commit(&circuit.s_sigma_3) * beta + P::G1::generator() * gamma) * alpha
        * (proof.openings.a + beta * proof.openings.s_sigma_1 + gamma)
        * (proof.openings.b + beta * proof.openings.s_sigma_2 + gamma) * proof.openings.z_shifted;

    D += (proof.commitments.z - P::G1::generator() * P::ScalarField::one()) * lagrange_1.evaluate(&zeta) * alpha.square();

    D -= (proof.commitments.t_lo + proof.commitments.t_mid * zeta.pow([n as u64]) + proof.commitments.t_hi * zeta.pow([n as u64 * 2])) * Zh.evaluate(&zeta);

    D
}

fn derive_challenges<P: Pairing>(
    omega: P::ScalarField,
    proof: &Proof<P>,
    public_input_vector: &[P::ScalarField],
) -> (P::ScalarField, P::ScalarField, P::ScalarField, P::ScalarField, P::ScalarField, P::ScalarField) {
    let mut transcript = TranscriptProtocol::<P>::new(omega, public_input_vector);

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

    (beta, gamma, alpha, zeta, vi, u)
}

// pub fn compute_linearized_commitment_2<P: Pairing>(
//     kzg: &KZG<P>,
//     proof: &Proof<P>,
//     circuit: &CompiledCircuit<P::ScalarField>,
//     beta: P::ScalarField,
//     gamma: P::ScalarField,
//     zeta: P::ScalarField,
//     alpha: P::ScalarField,
//     k1: P::ScalarField,
//     k2: P::ScalarField,
//     Zh: DensePolynomial<P::ScalarField>,
//     n: usize,
//     lagrange_1: &DensePolynomial<P::ScalarField>,
//     u: P::ScalarField,
// ) -> P::G1 {
//     let gates = kzg.commit(&solution.qm) * proof.openings.a * proof.openings.b
//         + kzg.commit(&solution.ql) * proof.openings.a
//         + kzg.commit(&solution.qr) * proof.openings.b
//         + kzg.commit(&solution.qo) * proof.openings.c
//         + kzg.commit(&solution.qc);
//
//     let numerator = proof.commitments.z
//         * (
//             (proof.openings.a + beta * zeta + gamma)
//             * (proof.openings.b + beta * zeta * k1 + gamma)
//             * (proof.openings.c + beta * k2 * zeta + gamma)
//             * alpha
//             + lagrange_1.evaluate(&zeta) * alpha.square()
//             + u
//         );
//     let denominator = kzg.commit(&solution.s_sigma_3)
//         * (proof.openings.a + beta * proof.openings.s_sigma_1 + gamma)
//         * (proof.openings.b + beta * proof.openings.s_sigma_2 + gamma)
//         * alpha
//         * beta
//         * proof.openings.z_shifted;
//
//     let vanishing_t = (
//         proof.commitments.t_lo
//             + proof.commitments.t_mid * zeta.pow([n as u64])
//             + proof.commitments.t_hi * zeta.pow([2 * n as u64])
//     ) * Zh.evaluate(&zeta);
//
//     gates + (numerator - denominator) - vanishing_t
// }