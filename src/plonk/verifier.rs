use ark_ec::pairing::Pairing;
use ark_ec::PrimeGroup;
use ark_ff::Field;
use ark_poly::Polynomial;
use ark_poly::univariate::DensePolynomial;
use ark_std::{One, Zero};
use ark_test_curves::bls12_381::{Fr, G1Projective};
use crate::kzg;
use crate::kzg::{setup, BatchOpening, MultipointOpening, KZG};
use crate::plonk::circuit::{CompiledCircuit, Solution};
use crate::plonk::proof::Proof;
use crate::plonk::prover::{generate_alpha, generate_beta_gamma, generate_u, generate_vi, generate_zeta, pick_coset_shifters};
use crate::poly_utils::{const_poly, generate_lagrange_basis_polys, generate_multiplicative_subgroup, MultiplicativeSubgroup};

pub fn verify<P: Pairing>(
    circuit: &CompiledCircuit<P::ScalarField>,
    domain: &MultiplicativeSubgroup<P::ScalarField>,
    kzg: &KZG<P>,
    proof: Proof<P>,
) {
    let (beta, gamma) = generate_beta_gamma();
    let alpha = generate_alpha::<P::ScalarField>();
    let zeta = generate_zeta();
    let vi = generate_vi::<P::ScalarField>();
    let u = generate_u::<P::ScalarField>();
    let omega = domain.generator();

    let vanishing_poly = domain.get_vanishing_polynomial();

    let lagrange_1 = &generate_lagrange_basis_polys(&domain)[0];
    let lagrange_eval = lagrange_1.evaluate(&zeta);

    let (k1, k2) = pick_coset_shifters(&domain);

    let solution = circuit.get_solution(&domain, k1, k2);
    let pi = solution.pi.clone();

    let pi_zeta = pi.evaluate(&zeta);

    let r0 = pi_zeta
        + alpha.square() * lagrange_eval
        - alpha
            * (proof.openings.a + beta * proof.openings.s_sigma_1 + gamma)
            * (proof.openings.b + beta * proof.openings.s_sigma_2 + gamma)
            * (proof.openings.c + gamma)
            * proof.openings.z_shifted;

    let D = compute_linearized_commitment(
        &kzg,
        &proof,
        &solution,
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
                commitments: &[D, proof.commitments.a, proof.commitments.b, proof.commitments.c, kzg.commit(&solution.s_sigma_1), kzg.commit(&solution.s_sigma_2)],
                batch_opening: &BatchOpening::new(
                    vec![P::ScalarField::zero(), proof.openings.a, proof.openings.b, proof.openings.c, proof.openings.s_sigma_1, proof.openings.s_sigma_2],
                    proof.opening_proofs.w_zeta,
                ),
                linearization_scalar: &vi
            },
            MultipointOpening {
                opening_point: &(zeta * omega),
                batch_opening: &BatchOpening::new(
                    vec![proof.openings.z_shifted],
                    proof.opening_proofs.w_zeta_omega,
                ),
                linearization_scalar: &P::ScalarField::one(),
                commitments: &[proof.commitments.z],
            },
        ],
        &u,
    );

    assert!(is_valid);
}

pub fn compute_linearized_commitment<P: Pairing>(
    kzg: &KZG<P>,
    proof: &Proof<P>,
    solution: &Solution<P::ScalarField>,
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
    let mut D = kzg.commit(&solution.qm) * proof.openings.a * proof.openings.b
        + kzg.commit(&solution.ql) * proof.openings.a
        + kzg.commit(&solution.qr) * proof.openings.b
        + kzg.commit(&solution.qo) * proof.openings.c
        + kzg.commit(&const_poly(solution.pi.evaluate(&zeta)))
        + kzg.commit(&solution.qc);

    D += proof.commitments.z * alpha * (proof.openings.a + beta * zeta + gamma)
        * (proof.openings.b + beta * zeta * k1 + gamma)
        * (proof.openings.c + beta * zeta * k2 + gamma);

    D -= (P::G1::generator() * proof.openings.c + kzg.commit(&solution.s_sigma_3) * beta + P::G1::generator() * gamma) * alpha
        * (proof.openings.a + beta * proof.openings.s_sigma_1 + gamma)
        * (proof.openings.b + beta * proof.openings.s_sigma_2 + gamma) * proof.openings.z_shifted;

    D += (proof.commitments.z - P::G1::generator() * P::ScalarField::one()) * lagrange_1.evaluate(&zeta) * alpha.square();

    D -= (proof.commitments.t_lo + proof.commitments.t_mid * zeta.pow([n as u64]) + proof.commitments.t_hi * zeta.pow([n as u64 * 2])) * Zh.evaluate(&zeta);

    D
}

pub fn compute_linearized_commitment_2<P: Pairing>(
    kzg: &KZG<P>,
    proof: &Proof<P>,
    solution: &Solution<P::ScalarField>,
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
    let gates = kzg.commit(&solution.qm) * proof.openings.a * proof.openings.b
        + kzg.commit(&solution.ql) * proof.openings.a
        + kzg.commit(&solution.qr) * proof.openings.b
        + kzg.commit(&solution.qo) * proof.openings.c
        + kzg.commit(&solution.qc);

    let numerator = proof.commitments.z
        * (
            (proof.openings.a + beta * zeta + gamma)
            * (proof.openings.b + beta * zeta * k1 + gamma)
            * (proof.openings.c + beta * k2 * zeta + gamma)
            * alpha
            + lagrange_1.evaluate(&zeta) * alpha.square()
            + u
        );
    let denominator = kzg.commit(&solution.s_sigma_3)
        * (proof.openings.a + beta * proof.openings.s_sigma_1 + gamma)
        * (proof.openings.b + beta * proof.openings.s_sigma_2 + gamma)
        * alpha
        * beta
        * proof.openings.z_shifted;

    let vanishing_t = (
        proof.commitments.t_lo
            + proof.commitments.t_mid * zeta.pow([n as u64])
            + proof.commitments.t_hi * zeta.pow([2 * n as u64])
    ) * Zh.evaluate(&zeta);

    gates + (numerator - denominator) - vanishing_t
}