use std::hash::{DefaultHasher, Hash, Hasher};
use ark_ec::AffineRepr;
use ark_ff::Field;
use ark_test_curves::bls12_381::{G1Affine as G1, Fr, G1Projective};

pub fn run_schnorr() {
    run_exploit_fs_weak();
}

pub fn run_exploit_fs_weak() {
    let R = G1::generator() * Fr::from(909090);
    let siglet = Fr::from(8989892);
    let challenge = challenge_fs(R);
    let X = (G1::generator() * siglet - R) * challenge.inverse().unwrap();
    let is_valid = verifier_verify(siglet, X, R, challenge);

    println!("{}", is_valid);
}

pub fn run_schnorr_fs_weak() {
    // step 1
    let (x, X) = prover_data();
    let (r, R) = prover_blinder();

    // step 2
    let challenge = challenge_fs(R);

    // step 3
    let siglet = prover_prove(x, r, challenge);

    // step 4
    let is_valid = verifier_verify(siglet, X, R, challenge);

    println!("{}", is_valid);
}

pub fn run_schnorr_interactive() {
    // step 1
    let (x, X) = prover_data();
    let (r, R) = prover_blinder();

    // step 2
    let challenge = verifier_challenge();

    // step 3
    let siglet = prover_prove(x, r, challenge);

    // step 4
    let is_valid = verifier_verify(siglet, X, R, challenge);

    println!("{}", is_valid);
}

fn prover_data() -> (Fr, G1Projective) {
    let x = Fr::from(789);
    let X = G1::generator() * x;

    (x, X)
}

fn prover_blinder() -> (Fr, G1Projective) {
    let r = Fr::from(998877);
    let R = G1::generator() * r;

    (r, R)
}

fn prover_prove(
    x: Fr,
    r: Fr,
    challenge: Fr,
) -> Fr {
    x * challenge + r
}

fn prover_fake_proof(
    r: Fr,
    challenge: Fr,
) -> Fr {
    r + challenge
}

fn verifier_challenge() -> Fr {
    Fr::from(1122334458)
}

fn verifier_verify(
    siglet: Fr,
    X: G1Projective,
    R: G1Projective,
    challenge: Fr,
) -> bool {
    G1::generator() * siglet == X * challenge + R
}

fn challenge_fs(R: G1Projective) -> Fr {
    let mut hasher = DefaultHasher::new();
    R.x.hash(&mut hasher);
    Fr::from(hasher.finish())
}