use std::ops::{Div, Mul, Sub};
use ark_ec::{AffineRepr, CurveGroup};
use ark_ec::pairing::Pairing;
use ark_std::{UniformRand, Zero};
use ark_test_curves::bls12_381::{Bls12_381, Fr, G1Affine as G1, G2Affine as G2};
use ark_ff::{One};
use ark_poly::univariate::{DensePolynomial};
use ark_ec::short_weierstrass::{Affine};
use ark_ff::Field;
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_test_curves::bls12_381::g1::Config as G1Config;
use ark_test_curves::bls12_381::g2::Config as G2Config;

pub fn run_kzg() {
    let polynomial = DensePolynomial::from_coefficients_vec(
        vec![Fr::from(11), Fr::from(12), Fr::from(345)]
    );

    let (srs, tau_g2) = setup();
    let commitment = commit_to_poly(&polynomial, &srs);

    let random_a = Fr::from(777);

    let (f_a, evaluation_proof) = open_commitment(&srs, &random_a, &polynomial);

    let is_valid = verify_commitment(&commitment, &tau_g2, &random_a, &f_a, &evaluation_proof);

    println!("is_valid = {}", is_valid);
}

fn setup() -> (Vec<Affine<G1Config>>, Affine<G2Config>) {
    let max_degree = 8;
    let rng = &mut ark_std::test_rng();

    let tau = Fr::rand(rng);

    let tau_g2 = (G2::generator() * tau).into_affine();
    let srs = (0..max_degree).map(|d| {
        let tau_times_d  = tau.pow(&[d]);
        G1::generator().mul(tau_times_d).into_affine()
    }).collect::<Vec<_>>();

    (srs, tau_g2)
}

fn commit_to_poly(
    poly: &DensePolynomial<Fr>,
    srs: &Vec<Affine<G1Config>>
) -> Affine<G1Config> {
    let comm_f = eval_poly_on_srs(poly, srs);

    comm_f
}

fn open_commitment(
    srs: &Vec<Affine<G1Config>>,
    a: &Fr,
    poly: &DensePolynomial<Fr>,
) -> (Fr, Affine<G1Config>) {
    let f_a = poly.evaluate(a);

    let f_a_poly = DensePolynomial::from_coefficients_vec(vec![f_a]);
    let x_minus_a_poly = DensePolynomial::from_coefficients_vec(vec![-*a, Fr::one()]);
    let q_x = poly.sub(&f_a_poly).div(&x_minus_a_poly);

    let evaluation_proof = eval_poly_on_srs(&q_x, srs);

    (f_a, evaluation_proof)
}

fn verify_commitment(
    commitment: &Affine<G1Config>,
    tau_g2: &Affine<G2Config>,
    a: &Fr,
    f_a: &Fr,
    evaluation_proof: &Affine<G1Config>,
) -> bool {
    let l_pairing = Bls12_381::pairing(
        evaluation_proof,
        *tau_g2 - G2::generator() * a,
    );
    let r_pairing = Bls12_381::pairing(
        *commitment - G1::generator() * f_a,
        G2::generator(),
    );

    l_pairing == r_pairing
}

fn eval_poly_on_srs(
    poly: &DensePolynomial<Fr>,
    srs: &[Affine<G1Config>]
) -> Affine<G1Config> {
    if srs.len() < poly.degree() {
        panic!("srs lower than poly degree");
    }

    poly.coeffs
        .iter()
        .zip(srs)
        .map(|(x_i, h)| h.mul(x_i))
        .reduce(|a, b| a + b)
        .unwrap()
        .into_affine()
}