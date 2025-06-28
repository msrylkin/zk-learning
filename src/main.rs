mod kzg;
mod schnorr;
mod gkr;
mod sumcheck;

use std::ops::{Div, Mul, Sub};
use std::process::id;
use ark_ec::{AdditiveGroup, PrimeGroup, AffineRepr, CurveGroup};
use ark_ec::pairing::Pairing;
use ark_std::{Zero, UniformRand};
use ark_test_curves::bls12_381::{Bls12_381, G1Affine as G1, G2Affine as G2, Fr, Fq};
use ark_ff::{FftField, Fp256, One};
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial};
use ark_ec::short_weierstrass::{Affine, SWCurveConfig};
use ark_ff::Field;
use ark_poly::Polynomial;
use ark_test_curves::bls12_381::g1::Config as G1Config;
use ark_test_curves::bls12_381::g2::Config as G2Config;

use ark_poly::DenseUVPolynomial;

fn run_kzg() {
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

fn main() {
    // run_kzg();
    // schnorr::run_schnorr();
    // sumcheck::test_sc();
    gkr::test_gkr();
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

fn multi_open() {
    
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

fn kzg() {
    let rng = &mut ark_std::test_rng();

    let tau = Fr::rand(rng);

    let tau_pow0 = tau.pow(&[0]);
    let tau_pow1 = tau.pow(&[1]);
    let tau_pow2 = tau.pow(&[2]);
    // 
    // println!("tau {}\ntau 0 {}\ntau 1 {}\ntau 2 {}\n", tau, tau_pow0, tau_pow1, tau_pow2);
    // println!("gen {}\ngen * tau 0 {}\n gent * tau 1 {}\n gen * tau 2 {}\n",
    //          G1::generator(),
    //          G1::generator() * tau_pow0,
    //          G1::generator() * tau_pow1,
    //          G1::generator() * tau_pow2
    // );

    let f_x = DensePolynomial {
        coeffs: vec![Fr::from(5), Fr::from(7), Fr::from(17)],
    };
    


    let srs = (0..8).map(|d| {
        let tau_times_d  = tau.pow(&[d]);
        G1::generator().mul(tau_times_d).into_affine()
    }).collect::<Vec<_>>();

    let tau_g2 = G2::generator() * tau;

    // println!("srs {:?}", srs);

    let comm_f = f_x.coeffs.iter().zip(&srs).map(|(x_i, h)| {
        h.mul(x_i)
    }).reduce(|a, b| a + b).unwrap();

    let a = Fr::from(777);
    let f_a = f_x.evaluate(&a);
    
    let f_a_poly = DensePolynomial {
        coeffs: vec![Fr::from(f_a)],
    };
    let x_minus_f_a_poly = DensePolynomial {
        coeffs: vec![-a, Fr::one()],
    };

    let q_x = f_x.clone().sub(&f_a_poly).div(&x_minus_f_a_poly);

    // let evaluation_proof = q_x.coeffs.iter().zip(&srs).map(|(x_i, h)| {
    //     h.mul(x_i)
    // }).reduce(|a, b| a + b).unwrap();
    let evaluation_proof = eval_poly_on_srs(&q_x, &srs);

    let l_pairing = Bls12_381::pairing(
        evaluation_proof,
        tau_g2 - G2::generator() * a
    );
    let r_pairing = Bls12_381::pairing(
        comm_f - G1::generator() * f_a,
        G2::generator()
    );
    
    // println!("l {:?}", l_pairing);
    // println!("r {:?}", r_pairing);

    // println!("comm_f, {}\nf_a {}\nq_x {:?}\n eval_proof {}", comm_f, f_a, q_x, evaluation_proof);

    let (q, r) = DenseOrSparsePolynomial::divide_with_q_and_r(
        &f_x.clone().sub(&f_a_poly).clone().into(),
        &x_minus_f_a_poly.clone().into()
    ).unwrap();

    println!("\nf_x {:?}\n", f_x);
    println!("f_a_poly {:?}\n", f_a_poly);
    println!("x - a poly {:?}", x_minus_f_a_poly);
    println!("\nf_x - f_a_poly {:?}", f_x.clone().sub(&f_a_poly));
    println!("\nqoutient poly {:?}", (&f_x).sub(&f_a_poly).div(&x_minus_f_a_poly));
    // println!("f_x {:?}\nq_x sub f_a {:?}\nf_a - 5 {}", f_x, f_x.clone().sub(&DensePolynomial {
    //     coeffs: vec![Fr::from(f_a)],
    // }), -f_a + Fr::from(5));

    println!("l == r {}", l_pairing == r_pairing);
    
    println!("\nq{:?}\n\nr {:?}", q, r)
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
