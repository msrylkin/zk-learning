use ark_ff::{FftField, Field, PrimeField};
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use crate::plonk::circuit::solution::Solution;

pub fn blind_solution<F: PrimeField + FftField>(
    solution: Solution<F>,
    Zh: &SparsePolynomial<F>
)  -> Solution<F> {
    let [b1, b2, b3, b4, b5, b6] = generate_solution_blinders::<F>();

    Solution {
        a: blind_poly(&solution.a, Zh, &[b1, b2]),
        b: blind_poly(&solution.b, Zh, &[b3, b4]),
        c: blind_poly(&solution.c, Zh, &[b5, b6]),
        ..solution
    }
}

pub fn blind_z_poly<F: PrimeField + FftField>(z: &DensePolynomial<F>, Zh: &SparsePolynomial<F>) -> DensePolynomial<F> {
    blind_poly(z, Zh, &generate_z_blinders())
}

pub fn blind_splitted_t<F: FftField + PrimeField>(
    t_lo: &DensePolynomial<F>,
    t_mid: &DensePolynomial<F>,
    l_hi: &DensePolynomial<F>,
    n: usize,
) -> (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>) {
    let [b10, b11] = generate_t_blinders();

    (
        t_lo + DensePolynomial::from(SparsePolynomial::from_coefficients_slice(&[(n, b10)])),
        t_mid + DensePolynomial::from(SparsePolynomial::from_coefficients_slice(&[(0, -b10), (n, b11)])),
        l_hi - DensePolynomial::from_coefficients_slice(&[b11]),
    )
}

fn generate_solution_blinders<F: Field>() -> [F; 6] {
    [
        F::from(11),
        F::from(22),
        F::from(33),
        F::from(44),
        F::from(55),
        F::from(66),
    ]
}

fn generate_z_blinders<F: Field>() -> [F; 3] {
    [
        F::from(77),
        F::from(88),
        F::from(99),
    ]
}

fn generate_t_blinders<F: Field>() -> [F; 2] {
    [
        F::from(110),
        F::from(111),
    ]
}

fn blind_poly<F: FftField>(
    poly: &DensePolynomial<F>,
    Zh: &SparsePolynomial<F>,
    blinders: &[F],
) -> DensePolynomial<F> {
    poly + DensePolynomial::from(Zh.clone()) * DensePolynomial::from_coefficients_slice(blinders)
}