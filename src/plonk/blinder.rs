use ark_ff::{FftField, Field, PrimeField};
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use crate::plonk::big_quotient_poly::BigQuotientPoly;
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

pub fn blind_big_quotient<F: PrimeField>(
    quotient_poly: BigQuotientPoly<F>,
    n: usize,
) -> BigQuotientPoly<F> {
    let [b10, b11] = generate_t_blinders();

    let (t_lo, t_mid, t_hi) = quotient_poly.into_splitted_polys();

    BigQuotientPoly::create_from_splitted_parts(
        t_lo + DensePolynomial::from(SparsePolynomial::from_coefficients_slice(&[(n, b10)])),
        t_mid + DensePolynomial::from(SparsePolynomial::from_coefficients_slice(&[(0, -b10), (n, b11)])),
        t_hi - DensePolynomial::from_coefficients_slice(&[b11]),
        n,
    )
}

// fixed values for simplicity
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

#[cfg(test)]
mod tests {
    use ark_ff::Field;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use ark_poly::univariate::DensePolynomial;
    use ark_test_curves::bls12_381::Fr;
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::plonk::big_quotient_poly::BigQuotientPoly;
    use crate::plonk::blinder::blind_big_quotient;
    use crate::poly_utils::{split_poly, to_f};

    #[test]
    fn test_t_blinding() {
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 3 }, Fr>();
        let poly = DensePolynomial::from_coefficients_vec(to_f::<Fr>(vec![
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        ]));
        let (lo, mid, hi) = split_poly(&poly, 4);
        let quotient_poly = BigQuotientPoly::create_from_splitted_parts(lo, mid, hi, 4);
        let (lo, mid, hi) = blind_big_quotient(quotient_poly, 4).into_splitted_polys();

        for w in &domain {
            assert_eq!(
                poly.evaluate(&w),
                lo.evaluate(&w) + w.pow([4]) * mid.evaluate(&w) + w.pow([8]) * hi.evaluate(&w),
            );
        }
    }
}