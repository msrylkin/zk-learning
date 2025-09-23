use ark_ff::PrimeField;
use ark_poly::DenseUVPolynomial;
use ark_poly::univariate::DensePolynomial;
use crate::poly_utils::split_poly;

pub struct BigQuotientPoly<F: PrimeField> {
    t_lo: DensePolynomial<F>,
    t_mid: DensePolynomial<F>,
    t_hi: DensePolynomial<F>,
    n: usize,
}

pub fn shift_coefficients<F: PrimeField>(mut values: Vec<F>, n: usize) -> Vec<F> {
    let orig_len = values.len();
    values.resize(n + values.len(), F::zero());

    for i in 0..orig_len {
        values[i + n] = values[i];
    }

    for i in 0..n {
        values[i] = F::zero();
    }

    values
}

impl<F: PrimeField> BigQuotientPoly<F> {
    pub fn create_for_domain(t: DensePolynomial<F>, n: usize) -> Self {
        println!("orig t {:?}", t);
        let (t_lo, t_mid, t_hi) = split_poly(&t, n);

        Self {
            // t,
            t_lo,
            t_mid,
            t_hi,
            n
        }
    }

    pub fn create_from_splitted_parts(
        t_lo: DensePolynomial<F>,
        t_mid: DensePolynomial<F>,
        t_hi: DensePolynomial<F>,
        n: usize,
    ) -> Self {
        Self {
            t_lo,
            t_mid,
            t_hi,
            n
        }
    }

    pub fn get_splitted_polys(&self) -> (&DensePolynomial<F>, &DensePolynomial<F>, &DensePolynomial<F>) {
        (&self.t_lo, &self.t_mid, &self.t_hi)
    }

    pub fn to_splitted_polys(self) -> (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>) {
        (self.t_lo, self.t_mid, self.t_hi)
    }

    pub fn linearize_on_shifts(&self, x_shift: &F) -> DensePolynomial<F> {
        &self.t_lo
            + &self.t_mid * x_shift.pow([self.n as u64])
            + &self.t_hi * x_shift.pow([self.n as u64 * 2])
    }

    pub fn to_combined(self) -> DensePolynomial<F> {
        self.t_lo
            + DensePolynomial::from_coefficients_vec(shift_coefficients(self.t_mid.coeffs, self.n))
            + DensePolynomial::from_coefficients_vec(shift_coefficients(self.t_hi.coeffs, self.n * 2))
    }
}