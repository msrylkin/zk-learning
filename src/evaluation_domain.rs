use std::ops::Deref;
use ark_ff::{Field, PrimeField};
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use crate::poly_utils::{generate_lagrange_basis_polys, interpolate_on_lagrange_basis_polys};

#[derive(Debug)]
pub struct MultiplicativeSubgroup<F: PrimeField> {
    subgroup: Vec<F>,
    generator: F,
    lagrange_polys: Vec<DensePolynomial<F>>,
}

impl<F: PrimeField> Deref for MultiplicativeSubgroup<F> {
    type Target = Vec<F>;
    fn deref(&self) -> &Self::Target {
        &self.subgroup
    }
}

impl<F: PrimeField> MultiplicativeSubgroup<F> {
    pub fn new(subgroup: Vec<F>, generator: F, lagrange_polys: Vec<DensePolynomial<F>>) -> Self {
        Self {
            subgroup,
            generator,
            lagrange_polys,
        }
    }

    pub fn order(&self) -> usize {
        self.subgroup.len()
    }

    pub fn get_vanishing_polynomial(&self) -> SparsePolynomial<F> {
        SparsePolynomial::from_coefficients_vec(vec![(0, F::from(-1)), (self.order(), F::one())])
    }

    pub fn iter(&self) -> std::slice::Iter<F> {
        self.subgroup.iter()
    }

    pub fn generator(&self) -> F {
        self.generator
    }

    pub fn lagrange_polys(&self) -> &[DensePolynomial<F>] {
        self.lagrange_polys.as_slice()
    }

    pub fn interpolate_univariate(&self, values: &[F]) -> DensePolynomial<F> {
        interpolate_on_lagrange_basis_polys(&self.lagrange_polys.as_slice(), values)
    }
}

impl<'a, F: PrimeField> IntoIterator for &'a MultiplicativeSubgroup<F> {
    type Item = &'a F;
    type IntoIter = std::slice::Iter<'a, F>;

    fn into_iter(self) -> Self::IntoIter {
        self.subgroup.iter()
    }
}

pub fn generate_multiplicative_subgroup<const N: u64, F: PrimeField>() -> MultiplicativeSubgroup<F> {
    let subgroup_order = N;
    let field_multiplicative_order = F::zero() - F::one();

    let cofactor = field_multiplicative_order.div(F::from(subgroup_order));
    let remainder = field_multiplicative_order - F::from(subgroup_order).mul(cofactor);

    assert!(remainder.is_zero());

    let mut ii = 3;
    let h = loop {
        let g = F::from(ii);
        let h = g.pow(cofactor.into_bigint());

        let mut has_inf = false;

        for i in 1..subgroup_order {
            if h.pow([i]) == F::one() {
                has_inf = true;
            }
        }

        if !has_inf {
            break h;
        }

        ii += 1;
    };

    let H = (0..subgroup_order)
        .into_iter()
        .map(|i| h.pow([i]))
        .collect::<Vec<_>>();

    let lagrange_polys = generate_lagrange_basis_polys(&H);

    MultiplicativeSubgroup::new(H, h, lagrange_polys)
}

#[cfg(test)]
mod tests {
    use ark_ff::Field;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use ark_poly::univariate::DensePolynomial;
    use ark_std::{One, Zero};
    use ark_test_curves::bls12_381::Fr;
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::poly_utils::{to_f};

    #[test]
    fn vanishing_poly_test() {
        const SUBGROUP_ORDER: u64 = 16;
        let subgroup = generate_multiplicative_subgroup::<{ SUBGROUP_ORDER }, Fr>();

        let Z_h = subgroup.get_vanishing_polynomial();

        for h in &subgroup {
            assert_eq!(Z_h.evaluate(h), Fr::zero());
        }

        assert_ne!(Z_h.evaluate(&Fr::from(123)), Fr::zero());
    }

    #[test]
    fn generate_multiplicative_subgroup_test() {
        const SUBGROUP_ORDER: u64 = 16;
        let subgroup = generate_multiplicative_subgroup::<{ SUBGROUP_ORDER }, Fr>();

        assert_eq!(subgroup.iter().sum::<Fr>(), Fr::zero());

        for h in &subgroup {
            assert_eq!(h.pow([SUBGROUP_ORDER]), Fr::one());
        }

        let poly = DensePolynomial::from_coefficients_vec(vec![
            Fr::from(123),
            Fr::from(22),
            Fr::from(33),
            Fr::from(456),
            Fr::from(9),
        ]);

        let mut subgroup_evals_sum = Fr::zero();

        for h in &subgroup {
            subgroup_evals_sum += poly.evaluate(h);
        }

        assert_eq!(poly.evaluate(&Fr::zero()) * Fr::from(SUBGROUP_ORDER), subgroup_evals_sum);
    }

    #[test]
    fn test_interpolate_univariate_mul_domain() {
        let values = to_f::<Fr>(vec![8,10,15]);
        let domain = generate_multiplicative_subgroup::<3, Fr>();

        let poly = domain.interpolate_univariate(&values);

        for (y, x) in values.into_iter().zip(domain.into_iter()) {
            assert_eq!(poly.evaluate(&x), y);
        }
    }
}