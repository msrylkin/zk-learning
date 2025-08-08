use std::fmt::Debug;
use std::ops::{Add, Deref};
use ark_ff::{Field, PrimeField};
use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension, Polynomial};
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};

pub fn get_reversed_vars_poly<F: Field>(poly: &DenseMultilinearExtension<F>) -> DenseMultilinearExtension<F> {
    let n = poly.num_vars();
    let mut res_poly = poly.clone();


    if n == 0 || n == 1 {
        return res_poly;
    }

    for i in 0..n / 2 {
        res_poly.relabel_in_place(i, n - i - 1,1);
    }

    res_poly
}

pub fn get_bits(k: usize, bit_len: usize) -> Vec<usize> {
    let mut result = vec![];
    for i in (0..bit_len).rev() {
        let kk = k >> i;
        let res = kk & 1;
        result.push(res);
    }

    result
}

pub fn reverse_bits(mut n: usize, k: usize) -> usize {
    let mut res = 0;
    for _ in 0..k {
        res <<= 1;
        res |= n & 1;
        n >>= 1;
    }

    res
}

// f(x_1, x_2) | l(t) = f_l(x_1(t), x_2(t))
pub fn restrict_poly<F: Field>(
    line: &Vec<DensePolynomial<F>>,
    w: DenseMultilinearExtension<F>,
) -> DensePolynomial<F> {
    let num_vars = w.num_vars();
    assert_eq!(line.len(), num_vars);

    let line_rev = line.clone().into_iter().rev().collect::<Vec<_>>();

    let mut res = DensePolynomial::from_coefficients_slice(&[F::zero()]);

    let mut evals_to_map = w.to_evaluations();
    remap_to_reverse_bits_indexing(&mut evals_to_map, num_vars);

    for (i, term) in evals_to_map.into_iter().enumerate() {
        let mut restricted = DensePolynomial::from_coefficients_slice(&[term]);
        for k in 0..num_vars {
            let bit = (i >> k) & 1;

            if bit == 0 {
                let one_poly = DensePolynomial::from_coefficients_slice(&[F::one()]);
                let mul_res = restricted.naive_mul(&(&one_poly - &line_rev[k]));

                restricted = mul_res;
            } else {
                let mul_res = restricted.naive_mul(&line_rev[k]);

                restricted = mul_res;
            }
        }

        res = res.add(restricted);
    }

    res
}

pub fn pad_with_zeroes<F: Field>(evaluations: &[F]) -> Vec<F> {
    let n = (evaluations.len() as f64).log2().ceil() as usize;
    let padded_len = 1 << n;
    let mut evaluations = Vec::from(evaluations);
    evaluations.resize(padded_len, F::zero());

    evaluations
}

pub fn interpolate<F: Field>(evaluations: &[F]) -> DenseMultilinearExtension<F> {
    let mut evaluations = pad_with_zeroes(evaluations);
    let bitlen = evaluations.len().ilog2() as usize;
    // TODO: maybe get rid of this and create solution object directly with reversed bits order
    remap_to_reverse_bits_indexing(&mut evaluations, bitlen);

    DenseMultilinearExtension::from_evaluations_slice(
        bitlen,
        &evaluations,
    )
}

pub fn remap_to_reverse_bits_indexing<T: Debug>(vec: &mut [T], k: usize) {
    for i in 0..vec.len() {
        let i_reversed = reverse_bits(i, k);
        if i_reversed > i {
            vec.swap(i, i_reversed);
        }
    }
}

// l(t) = (1 - t) * b + t * c = b + t * (c - b)
pub fn line<F: Field>(b: &[F], c: &[F]) -> Vec<DensePolynomial<F>> {
    assert_eq!(b.len(), c.len());

    let mut polys = vec![];
    for (b, c) in b.iter().zip(c.iter()) {
        polys.push(DensePolynomial::from_coefficients_slice(&[*b, *c - *b]));
    }

    polys
}

pub fn to_two_or_one_degree<F: Field>(
    poly: &DenseMultilinearExtension<F>,
    mask: usize,
    nvars: usize,
) -> DensePolynomial<F> {
    interpolate_or_const(&get_evaluations_by_mask(poly, mask, nvars))
}

pub fn interpolate_or_const<F: Field>(evals: &[F]) -> DensePolynomial<F> {
    match evals.len() {
        1 => DensePolynomial::from_coefficients_slice(evals),
        2 => interpolate_univariate_on_evals(evals.try_into().unwrap()),
        _ => panic!("wrong evaluations len"),
    }
}

// f(x) = (1 - x) * a + x * b ==> a - x * a + x * b ==> a + x * (b - a)
pub fn interpolate_univariate_on_evals<F: Field>(
    evals: &[F; 2]
) -> DensePolynomial<F> {
    let a = evals[0];
    let b = evals[1];

    DensePolynomial::from_coefficients_vec(vec![a, b - a])
}

pub fn get_evaluations_by_mask<F: Field>(
    poly: &DenseMultilinearExtension<F>,
    mask: usize,
    nvars: usize,
) -> Vec<F> {
    if nvars == poly.num_vars() {
        let bits = get_bits(mask, nvars).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
        return vec![Polynomial::evaluate(poly, &bits)];
    }

    if nvars != poly.num_vars() - 1 {
        todo!("only single variable supported for now");
    }

    let bit_set = 1 << nvars;
    let bit_set_negated = !bit_set;
    let mut mle_evals = poly.to_evaluations();
    remap_to_reverse_bits_indexing(&mut mle_evals, poly.num_vars());

    let index_1 = bit_set_negated & mask;
    let index_2 = bit_set | mask;

    let mut new_evals = vec![];

    new_evals.push(mle_evals[index_1]);
    new_evals.push(mle_evals[index_2]);

    new_evals
}

#[derive(Debug)]
struct MultiplicativeSubgroup<F: PrimeField>(Vec<F>, F);

impl<F: PrimeField> Deref for MultiplicativeSubgroup<F> {
    type Target = Vec<F>;
    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<F: PrimeField> MultiplicativeSubgroup<F> {
    pub fn new(vec: Vec<F>, generator: F) -> Self {
        Self(vec, generator)
    }

    pub fn order(&self) -> usize {
        self.0.len()
    }

    pub fn get_vanishing_polynomial(&self) -> SparsePolynomial<F> {
        SparsePolynomial::from_coefficients_vec(vec![(0, F::from(-1)), (self.order(), F::one())])
    }

    pub fn iter(&self) -> std::slice::Iter<F> {
        self.0.iter()
    }
}

impl<'a, F: PrimeField> IntoIterator for &'a MultiplicativeSubgroup<F> {
    type Item = &'a F;
    type IntoIter = std::slice::Iter<'a, F>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

fn generate_multiplicative_subgroup<const N: u64, F: PrimeField>() -> MultiplicativeSubgroup<F> {
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

    MultiplicativeSubgroup::new(H, h)
}

pub fn to_f<F: Field>(vals: Vec<i64>) -> Vec<F> {
    vals.into_iter().map(|e| F::from(e)).collect()
}

#[cfg(test)]
mod tests {
    use ark_ff::{One, Zero};
    use super::*;
    use ark_poly::DenseMultilinearExtension;
    use ark_test_curves::bls12_381::Fr;

    fn test_3_var_poly() -> DenseMultilinearExtension<Fr> {
        let nvars = 3;
        let mut evals = [1,2,3,4,5,6,7,8].into_iter().map(|e| Fr::from(e as u64)).collect::<Vec<_>>();
        remap_to_reverse_bits_indexing(&mut evals, nvars);
        DenseMultilinearExtension::from_evaluations_vec(nvars, evals)
    }

    #[test]
    fn get_evals_by_mask_one_3_var() {
        let poly = test_3_var_poly();
        let nvars = 3;
        let cases = vec![(0b010, 3), (0b011, 4), (0b110, 7)];

        for (mask, expected_val) in cases {
            let res = get_evaluations_by_mask(&poly, mask, nvars);

            assert_eq!(res.len(), 1);
            assert_eq!(res[0], Fr::from(expected_val));
        }
    }

    #[test]
    fn get_evals_by_mask_two_3_var() {
        let poly = test_3_var_poly();
        let nvars = 3;

        let cases = vec![
            (0b010, (3, 7)),
            (0b011, (4, 8)),
            (0b000, (1, 5)),
        ];

        for (mask, (coef_1, coef_2)) in cases {
            let res = get_evaluations_by_mask(&poly, mask, nvars - 1);

            assert_eq!(res.len(), 2);
            assert_eq!(res, [Fr::from(coef_1), Fr::from(coef_2)]);
        }
    }

    #[test]
    fn to_two_degree() {
        let poly = test_3_var_poly();
        let nvars = 3;

        let cases = vec![
            (0b010, (3, 4)),
            (0b011, (4, 4)),
            (0b000, (1, 4)),
        ];

        for (mask, (coef_1, coef_2)) in cases {
            let res = to_two_or_one_degree(&poly, mask, nvars - 1);

            assert_eq!(res.coeffs.len(), 2);
            assert_eq!(res.coeffs, [Fr::from(coef_1), Fr::from(coef_2)]);
        }
    }

    #[test]
    fn restrict_to_line_test() {
        let b = &[Fr::from(3), Fr::from(7)];
        let c = &[Fr::from(4), Fr::from(5)];
        let l = line(b, c);
        let w = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![10, 20, 200, 300].into_iter().map(Fr::from).collect(),
        );

        let w_restricted = restrict_poly(&l, w.clone());
        let r_star = Fr::from(12);

        let l_evals = l.iter().map(|li| li.evaluate(&r_star)).collect::<Vec<_>>();

        assert_eq!(l_evals.len(), 2);
        assert_eq!(l_evals[0], Fr::from(15));
        assert_eq!(l_evals[1], Fr::from(-17));

        let w_r_star = ark_poly::Polynomial::evaluate(&w, &l_evals);
        assert_eq!(w_r_star, Fr::from(-26020));

        let w_restricted_r_star = w_restricted.evaluate(&r_star);
        assert_eq!(w_restricted_r_star, w_r_star);
    }

    #[test]
    fn line_test() {
        let b = &[Fr::from(3), Fr::from(7)];
        let c = &[Fr::from(4), Fr::from(5)];
        let l = line(b, c);

        assert_eq!(l.len(), 2);
        assert_eq!(l[0].coeffs, [Fr::from(3), Fr::one()]);
        assert_eq!(l[1].coeffs, [Fr::from(7), Fr::from(-2)]);
    }

    #[test]
    fn test_interpolate_padded() {
        let evals = to_f::<Fr>(vec![0,1,2,3,4,5]);
        
        let poly = interpolate(&evals);
        
        assert_eq!(poly.num_vars(), 3);
        
        assert_eq!(poly.evaluate(&to_f(vec![0,0,0])), Fr::from(0));
        assert_eq!(poly.evaluate(&to_f(vec![0,0,1])), Fr::from(1));
        assert_eq!(poly.evaluate(&to_f(vec![0,1,0])), Fr::from(2));
        assert_eq!(poly.evaluate(&to_f(vec![0,1,1])), Fr::from(3));
        assert_eq!(poly.evaluate(&to_f(vec![1,0,0])), Fr::from(4));
        assert_eq!(poly.evaluate(&to_f(vec![1,0,1])), Fr::from(5));
        assert_eq!(poly.evaluate(&to_f(vec![1,1,0])), Fr::from(0));
        assert_eq!(poly.evaluate(&to_f(vec![1,1,1])), Fr::from(0));
    }

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
}