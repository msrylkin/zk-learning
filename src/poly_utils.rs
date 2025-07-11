use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};

pub fn get_reversed_vars_poly<F: Field>(poly: &DenseMultilinearExtension<F>) -> DenseMultilinearExtension<F> {
    let n = MultilinearExtension::num_vars(poly);
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