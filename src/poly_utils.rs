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