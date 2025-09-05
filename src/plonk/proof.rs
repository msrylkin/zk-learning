use ark_ff::{FftField, PrimeField};

pub struct Openings<F: PrimeField + FftField> {
    pub a: F,
    pub b: F,
    pub c: F,
    pub s_sigma_1: F,
    pub s_sigma_2: F,
    pub z_shifted: F,
}
