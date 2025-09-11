use ark_ff::{FftField, Field, PrimeField};
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial, SparsePolynomial};
use ark_std::iterable::Iterable;
use ark_std::Zero;
use crate::plonk::blinder::{blind_solution, blind_splitted_t, blind_z_poly};
use crate::plonk::circuit::{CompiledCircuit, Solution};
use crate::plonk::permutation::PermutationArgument;
use crate::plonk::proof::Openings;
use crate::poly_utils::{generate_lagrange_basis_polys, generate_multiplicative_subgroup, interpolate_univariate, split_poly};

fn prove<F: PrimeField + FftField>(
    circuit: &CompiledCircuit<F>,
) {

    let domain = generate_multiplicative_subgroup::<16, F>();
    let omega = domain.generator();
    let Zh = domain.get_vanishing_polynomial();
    let (k1, k2) = pick_coset_shifters(&domain);
    let solution = circuit.get_solution(&domain, k1, k2);
    let solution = blind_solution(solution, &Zh);
    let lagrange_1 = &generate_lagrange_basis_polys(&domain)[0];

    let (beta, gamma) = generate_beta_gamma::<F>();

    let permutation_argument = PermutationArgument::new(&domain, &beta, &gamma, &solution);
    let z_poly = permutation_argument.z_poly();
    let z_poly = blind_z_poly(&z_poly, &Zh);

    let z_shifted = shift_poly(&z_poly, &omega);

    let alpha = generate_alpha();

    let big_t = compute_big_quotient(&solution, &z_poly, &z_shifted, &permutation_argument, &Zh, lagrange_1, &alpha);

    let (lo_t, mid_t, hi_t) = split_poly(&big_t, domain.len());
    let (lo_t, mid_t, hi_t) = blind_splitted_t(&lo_t, &mid_t, &hi_t, domain.len());

    let zeta = generate_zeta();

    let openings = compute_openings(&solution, &z_poly, &zeta);

    let r_poly = linearization_poly(
        &openings,
        &solution,
        &permutation_argument,
        &z_poly,
        &DensePolynomial::from(Zh),
        lagrange_1,
        &zeta,
        &alpha,
        &lo_t,
        &mid_t,
        &hi_t,
        domain.len(),
        &domain,
    );
}

fn compute_gate_check_poly<F: PrimeField + FftField>(solution: &Solution<F>) -> DensePolynomial<F> {
    &solution.a * &solution.b * &solution.qm
        + &solution.a * &solution.ql
        + &solution.b * &solution.qr
        + &solution.c * &solution.qo
        + &solution.pi + &solution.qc
}

pub fn generate_beta_gamma<F: FftField + PrimeField>() -> (F, F) {
    (F::from(43), F::from(35))
}

pub fn generate_alpha<F: FftField + PrimeField>() -> F {
    F::from(4004)
}

pub fn generate_vi<F: FftField + PrimeField>() -> F {
    F::from(4234)
}

pub fn generate_zeta<F: FftField + PrimeField>() -> F {
    F::from(776655)
}


pub fn pick_coset_shifters<F: Field>(domain: &[F]) -> (F, F) {
    let mut i = 2;
    let n = domain.len();
    let (k1, k2) = loop {
        let k1 = F::from(i);
        let k2 = k1.square();

        let k1_n = k1.pow(&[n as u64]);
        let k2_n  = k2.pow(&[n as u64]);
        let k1_over_k2 = (k1 / k2).pow(&[n as u64]);

        if !k1_n.is_one() && !k2_n.is_one() && !k1_over_k2.is_one() {
            break (k1, k2);
        }

        i += 1;
    };

    (k1, k2)
}

fn divide_by_vanishing<F: FftField + PrimeField>(poly: &DensePolynomial<F>, Zh: &SparsePolynomial<F>) -> DensePolynomial<F> {
    let res = DenseOrSparsePolynomial::from(poly).divide_with_q_and_r(&DenseOrSparsePolynomial::from(Zh));

    let (q, r) = res.unwrap();

    if !r.is_zero() {
        panic!("Remainder is not zero");
    }

    q
}

fn shift_poly<F: Field>(poly: &DensePolynomial<F>, scalar: &F) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(
        poly.coeffs
            .iter()
            .enumerate()
            .map(|(i, c)| *c * scalar.pow(&[i as u64]))
            .collect(),
    )
}

fn compute_big_quotient<F: FftField + PrimeField>(
    solution: &Solution<F>,
    z: &DensePolynomial<F>,
    z_shifted: &DensePolynomial<F>,
    perm_argument: &PermutationArgument<F>,
    Zh: &SparsePolynomial<F>,
    lagrange_base_1: &DensePolynomial<F>,
    alpha: &F,
) -> DensePolynomial<F> {
    let gate_check_poly = compute_gate_check_poly(&solution);

    let perm_numerator_poly = perm_argument.numerator_poly() * z;
    let perm_denominator_poly = perm_argument.denominator_poly() * z_shifted;

    let z_poly_m1 = (z - DensePolynomial::from_coefficients_slice(&[F::one()])) * lagrange_base_1;

    divide_by_vanishing(&gate_check_poly, Zh)
        + divide_by_vanishing(&(perm_numerator_poly - perm_denominator_poly), Zh) * *alpha
        + divide_by_vanishing(&z_poly_m1, Zh) * alpha.square()
}

fn compute_openings<F: FftField + PrimeField>(
    solution: &Solution<F>,
    z_shifted: &DensePolynomial<F>,
    zeta: &F,
) -> Openings<F> {
    Openings {
        a: solution.a.evaluate(zeta),
        b: solution.b.evaluate(zeta),
        c: solution.c.evaluate(zeta),
        s_sigma_1: solution.s_sigma_1.evaluate(zeta),
        s_sigma_2: solution.s_sigma_2.evaluate(zeta),
        z_shifted: z_shifted.evaluate(zeta),
    }
}

fn linearization_poly<F: FftField + PrimeField>(
    openings: &Openings<F>,
    solution: &Solution<F>,
    permutation_argument: &PermutationArgument<F>,
    z_poly: &DensePolynomial<F>,
    Zh: &DensePolynomial<F>,
    lagrange_1: &DensePolynomial<F>,
    zeta: &F,
    alpha: &F,
    t_lo: &DensePolynomial<F>,
    t_mid: &DensePolynomial<F>,
    t_hi: &DensePolynomial<F>,
    n: usize,
    domain: &[F],
) -> DensePolynomial<F> {
    let gate_check_poly_linearized = &solution.qm * openings.a * openings.b
        + &solution.ql * openings.a
        + &solution.qr * openings.b
        + &solution.qo * openings.c
        + const_poly(solution.pi.evaluate(&zeta)) + &solution.qc;

    let perm_numerator_poly_linearized = z_poly * permutation_argument.numerator_poly().evaluate(&zeta);

    let sigma_a_linearized = openings.a + permutation_argument.beta() * openings.s_sigma_1 + permutation_argument.gamma();
    let sigma_b_linearized = openings.b + permutation_argument.beta() * openings.s_sigma_2 + permutation_argument.gamma();
    let sigma_c_lin_poly = const_poly(openings.c)
        + &solution.s_sigma_3 * permutation_argument.beta()
        + const_poly(permutation_argument.gamma());
    let perm_denominator_poly_linearized = sigma_c_lin_poly * sigma_a_linearized * sigma_b_linearized * openings.z_shifted;

    let perm_start_linearized = (z_poly - const_poly(F::one())) * lagrange_1.evaluate(zeta);

    let linearized_vanishing_t = (
        t_lo
        + t_mid * zeta.pow([n as u64])
        + t_hi * zeta.pow([n as u64 * 2])
    ) * Zh.evaluate(zeta);

    let r_poly = gate_check_poly_linearized
        + perm_numerator_poly_linearized * *alpha
        - perm_denominator_poly_linearized * *alpha
        + perm_start_linearized * alpha.square()
        - linearized_vanishing_t;

    r_poly
}

fn const_poly<F: Field>(c: F) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_slice(&[c])
}

#[cfg(test)]
mod tests {
    use ark_ff::Field;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use ark_poly::univariate::DensePolynomial;
    use ark_std::{One, Zero};
    use ark_test_curves::bls12_381::Fr;
    use crate::kzg::{setup, KZG};
    use crate::plonk::blinder::blind_splitted_t;
    use crate::plonk::circuit::get_test_circuit;
    use crate::plonk::permutation::PermutationArgument;
    use crate::plonk::prover::{compute_big_quotient, compute_gate_check_poly, compute_openings, interpolate_univariate, linearization_poly, pick_coset_shifters, shift_poly};
    use crate::poly_utils::{generate_lagrange_basis_polys, generate_multiplicative_subgroup, split_poly, to_f};

    #[test]
    fn pick_coset_shifters_test() {
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 6 }, Fr>();
        let (k1, k2) = pick_coset_shifters(&domain);

        let coset_k1 = domain.iter().map(|e| k1 * e).collect::<Vec<_>>();
        let coset_k2 = domain.iter().map(|e| k2 * e).collect::<Vec<_>>();

        for h in &domain {
            for k1h in &coset_k1 {
                for k2h in &coset_k2 {
                    assert_ne!(h, k1h);
                    assert_ne!(k1h, k2h);
                }
            }
        }
    }

    #[test]
    fn test_z_poly() {
        let test_circuit = get_test_circuit();
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 3 }, Fr>();
        let (k1, k2) = pick_coset_shifters(&domain);
        let (a, b, c) = test_circuit.get_abc_vectors();
        let (a, b, c) = (
            interpolate_univariate(&domain, &a),
            interpolate_univariate(&domain, &b),
            interpolate_univariate(&domain, &c),
        );

        let beta = Fr::from(43);
        let gamma = Fr::from(35);

        let solution = test_circuit.get_solution(&domain, k1, k2);
        let permutation = PermutationArgument::new(&domain, &beta, &gamma, &solution);

        let z_poly = permutation.z_poly();

        assert_eq!(z_poly.evaluate(&domain[domain.len() - 1]), Fr::one());
        assert_eq!(z_poly.evaluate(&Fr::one()), Fr::one());

        let omega= domain.generator();
        let solution = test_circuit.get_solution(&domain, k1, k2);
        let permutation = PermutationArgument::new(&domain, &beta, &gamma, &solution);
        let num_poly = permutation.numerator_poly();
        let denom_poly = permutation.denominator_poly();

        for x in &domain {
            assert_eq!(
                z_poly.evaluate(&(omega * x)) * denom_poly.evaluate(x),
                z_poly.evaluate(x) * num_poly.evaluate(x),
            );
        }
    }

    #[test]
    fn compute_big_quotient_test() {
        let test_circuit = get_test_circuit();
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 3 }, Fr>();
        let Zh = domain.get_vanishing_polynomial();
        let (k1, k2) = pick_coset_shifters(&domain);
        let beta = Fr::from(43);
        let gamma = Fr::from(35);

        let solution = test_circuit.get_solution(&domain, k1, k2);
        let perm_argument = PermutationArgument::new(&domain, &beta, &gamma, &solution);
        let z = perm_argument.z_poly();
        let z_shifted = shift_poly(&z, &domain.generator());
        let alpha = Fr::from(123);

        let big_q = compute_big_quotient(
            &solution,
            &z,
            &z_shifted,
            &perm_argument,
            &Zh,
            &generate_lagrange_basis_polys(&domain)[0],
            &alpha,
        );

        let gate_check_poly = compute_gate_check_poly(&solution);

        let z_shifted = shift_poly(&z, &domain.generator());
        let perm_numerator_poly = perm_argument.numerator_poly() * z;
        let z = perm_argument.z_poly();
        let perm_denominator_poly = perm_argument.denominator_poly() * z_shifted;
        let lagrange_base_1 = &generate_lagrange_basis_polys(&domain)[0];
        let z_poly_m1 = (z - DensePolynomial::from_coefficients_slice(&[Fr::one()])) * lagrange_base_1;
        let alpha = Fr::from(123);

        assert_eq!(
            big_q * DensePolynomial::from(Zh),
            gate_check_poly
                + (perm_numerator_poly - perm_denominator_poly) * alpha
                + z_poly_m1 * alpha.square()
        );
    }

    #[test]
    fn test_t_blinding() {
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 3 }, Fr>();
        let poly = DensePolynomial::from_coefficients_vec(to_f::<Fr>(vec![
            1, 2, 3, 4, 5, 6, 7, 8, 9, 10
        ]));
        let (lo, mid, hi) = split_poly(&poly, 4);
        let (lo, mid, hi) = blind_splitted_t(&lo, &mid, &hi, 4);

        for w in &domain {
            assert_eq!(
                poly.evaluate(&w),
                lo.evaluate(&w) + w.pow([4]) * mid.evaluate(&w) + w.pow([8]) * hi.evaluate(&w),
            );
        }
    }

    #[test]
    fn test_linearization_poly() {
        let test_circuit = get_test_circuit();
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 3 }, Fr>();
        let omega = domain.generator();
        let Zh = domain.get_vanishing_polynomial();
        let (k1, k2) = pick_coset_shifters(&domain);
        let beta = Fr::from(43);
        let gamma = Fr::from(35);

        let solution = test_circuit.get_solution(&domain, k1, k2);
        let perm_argument = PermutationArgument::new(&domain, &beta, &gamma, &solution);
        let z = perm_argument.z_poly();
        let z_shifted = shift_poly(&z, &domain.generator());
        let alpha = Fr::from(123);
        let zeta = Fr::from(999);

        let lagrange_base_1 = &generate_lagrange_basis_polys(&domain)[0];

        let big_q = compute_big_quotient(
            &solution,
            &z,
            &z_shifted,
            &perm_argument,
            &Zh,
            &lagrange_base_1,
            &alpha,
        );
        let (lo, mid, hi) = split_poly(&big_q, domain.len());

        let openings = compute_openings(&solution, &z_shifted, &zeta);
        let r_poly = linearization_poly(
            &openings,
            &solution,
            &perm_argument,
            &z,
            &DensePolynomial::from(Zh),
            &lagrange_base_1,
            &zeta,
            &alpha,
            &lo,
            &mid,
            &hi,
            domain.len(),
            &domain,
        );

        assert_eq!(r_poly.evaluate(&zeta), Fr::zero());
    }
}