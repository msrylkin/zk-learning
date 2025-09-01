use ark_ff::{FftField, Field, PrimeField};
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial, SparsePolynomial};
use ark_std::iterable::Iterable;
use ark_std::Zero;
use ark_test_curves::bls12_381::Fr;
use crate::plonk::circuit::CompiledCircuit;
use crate::plonk::permutation::PermutationArgument;
use crate::poly_utils::generate_multiplicative_subgroup;

fn prove<F: FftField + PrimeField>(
    A: &[F],
    B: &[F],
    C: &[F],
    ql: DensePolynomial<F>,
    qr: DensePolynomial<F>,
    qm: DensePolynomial<F>,
    qo: DensePolynomial<F>,
    qc: DensePolynomial<F>,
) {
    let [b1, b2, b3, b4, b5, b6, b7, b8, b9] = generate_blinders::<F>();
    let domain = generate_multiplicative_subgroup::<16, F>();
    let Zh = domain.get_vanishing_polynomial();

    let a = interpolate_univariate(&domain, A);
    let b = interpolate_univariate(&domain, A);
    let c = interpolate_univariate(&domain, A);

    let a_blinded = blind_poly(&a, &Zh, &[b1, b2]);
    let b_blinded = blind_poly(&b, &Zh, &[b3, b4]);
    let c_blinded = blind_poly(&c, &Zh, &[b5, b6]);
}

fn blind_poly<F: FftField>(
    poly: &DensePolynomial<F>,
    Zh: &SparsePolynomial<F>,
    blinders: &[F],
) -> DensePolynomial<F> {
    poly + DensePolynomial::from(Zh.clone()) * DensePolynomial::from_coefficients_slice(blinders)
}

fn generate_blinders<F: Field>() -> [F; 9] {
    [
        F::from(11),
        F::from(22),
        F::from(33),
        F::from(44),
        F::from(55),
        F::from(66),
        F::from(77),
        F::from(88),
        F::from(99),
    ]
}

pub fn interpolate_univariate<F: Field>(domain: &[F], values : &[F]) -> DensePolynomial<F> {
    assert!(domain.len() >= values.len());

    generate_lagrange_basis_polys(&domain[0..values.len()])
        .into_iter()
        .zip(values)
        .map(|(lp, y)| lp * *y)
        .fold(DensePolynomial::zero(), |acc, lp| acc + lp)
}

fn generate_lagrange_basis_polys<F: Field>(domain: &[F]) -> Vec<DensePolynomial<F>> {
    let mut lagrange_polys = vec![];

    for (i, x_i) in domain.iter().enumerate() {
        // calculate denominator
        let mut full_denom = F::one();

        for (j, x_j) in domain.iter().enumerate() {
            if j != i {
                full_denom *= (*x_i - *x_j).inverse().unwrap();
            }
        }

        let remapped = domain.iter()
            .enumerate()
            .filter_map(|(j, e)| {
                if i != j {
                    Some(*e)
                } else {
                    None
                }
            }).collect::<Vec<_>>();

        let sums = symmetric_sums(&remapped, remapped.len());

        let mut res = vec![];
        let mut sign = F::one();

        for degree_index in 0..sums.len() {
            res.push(full_denom * sums[degree_index] * sign);

            sign *= F::from(-1);
        }

        res.reverse();

        lagrange_polys.push(DensePolynomial::from_coefficients_vec(res));
    }

    lagrange_polys
}

pub fn pick_coset_shifters<F: Field>(domain: &[F]) -> (F, F) {
    let mut i = 2;
    let n = domain.len();
    let (k1, k2) = loop {
        let k1 = F::from(i);
        let k2 = k1.square();

        let k1_n = k1.pow(&[n as u64]);
        let k2_n  = k2.pow(&[n as u64]);
        let k1_over_k2 = (k1 * k2.inverse().unwrap()).pow(&[n as u64]);

        if !k1_n.is_one() && !k2_n.is_one() && !k1_over_k2.is_one() {
            break (k1, k2);
        }

        i += 1;
    };

    (k1, k2)
}

fn symmetric_sums<F: Field>(values: &[F], max_k: usize) -> Vec<F> {
    assert!(values.len() >= max_k);

    let mut dp = vec![vec![F::zero(); max_k + 1]; values.len() + 1];

    dp[0][0] = F::one();

    for x in &mut dp {
        x[0] = F::one();
    }

    for j in 1..dp.len() {
        let x = values[j - 1];
        for k in 1..max_k + 1 {
            dp[j][k] = x * dp[j - 1][k - 1] + dp[j - 1][k];
        }
    }

    dp.remove(dp.len() - 1)
}

// fn permutation_product_acc<F: FftField + PrimeField>(
//     a: &DensePolynomial<F>,
//     b: &DensePolynomial<F>,
//     c: &DensePolynomial<F>,
//     perm_poly_1: &DensePolynomial<F>,
//     perm_poly_2: &DensePolynomial<F>,
//     perm_poly_3: &DensePolynomial<F>,
//     beta: &F,
//     gamma: &F,
//     domain: &[F],
//     max_i: usize,
// ) -> F {
//     domain
//         .iter()
//         .take(max_i)
//         .map(|w|permutation_product(a, b, c, perm_poly_1, perm_poly_2, perm_poly_3, beta, gamma, &w))
//         .product()
// }

// fn permutation_product<F: FftField + PrimeField>(
//     a: &DensePolynomial<F>,
//     b: &DensePolynomial<F>,
//     c: &DensePolynomial<F>,
//     perm_poly_1: &DensePolynomial<F>,
//     perm_poly_2: &DensePolynomial<F>,
//     perm_poly_3: &DensePolynomial<F>,
//     beta: &F,
//     gamma: &F,
//     point: &F,
// ) -> F {
//     hash_permutation(a, &perm_poly_1, point, &beta, &gamma)
//     * hash_permutation(b, &perm_poly_2, point, &beta, &gamma)
//     * hash_permutation(c, &perm_poly_3, point, &beta, &gamma)
// }

// fn numerator_poly<F: FftField + PrimeField>(
//     a: &DensePolynomial<F>,
//     b: &DensePolynomial<F>,
//     c: &DensePolynomial<F>,
//     compiled_circuit: &CompiledCircuit<F>,
//     beta: &F,
//     gamma: &F,
//     domain: &[F],
//     k1: F,
//     k2: F,
// ) -> DensePolynomial<F> {
//     let (sid_1, sid_2, sid_3) = compiled_circuit.get_s_id_polys(domain, k1, k2);
//     let values = domain.iter().map(|w| {
//         permutation_product(a, b, c, &sid_1, &sid_2, &sid_3, &beta, &gamma, &w)
//     }).collect::<Vec<_>>();
//
//     interpolate_univariate(domain, &values)
// }

// fn denominator_poly<F: FftField + PrimeField>(
//     a: &DensePolynomial<F>,
//     b: &DensePolynomial<F>,
//     c: &DensePolynomial<F>,
//     compiled_circuit: &CompiledCircuit<F>,
//     beta: &F,
//     gamma: &F,
//     domain: &[F],
//     k1: F,
//     k2: F,
// ) -> DensePolynomial<F> {
//     let (sigma_1, sigma_2, sigma_3) = compiled_circuit.get_sigma_polys(domain, k1, k2);
//     let values = domain.iter().map(|w| {
//         permutation_product(a, b, c, &sigma_1, &sigma_2, &sigma_3, &beta, &gamma, &w)
//     }).collect::<Vec<_>>();
//
//     interpolate_univariate(domain, &values)
// }

// fn hash_permutation<F: PrimeField + FftField>(
//     f: &DensePolynomial<F>,
//     perm_poly: &DensePolynomial<F>,
//     point: &F,
//     beta: &F,
//     gamma: &F,
// ) -> F {
//     f.evaluate(&point) + *beta * perm_poly.evaluate(&point) + *gamma
// }
//
// fn numerator_acc<F: FftField + PrimeField>(
//     a: &DensePolynomial<F>,
//     b: &DensePolynomial<F>,
//     c: &DensePolynomial<F>,
//     compiled_circuit: &CompiledCircuit<F>,
//     beta: &F,
//     gamma: &F,
//     domain: &[F],
//     k1: F,
//     k2: F,
//     max_i: usize,
// ) -> F {
//     let (sid_1, sid_2, sid_3) = compiled_circuit.get_s_id_polys(domain, k1, k2);
//
//     permutation_product_acc(
//         a,
//         b,
//         c,
//         &sid_1,
//         &sid_2,
//         &sid_3,
//         beta,
//         gamma,
//         domain,
//         max_i,
//     )
// }
//
// fn denominator_acc<F: FftField + PrimeField>(
//     a: &DensePolynomial<F>,
//     b: &DensePolynomial<F>,
//     c: &DensePolynomial<F>,
//     compiled_circuit: &CompiledCircuit<F>,
//     beta: &F,
//     gamma: &F,
//     domain: &[F],
//     k1: F,
//     k2: F,
//     max_i: usize,
// ) -> F {
//     let (s_sigma_1, s_sigma_2, s_sigma_3) = compiled_circuit.get_sigma_polys(domain, k1, k2);
//
//     permutation_product_acc(
//         a,
//         b,
//         c,
//         &s_sigma_1,
//         &s_sigma_2,
//         &s_sigma_3,
//         beta,
//         gamma,
//         domain,
//         max_i,
//     )
// }

fn divide_by_vanishing<F: FftField + PrimeField>(poly: &DensePolynomial<F>, Zh: &SparsePolynomial<F>) -> DensePolynomial<F> {
    let res = DenseOrSparsePolynomial::from(poly).divide_with_q_and_r(&DenseOrSparsePolynomial::from(Zh));

    let (q, r) = res.unwrap();

    if !r.is_zero() {
        panic!("Remainder is not zero");
    }

    q
}

fn z_poly<F: FftField + PrimeField>(
    a: &DensePolynomial<F>,
    b: &DensePolynomial<F>,
    c: &DensePolynomial<F>,
    compiled_circuit: &CompiledCircuit<F>,
    beta: &F,
    gamma: &F,
    domain: &[F],
    k1: F,
    k2: F,
) -> DensePolynomial<F> {
    let solution = compiled_circuit.get_solution(domain, k1, k2);
    let permutation_argument = PermutationArgument::new(domain, beta, gamma, &solution);
    let values = (0..domain.len())
        .map(|i| {
            let num = permutation_argument.numerator_acc(i);
            let denom = permutation_argument.denominator_acc(i);
            println!("i {i} num / denom {}", num / denom);

            num / denom
        })
        .collect::<Vec<_>>();

    interpolate_univariate(domain, &values)
}

fn shift_poly<F: Field>(poly: &DensePolynomial<F>, scalar: F) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(
        poly.coeffs
            .iter()
            .enumerate()
            .map(|(i, c)| *c * scalar.pow(&[i as u64]))
            .collect(),
    )
}

fn compute_big_quotient<F: FftField + PrimeField>(
    a: &DensePolynomial<F>,
    b: &DensePolynomial<F>,
    c: &DensePolynomial<F>,
    z: &DensePolynomial<F>,
    public_input_poly: &DensePolynomial<F>,
    beta: &F,
    gamma: &F,
    k1: F,
    k2: F,
    Zh: &SparsePolynomial<F>,
    circuit: &CompiledCircuit<F>,
    domain: &[F],
    omega: F,
) -> DensePolynomial<F> {
    let (ql, qr, qm , qo, qc) = circuit.get_selectors(domain);
    let gate_check_poly = a * b * qm + a * ql + b * qr + c * qo + public_input_poly + qc;
    let (ql, qr, qm , qo, qc) = circuit.get_selectors(domain);
    let z_shifted = shift_poly(&z, omega);

    let solution = circuit.get_solution(domain, k1, k2);
    let perm_argument = PermutationArgument::new(domain, beta, gamma, &solution);

    let perm_numerator_poly = perm_argument.numerator_poly() * z;
    let perm_denominator_poly = perm_argument.denominator_poly() * z_shifted;
    let lagrange_base_1 = &generate_lagrange_basis_polys(domain)[0];
    let z_poly_m1 = (z - DensePolynomial::from_coefficients_slice(&[F::one()])) * lagrange_base_1;
    let alpha = F::from(123);

     for w in domain {
         println!("\na {} b {} c {}", a.evaluate(&w), b.evaluate(&w), c.evaluate(&w));
         println!(
             "ql {} qr {} qm {} qo {} qc {} pi {}",
             ql.evaluate(&w),
             qr.evaluate(&w),
             qm.evaluate(&w),
             qo.evaluate(&w),
             qc.evaluate(&w),
             public_input_poly.evaluate(&w),
         );
         println!("eval {} | c - {}", gate_check_poly.evaluate(&w), c.evaluate(&w));
     }

    for w in domain {
        println!("\nperm num {}", perm_numerator_poly.evaluate(&w));
        println!("perm denom {}", perm_denominator_poly.evaluate(&w));
    }
    let z_shifted = z * domain[1];


    for w in domain {
        println!("\nz {}", z.evaluate(&w));
        println!("z shifted {}", z_shifted.evaluate(&w));
    }

    divide_by_vanishing(&gate_check_poly, Zh)
        + divide_by_vanishing(&(perm_numerator_poly - perm_denominator_poly), Zh) * alpha
        + divide_by_vanishing(&z_poly_m1, Zh) * alpha.square()
}

#[cfg(test)]
mod tests {
    use ark_poly::{Polynomial};
    use ark_std::One;
    use ark_test_curves::bls12_381::Fr;
    use crate::plonk::circuit::get_test_circuit;
    use crate::plonk::permutation::PermutationArgument;
    use crate::plonk::prover::{compute_big_quotient, interpolate_univariate, pick_coset_shifters, symmetric_sums, z_poly};
    use crate::poly_utils::{generate_multiplicative_subgroup, to_f};

    #[test]
    fn test_symmetric_sums_4() {
        let sums = symmetric_sums(&[
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(4),
        ], 4);

        assert_eq!(sums, to_f::<Fr>(vec![1, 10, 35, 50, 24]));
    }

    #[test]
    fn test_symmetric_sums_5() {
        let sums = symmetric_sums(&[
            Fr::from(1),
            Fr::from(2),
            Fr::from(3),
            Fr::from(4),
            Fr::from(5),
        ], 5);

        assert_eq!(sums, to_f::<Fr>(vec![1, 15, 85, 225, 274, 120]));
    }

    #[test]
    fn test_interpolate_univariate() {
        let values = to_f::<Fr>(vec![8,10,15]);
        let domain = to_f::<Fr>(vec![1,2,3]);

        let poly = interpolate_univariate(&domain, &values);

        for (y, x) in values.into_iter().zip(domain) {
            assert_eq!(poly.evaluate(&x), y);
        }
    }

    #[test]
    fn test_interpolate_univariate_mul_domain() {
        let values = to_f::<Fr>(vec![8,10,15]);
        let domain = generate_multiplicative_subgroup::<3, Fr>();

        let poly = interpolate_univariate(&domain, &values);

        for (y, x) in values.into_iter().zip(domain.into_iter()) {
            assert_eq!(poly.evaluate(&x), y);
        }
    }

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

        let z_poly = z_poly(&a, &b, &c, &test_circuit, &beta, &gamma, &domain, k1, k2);

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
        let (a, b, c) = test_circuit.get_abc_vectors();
        let (a, b, c) = (
            interpolate_univariate(&domain, &a),
            interpolate_univariate(&domain, &b),
            interpolate_univariate(&domain, &c),
        );
        let beta = Fr::from(43);
        let gamma = Fr::from(35);
        let z = z_poly(&a, &b, &c, &test_circuit, &beta, &gamma, &domain, k1, k2);
        let pi = test_circuit.get_public_input_poly(&domain);

        let big_q = compute_big_quotient(
            &a,
            &b,
            &c,
            &z,
            &pi,
            &beta,
            &gamma,
            k1,
            k2,
            &Zh,
            &test_circuit,
            &domain,
            domain.generator(),
        );
    }
}