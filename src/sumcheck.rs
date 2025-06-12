use ark_ff::One;
use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_std::Zero;
use ark_test_curves::bls12_381::Fr;

pub fn test_sc() {
    let poly = DenseMultilinearExtension::from_evaluations_vec(
        3,
        vec![91, 62, 13, 431, 98, 123, 2871, 7].into_iter().map(Fr::from).collect()
    );

    prove(
        poly,
    );
    
    println!("poly proved");
}

fn prove(
    poly: DenseMultilinearExtension<Fr>,
) {
    let mut current_poly = poly.clone();
    let random_points = vec![21,23,35,48].into_iter()
        .map(Fr::from)
        .collect::<Vec<_>>();
    let mut current_sum = current_poly.evaluations.iter().sum::<Fr>();
    let mut total_rounds = 0;

    for round in 0..poly.num_vars() {
        let partial_sum_poly = get_partial_sum_poly(&current_poly);

        if partial_sum_poly.degree() != 1 {
            panic!("degree does not match");
        }

        let eval_0 = partial_sum_poly.evaluate(&Fr::zero());
        let eval_1 = partial_sum_poly.evaluate(&Fr::one());

        let eval_sum = eval_0 + eval_1;

        if eval_sum != current_sum {
            panic!("sum does not match");
        }

        current_poly = current_poly.fix_variables(&[random_points[round]]);
        current_sum = current_poly.evaluations.iter().sum::<Fr>();
        total_rounds += 1;
    }

    let last_evaluation = current_poly.evaluations.pop().unwrap();
    let points = random_points.into_iter().take(total_rounds).collect::<Vec<_>>();
    let original_poly_evaluation = poly.evaluate(&points);

    if last_evaluation != original_poly_evaluation {
        panic!("last evaluation not match");
    }
}

pub fn get_partial_sum_poly(
    poly: &DenseMultilinearExtension<Fr>,
) -> DensePolynomial<Fr> {
    let sum_vars = poly.num_vars() - 1;
    let mut current_poly = poly.clone();

    for i in 0..sum_vars {
        let summed_poly = sum_over_last_variable(&current_poly);
        current_poly = summed_poly;
    }

    interpolate_univariate_on_evals(&current_poly.evaluations.try_into().unwrap())
}

fn sum_over_last_variable(poly: &DenseMultilinearExtension<Fr>) -> DenseMultilinearExtension<Fr> {
    let sums = poly
        .evaluations
        .iter()
        .enumerate()
        .fold(vec![Fr::zero(); poly.evaluations.len() / 2], |mut sums, (i, e)| {
            if i % 2 == 0 {
                sums[(i / 4) * 2] += e;
            } else {
                sums[(i / 4) * 2 + 1] += e;
            }
            sums
        });

    let res = DenseMultilinearExtension::from_evaluations_vec(poly.num_vars - 1, sums);

    res
}

// f(x) = (1 - x) * a + x * b ==> a - x * a + x * b ==> a + x * (b - a)
fn interpolate_univariate_on_evals(
    evals: &[Fr; 2]
) -> DensePolynomial<Fr> {
    let a = evals[0];
    let b = evals[1];

    DensePolynomial::from_coefficients_vec(vec![a, b - a])
}