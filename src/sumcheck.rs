use std::fmt::Debug;
use ark_ff::{Field, One};
use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_std::Zero;
use ark_test_curves::bls12_381::Fr;

pub fn test_sc() {
    let mut poly = DenseMultilinearExtension::from_evaluations_vec(
        4,
        vec![91, 62, 13, 431, 98, 123, 2871, 7, 512, 12, 63 ,982, 474, 2847, 912, 744].into_iter().map(Fr::from).collect()
    );

    // prove(poly);
    prove_2(poly);

    println!("poly proved");
}

impl<F: Field> SumCheckPoly<F> for DenseMultilinearExtension<F> {
    fn get_evaluations(&self) -> Vec<F> {
        self.to_evaluations()
    }

    fn num_vars(&self) -> usize {
        self.num_vars
    }

    fn get_partial_sum_poly(&self) -> DensePolynomial<F> {
        let sum_vars = SumCheckPoly::num_vars(self) - 1;
        let mut current_poly = self.clone();
    
        for i in 0..sum_vars {
            let summed_poly = sum_over_last_variable(&current_poly);
            current_poly = summed_poly;
        }
    
        interpolate_univariate_on_evals(&current_poly.evaluations.try_into().unwrap())
    }

    // fn get_partial_sum_poly(&self) -> DensePolynomial<F> {
    //     let mut current_sc_poly = self.clone();
    //     let num_vars = MultilinearExtension::num_vars(&current_sc_poly);
    // 
    //     let evals = self.get_evaluations();
    //     println!("evals {:?}", evals);
    // 
    //     let (zeroes, ones) = evals.iter().enumerate().fold((F::zero(), F::zero()), |(zeroes, ones), (i, eval)| {
    //         if i % 2 == 0 {
    //             (zeroes + eval, ones)
    //         } else {
    //             (zeroes, ones + eval)
    //         }
    //     });
    // 
    //     interpolate_univariate_on_evals(&[zeroes, ones])
    // }

    fn fix_variables(&self, partial_point: &[F]) -> Self {
        MultilinearExtension::fix_variables(self, partial_point)
    }

    fn evaluate(&self, point: &[F]) -> F {
        Polynomial::evaluate(self, &point.to_vec())
    }
}

pub trait SumCheckPoly<F: Field> {
    fn get_evaluations(&self) -> Vec<F>;

    fn num_vars(&self) -> usize;

    fn get_partial_sum_poly(&self) -> DensePolynomial<F>;

    fn fix_variables(&self, partial_point: &[F]) -> Self;

    fn evaluate(&self, point: &[F]) -> F;
}

fn prove_for_val<F: Field, S: SumCheckPoly<F> + Clone>(poly: S) {
    let mut current_poly = poly.clone();
    // let mut current_poly = poly;
    // TODO: Fiat-Shamir
    let random_points = vec![21,23,35,48].into_iter()
        .map(F::from)
        .collect::<Vec<_>>();
    let mut current_sum = current_poly.get_evaluations().iter().sum::<F>();
    let mut total_rounds = 0;

    for round in 0..poly.num_vars() {
        // let partial_sum_poly = get_partial_sum_poly(&current_poly);
        let partial_sum_poly = current_poly.get_partial_sum_poly();

        if partial_sum_poly.degree() != 1 {
            panic!("degree does not match");
        }

        let eval_0 = partial_sum_poly.evaluate(&F::zero());
        let eval_1 = partial_sum_poly.evaluate(&F::one());

        let eval_sum = eval_0 + eval_1;

        if eval_sum != current_sum {
            panic!("sum does not match");
        }

        current_poly = current_poly.fix_variables(&[random_points[round]]);
        current_sum = current_poly.get_evaluations().iter().sum::<F>();
        total_rounds += 1;
    }

    let last_evaluation = current_poly.get_evaluations().pop().unwrap();
    let points = random_points.into_iter().take(total_rounds).collect::<Vec<_>>();
    let original_poly_evaluation = poly.evaluate(&points);

    if last_evaluation != original_poly_evaluation {
        panic!("last evaluation not match");
    }
    // println!("")
}

pub fn prove_2<F: Field, S: SumCheckPoly<F> + Clone + Debug>(
    poly: S
) -> (Vec<F>, Vec<DensePolynomial<F>>) {
    let random_points = vec![32,23,35,48].into_iter()
        .map(F::from)
        .collect::<Vec<_>>();
    let H = poly.get_evaluations().iter().sum::<F>();
    println!("poly {:?} vnum {}", poly.get_evaluations(), poly.num_vars());
    println!("H {}", H);
    let num_vars = poly.num_vars();

    // round 1
    let partial_sum_poly = poly.get_partial_sum_poly();
    let g1_0 = partial_sum_poly.evaluate(&F::zero());
    let g1_1 = partial_sum_poly.evaluate(&F::one());

    println!("round 0\npartial_sum_poly {:?}\ng1_0 {}\ng1_1 {}", partial_sum_poly, g1_0, g1_1);

    // assert_eq!(partial_sum_poly.degree(), 2);
    assert_eq!(H, g1_0 + g1_1);

    let mut r_vals = vec![];
    let mut partial_sum_polys = vec![];
    
    r_vals.push(random_points[0].clone());
    
    let mut current_poly = poly.clone();
    current_poly = current_poly.fix_variables(&r_vals);
    
    let mut previous_partial_sum_poly = partial_sum_poly.clone();
    partial_sum_polys.push(previous_partial_sum_poly.clone());

    for i in 1..(num_vars - 1) {
        let partial_sum_poly = current_poly.get_partial_sum_poly();
        let g1_0 = partial_sum_poly.evaluate(&F::zero());
        let g1_1 = partial_sum_poly.evaluate(&F::one());
        println!("round {} partial sum poly {:?}", i, partial_sum_poly);
        println!("round {} r_vals {:?} g1_0 {} g1_1 {} g1_sum {}", i, r_vals, g1_0, g1_1, g1_0 + g1_1);
        let r = r_vals.last().unwrap().clone();
        let previous_poly_at_r = previous_partial_sum_poly.evaluate(&r);

        // assert_eq!(partial_sum_poly.degree(), 1);
        assert_eq!(g1_1 + g1_0, previous_poly_at_r);

        previous_partial_sum_poly = partial_sum_poly.clone();
        partial_sum_polys.push(previous_partial_sum_poly.clone());
        r_vals.push(random_points[i].clone());
        current_poly = current_poly.fix_variables(&[random_points[i].clone()]);
    }

    let last_partial_sum_poly = current_poly.get_partial_sum_poly();
    
    println!("last partial sum poly {:?}", last_partial_sum_poly);
    let gv_0 = last_partial_sum_poly.evaluate(&F::zero());
    let gv_1 = last_partial_sum_poly.evaluate(&F::one());
    println!("gv_0 {}", gv_0);
    println!("gv_1 {}", gv_1);
    let g_v_1_r = previous_partial_sum_poly.evaluate(r_vals.last().unwrap());
    println!("previous_partial_sum_poly {:?}", previous_partial_sum_poly);
    
    assert_eq!(gv_0 + gv_1, g_v_1_r);
    
    // -12719300 52435875175126190479447740508185965837690552500527637822603658699938568465213
    // -24282300 52435875175126190479447740508185965837690552500527637822603658699938556902213

    let last_r = F::from(100);
    r_vals.push(last_r);
    println!("r_vals {:?}", r_vals);
    previous_partial_sum_poly = last_partial_sum_poly.clone();
    partial_sum_polys.push(last_partial_sum_poly.clone());

    let last_eval = poly.evaluate(&r_vals.clone());

    println!("partial_sum_polys {:?}", partial_sum_polys);
    println!("current_poly varsnum {}", current_poly.num_vars());
    println!("las partial sum poly {:?}", last_partial_sum_poly);
    println!("las partial sum poly eval {}", last_partial_sum_poly.evaluate(&last_r));

    assert_eq!(previous_partial_sum_poly.evaluate(&last_r), last_eval);

    (r_vals, partial_sum_polys)
}

pub fn prove<F: Field, S: SumCheckPoly<F> + Clone>(
    // poly: DenseMultilinearExtension<Fr>,
    poly: S,
) -> (Vec<F>, Vec<DensePolynomial<F>>) {
    let mut current_poly = poly.clone();
    // let mut current_poly = poly;
    // TODO: Fiat-Shamir
    let random_points = vec![21,23,35,48].into_iter()
        .map(F::from)
        .collect::<Vec<_>>();
    let mut current_sum = current_poly.get_evaluations().iter().sum::<F>();
    let mut total_rounds = 0;

    let mut used_r = vec![];
    let mut used_polys = vec![];

    for round in 0..poly.num_vars() {
        // let partial_sum_poly = get_partial_sum_poly(&current_poly);
        let partial_sum_poly = current_poly.get_partial_sum_poly();
        println!("partial_sum_poly {:?}\n", partial_sum_poly);
        used_polys.push(partial_sum_poly.clone());

        if partial_sum_poly.degree() != 1 {
            panic!("degree does not match");
        }

        let eval_0 = partial_sum_poly.evaluate(&F::zero());
        let eval_1 = partial_sum_poly.evaluate(&F::one());

        let eval_sum = eval_0 + eval_1;

        if eval_sum != current_sum {
            panic!("sum does not match");
        }

        let r = random_points[round];

        current_poly = current_poly.fix_variables(&[r]);
        current_sum = current_poly.get_evaluations().iter().sum::<F>();
        println!("current_sum {}", current_sum);
        total_rounds += 1;
        used_r.push(r);
    }

    let last_evaluation = current_poly.get_evaluations().pop().unwrap();
    let points = random_points.into_iter().take(total_rounds).collect::<Vec<_>>();
    let original_poly_evaluation = poly.evaluate(&points);
    // println!("points {:?} used_r {:?}", points, used_r);

    if last_evaluation != original_poly_evaluation {
        panic!("last evaluation not match");
    }

    println!("last eval {:?}", last_evaluation);

    (used_r, used_polys)
}

pub fn get_partial_sum_poly(
    poly: &DenseMultilinearExtension<Fr>,
) -> DensePolynomial<Fr> {
    let sum_vars = SumCheckPoly::num_vars(poly) - 1;
    let mut current_poly = poly.clone();

    for i in 0..sum_vars {
        let summed_poly = sum_over_last_variable(&current_poly);
        current_poly = summed_poly;
    }

    interpolate_univariate_on_evals(&current_poly.evaluations.try_into().unwrap())
}

// 000 |
// 001 --> 00[1 + 0]
// 010 |
// 011 --> 01[1 + 0]
pub fn sum_over_last_variable<F: Field>(poly: &DenseMultilinearExtension<F>) -> DenseMultilinearExtension<F> {
    let sums = poly
        .evaluations
        .iter()
        .enumerate()
        .fold(vec![F::zero(); poly.evaluations.len() / 2], |mut sums, (i, e)| {
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
pub fn interpolate_univariate_on_evals<F: Field>(
    evals: &[F; 2]
) -> DensePolynomial<F> {
    let a = evals[0];
    let b = evals[1];

    DensePolynomial::from_coefficients_vec(vec![a, b - a])
}