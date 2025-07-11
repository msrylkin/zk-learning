pub(crate) mod gkr;
mod circuit;
mod prover;
mod verifier;
mod common;
mod test_utils;

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::Fr;
    use ark_ff::{Field, Zero};
    use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
    use crate::gkr::prover::LayerRoundPoly;
    use crate::poly_utils::{get_evaluations_by_mask, remap_to_reverse_bits_indexing, to_two_or_one_degree};
    use crate::sumcheck::SumCheckPoly;
    use super::*;

    // fn test_3_var_poly() -> DenseMultilinearExtension<Fr> {
    //     let nvars = 3;
    //     let mut evals = [1,2,3,4,5,6,7,8].into_iter().map(|e| Fr::from(e as u64)).collect::<Vec<_>>();
    //     remap_to_reverse_bits_indexing(&mut evals, nvars);
    //     DenseMultilinearExtension::from_evaluations_vec(nvars, evals)
    // }
    // 
    // fn get_test_round_poly_2_vars<F: Field>() -> LayerRoundPoly<F> {
    //     let add_i = DenseMultilinearExtension::from_evaluations_vec(
    //         2,
    //         vec![0,0,0,0].into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>(),
    //     );
    //     let mul_i = DenseMultilinearExtension::from_evaluations_vec(
    //         2,
    //         vec![0,0,1,0].into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>(),
    //     );
    //     let Wi_1 = DenseMultilinearExtension::from_evaluations_vec(
    //         1,
    //         vec![210, 320].into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>(),
    //     );
    //     // LayerRoundPoly {
    //     //     add_i,
    //     //     mul_i,
    //     //     Wi_1_a: Wi_1.clone(),
    //     //     Wi_1_b: Wi_1.clone(),
    //     // }
    // 
    //     LayerRoundPoly::new(add_i, mul_i, Wi_1.clone(), Wi_1.clone())
    // }
    // 
    // fn get_test_round_poly_4_vars<F: Field>() -> LayerRoundPoly<F> {
    //     let add_i = DenseMultilinearExtension::from_evaluations_vec(
    //         5,
    //         (0..32).into_iter().map(|e| F::from((e == 16 || e == 27) as u64)).collect::<Vec<_>>(),
    //     );
    //     let mul_i = DenseMultilinearExtension::from_evaluations_vec(
    //         5,
    //         (0..32).into_iter().map(|e| F::from(0)).collect::<Vec<_>>(),
    //     );
    //     let Wi_1 = DenseMultilinearExtension::from_evaluations_vec(
    //         2,
    //         vec![10,20,200,300].into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>(),
    //     );
    // 
    //     // LayerRoundPoly {
    //     //     add_i: MultilinearExtension::fix_variables(&add_i, &[F::from(3)]),
    //     //     mul_i: MultilinearExtension::fix_variables(&mul_i, &[F::from(3)]),
    //     //     Wi_1_a: Wi_1.clone(),
    //     //     Wi_1_b: Wi_1.clone(),
    //     // }
    //     LayerRoundPoly::new(
    //         MultilinearExtension::fix_variables(&add_i, &[F::from(3)]),
    //         MultilinearExtension::fix_variables(&mul_i, &[F::from(3)]),
    //         Wi_1.clone(),
    //         Wi_1.clone(),
    //     )
    // }

    // #[test]
    // fn get_evals_by_mask_one_3_var() {
    //     let poly = test_3_var_poly();
    //     let nvars = 3;
    //     let cases = vec![(0b010, 3), (0b011, 4), (0b110, 7)];
    // 
    //     for (mask, expected_val) in cases {
    //         let res = get_evaluations_by_mask(&poly, mask, nvars);
    // 
    //         assert_eq!(res.len(), 1);
    //         assert_eq!(res[0], Fr::from(expected_val));
    //     }
    // }

    // #[test]
    // fn get_evals_by_mask_two_3_var() {
    //     let poly = test_3_var_poly();
    //     let nvars = 3;
    // 
    //     let cases = vec![
    //         (0b010, (3, 7)),
    //         (0b011, (4, 8)),
    //         (0b000, (1, 5)),
    //     ];
    // 
    //     for (mask, (coef_1, coef_2)) in cases {
    //         let res = get_evaluations_by_mask(&poly, mask, nvars - 1);
    // 
    //         assert_eq!(res.len(), 2);
    //         assert_eq!(res, [Fr::from(coef_1), Fr::from(coef_2)]);
    //     }
    // }
    // 
    // #[test]
    // fn to_two_degree() {
    //     let poly = test_3_var_poly();
    //     let nvars = 3;
    // 
    //     let cases = vec![
    //         (0b010, (3, 4)),
    //         (0b011, (4, 4)),
    //         (0b000, (1, 4)),
    //     ];
    // 
    //     for (mask, (coef_1, coef_2)) in cases {
    //         let res = to_two_or_one_degree(&poly, mask, nvars - 1);
    // 
    //         println!("res {:?}", res);
    //         assert_eq!(res.coeffs.len(), 2);
    //         assert_eq!(res.coeffs, [Fr::from(coef_1), Fr::from(coef_2)]);
    //     }
    // }

    // #[test]
    // fn get_partial_sum_poly_2_vars() {
    //     let round_poly = get_test_round_poly_2_vars::<Fr>();
    //     let partial_sum_poly_1 = round_poly.get_partial_sum_poly();
    // 
    //     assert_eq!(partial_sum_poly_1.coeffs, [Fr::from(67200), Fr::from(-32000), Fr::from(-35200)]);
    // 
    //     let round_poly = round_poly.fix_variables(&[Fr::from(32)]);
    //     let partial_sum_poly_2 = round_poly.get_partial_sum_poly();
    // 
    //     assert_eq!(partial_sum_poly_2.coeffs, [Fr::zero(), Fr::from(-24282300), Fr::from(-12719300)])
    // }
    // 
    // #[test]
    // fn get_partial_sum_poly_4_vars() {
    //     let round_poly = get_test_round_poly_4_vars::<Fr>();
    //     let partial_sum_poly_1 = round_poly.get_partial_sum_poly();
    // 
    //     assert_eq!(partial_sum_poly_1.coeffs, [Fr::from(-420), Fr::from(1330), Fr::from(50)]);
    // 
    //     let round_poly = round_poly.fix_variables(&[Fr::from(4)]);
    //     let partial_sum_poly_2 = round_poly.get_partial_sum_poly();
    // 
    //     assert_eq!(partial_sum_poly_2.coeffs, [Fr::from(5700), Fr::from(4200), Fr::from(-9900)]);
    // 
    //     let round_poly = round_poly.fix_variables(&[Fr::from(5)]);
    //     let partial_sum_poly_3 = round_poly.get_partial_sum_poly();
    // 
    //     assert_eq!(partial_sum_poly_3.coeffs, [Fr::from(-72000), Fr::from(-74400), Fr::from(-2400)]);
    // 
    //     let round_poly = round_poly.fix_variables(&[Fr::from(8)]);
    //     let partial_sum_poly_4 = round_poly.get_partial_sum_poly();
    // 
    //     assert_eq!(partial_sum_poly_4.coeffs, [Fr::from(0), Fr::from(-624240), Fr::from(-196560)]);
    // }

    // #[test]
    // fn generate_round_poly_params() {
    //     let a_vars = 2;
    //     let b_vars = 2;
    // 
    //     let ((a_vars_fixing, b_vars_fixing), res) = generate_round_poly_eval_parameters(a_vars, b_vars);
    // 
    //     assert_eq!(res, [
    //         (0b000, 0b000, 0b000),
    //         (0b001, 0b000, 0b001),
    //         (0b010, 0b000, 0b010),
    //         (0b011, 0b000, 0b011),
    //         (0b100, 0b001, 0b000),
    //         (0b101, 0b001, 0b001),
    //         (0b110, 0b001, 0b010),
    //         (0b111, 0b001, 0b011),
    //     ]);
    //     assert_eq!(a_vars_fixing, 1);
    //     assert_eq!(b_vars_fixing, 2);
    // }

    // #[test]
    // fn generate_round_poly_params_one_var() {
    //     let a_vars = 1;
    //     let b_vars = 2;
    // 
    //     let ((a_vars_fixing, b_vars_fixing), res) = generate_round_poly_eval_parameters(a_vars, b_vars);
    // 
    //     assert_eq!(res, [
    //         (0b000, 0b000, 0b000),
    //         (0b001, 0b000, 0b001),
    //         (0b010, 0b000, 0b010),
    //         (0b011, 0b000, 0b011),
    //     ]);
    //     assert_eq!(a_vars_fixing, 0);
    //     assert_eq!(b_vars_fixing, 2);
    // }
    // 
    // #[test]
    // fn generate_round_poly_params_zero_var() {
    //     let a_vars = 0;
    //     let b_vars = 2;
    // 
    //     let ((a_vars_fixing, b_vars_fixing), res) = generate_round_poly_eval_parameters(a_vars, b_vars);
    // 
    //     assert_eq!(res, [
    //         (0b000, 0b000, 0b000),
    //         (0b001, 0b000, 0b001),
    //     ]);
    //     assert_eq!(a_vars_fixing, 0);
    //     assert_eq!(b_vars_fixing, 1);
    // }

    // #[test]
    // fn round_poly_fix_polys() {
    //     let round_poly = get_test_round_poly_4_vars::<Fr>();
    //     let mask = 0b0001;
    //     let a_mask = 0b0000;
    //     let b_mask = 0b0001;
    // 
    //     let add_uni = to_two_or_one_degree(&round_poly.add_i, mask, 3);
    //     let mul_uni = to_two_or_one_degree(&round_poly.mul_i, mask, 3);
    //     let Wi_1_a_uni = to_two_or_one_degree(&round_poly.Wi_1_a, a_mask, 1);
    //     let Wi_1_b_uni = to_two_or_one_degree(&round_poly.Wi_1_b, b_mask, 2);
    // 
    //     assert_eq!(add_uni.coeffs, [Fr::from(-2), Fr::from(2)]);
    //     assert_eq!(mul_uni.coeffs, []);
    //     assert_eq!(Wi_1_a_uni.coeffs, [Fr::from(10), Fr::from(10)]);
    //     assert_eq!(Wi_1_b_uni.coeffs, [Fr::from(200)]);
    // }
    // 
    // #[test]
    // fn round_poly_fix_polys_other_mask() {
    //     let round_poly = get_test_round_poly_4_vars::<Fr>();
    //     let mask = 0b0101;
    //     let a_mask = 0b0001;
    //     let b_mask = 0b0001;
    // 
    //     let add_uni = to_two_or_one_degree(&round_poly.add_i, mask, 3);
    //     let mul_uni = to_two_or_one_degree(&round_poly.mul_i, mask, 3);
    //     let Wi_1_a_uni = to_two_or_one_degree(&round_poly.Wi_1_a, a_mask, 1);
    //     let Wi_1_b_uni = to_two_or_one_degree(&round_poly.Wi_1_b, b_mask, 2);
    // 
    //     assert_eq!(add_uni.coeffs, []);
    //     assert_eq!(mul_uni.coeffs, []);
    //     assert_eq!(Wi_1_a_uni.coeffs, [Fr::from(200), Fr::from(100)]);
    //     assert_eq!(Wi_1_b_uni.coeffs, [Fr::from(200)]);
    // }

    // #[test]
    // fn round_poly_get_required_params() {
    //     let round_poly = get_test_round_poly_4_vars::<Fr>();
    //     let ((a_eval_params, b_eval_params), evals) = generate_round_poly_eval_parameters(
    //         round_poly.Wi_1_a.num_vars,
    //         round_poly.Wi_1_b.num_vars,
    //     );
    // 
    //     assert_eq!(a_eval_params, 1);
    //     assert_eq!(b_eval_params, 2);
    //     assert_eq!(evals[1], (0b0001usize, 0b0000usize, 0b0001usize));
    // }
    // 
    // #[test]
    // fn line_test() {
    //     let b = &[Fr::from(3), Fr::from(7)];
    //     let c = &[Fr::from(4), Fr::from(5)];
    //     let l = line(b, c);
    // 
    //     assert_eq!(l.len(), 2);
    //     assert_eq!(l[0].coeffs, [Fr::from(3), Fr::one()]);
    //     assert_eq!(l[1].coeffs, [Fr::from(7), Fr::from(-2)]);
    // }

    // #[test]
    // fn restrict_to_line_test() {
    //     let b = &[Fr::from(3), Fr::from(7)];
    //     let c = &[Fr::from(4), Fr::from(5)];
    //     let l = line(b, c);
    //     let w = DenseMultilinearExtension::from_evaluations_vec(
    //         2,
    //         vec![10, 20, 200, 300].into_iter().map(Fr::from).collect(),
    //     );
    // 
    //     let w_restricted = restrict_poly(&l, w.clone());
    //     let r_star = Fr::from(12);
    // 
    //     let l_evals = l.iter().map(|li| li.evaluate(&r_star)).collect::<Vec<_>>();
    // 
    //     assert_eq!(l_evals.len(), 2);
    //     assert_eq!(l_evals[0], Fr::from(15));
    //     assert_eq!(l_evals[1], Fr::from(-17));
    // 
    //     let w_r_star = ark_poly::Polynomial::evaluate(&w, &l_evals);
    //     assert_eq!(w_r_star, Fr::from(-26020));
    // 
    //     let w_restricted_r_star = w_restricted.evaluate(&r_star);
    //     assert_eq!(w_restricted_r_star, w_r_star);
    // }

    // #[test]
    // fn add_i() {
    //     let input_a = InputGate {
    //         index: 0,
    //         value: Fr::from(10)
    //     };
    //     let input_b = InputGate {
    //         index: 1,
    //         value: Fr::from(200)
    //     };
    //     let input_c = InputGate {
    //         index: 2,
    //         value: Fr::from(20),
    //     };
    //     let input_d = InputGate {
    //         index: 3,
    //         value: Fr::from(300),
    //     };
    // 
    //     let add_gate_1 = Gate {
    //         inputs: &GateInput::Inputs((&input_a, &input_b)),
    //         executor: ExecutorGateEnum::Add(AddGate {})
    //     };
    //     let add_gate_2 = Gate {
    //         inputs: &GateInput::Inputs((&input_c, &input_d)),
    //         executor: ExecutorGateEnum::Add(AddGate {})
    //     };
    //     let mul_gate = Gate {
    //         inputs: &GateInput::Gates((&add_gate_1, &add_gate_2)),
    //         executor: ExecutorGateEnum::Mul(MulGate {})
    //     };
    // 
    //     let layer_1 = Layer {
    //         gates: vec![&add_gate_1, &add_gate_2],
    //     };
    //     let layer_2 = Layer {
    //         gates: vec![&mul_gate],
    //     };
    // 
    //     let circuit = Circuit {
    //         layers: vec![layer_1, layer_2],
    //     };
    // 
    //     let add_i = circuit.add_i(0);
    //     println!("add_i evals {:?}", add_i.evaluations);
    //     println!("add_i varsnum {}", add_i.num_vars);
    //     println!("add_i(0,0,0,0,1) {}", ark_poly::Polynomial::evaluate(&add_i, &vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]));
    //     println!("add_i(1,1,0,1,1) {}", ark_poly::Polynomial::evaluate(&add_i, &vec![Fr::one(), Fr::one(), Fr::zero(), Fr::one(), Fr::one()]));
    //     let add_i = circuit.add_i(1);
    //     println!("add_i evals {:?}", add_i.evaluations);
    //     println!("add_i varsnum {}", add_i.num_vars);
    // 
    //     let mul_i = circuit.mul_i(0);
    //     println!("mul_i evals {:?}", mul_i.evaluations);
    //     println!("mul_i varsnum {}", mul_i.num_vars);
    //     let mul_i = circuit.mul_i(1);
    //     println!("mul_i evals {:?}", mul_i.evaluations);
    //     println!("mul_i varsnum {}", mul_i.num_vars);
    // }
}