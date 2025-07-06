use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::ops::{Add, AddAssign, Mul, MulAssign, Neg};
use ark_ff::{Field, One, Zero};
use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial, EvaluationDomain, GeneralEvaluationDomain, MultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_std::iterable::Iterable;
use ark_test_curves::bls12_381::Fr;
use crate::sumcheck;
use crate::sumcheck::{interpolate_univariate_on_evals, sum_over_last_variable, SumCheckPoly};
// use ark_poly::polynomial::univariate::lagrange_interpolate;

trait ExecuteGate<F: Field> {
    fn execute(&self, a: &F, b: &F) -> F;
}

#[derive(Debug, Hash, Eq, PartialEq,)]
enum ExecutorGateEnum {
    Add(AddGate),
    Mul(MulGate),
}

#[derive(Debug, Hash, Eq, PartialEq,)]
struct MulGate {}

#[derive(Debug, Hash, Eq, PartialEq,)]
struct AddGate {}

impl<F: Field> ExecuteGate<F> for MulGate {
    fn execute(&self, a: &F, b: &F) -> F {
        *a * *b
    }
}

impl<F: Field> ExecuteGate<F> for AddGate {
    fn execute(&self, a: &F, b: &F) -> F {
        *a + *b
    }
}

#[derive(Debug, Hash, Eq, PartialEq)]
enum GateInput<'a, F: Field> {
    Gates((&'a Gate<'a, F>, &'a Gate<'a, F>)),
    Inputs((&'a InputGate<F>, &'a InputGate<F>))
}

#[derive(Debug, Hash, Eq, PartialEq)]
struct Gate<'a, F: Field> {
    inputs: &'a GateInput<'a, F>,
    executor: ExecutorGateEnum,
}
#[derive(Debug, Hash, Eq, PartialEq)]
struct InputGate<F: Field> {
    index: usize,
    value: F
}
#[derive(Debug)]
struct Layer<'a, F: Field> {
    gates: Vec<&'a Gate<'a, F>>,
}
#[derive(Debug)]
struct Circuit<'a, F: Field> {
    layers: Vec<Layer<'a, F>>,
}

#[derive(Debug)]
struct Solution<F: Field> {
    evaluations: Vec<Vec<F>>,
    inputs: Vec<F>,
}

#[derive(Debug, Clone)]
struct RoundPoly<F: Field> {
    add_i: DenseMultilinearExtension<F>,
    mul_i: DenseMultilinearExtension<F>,
    Wi_1_a: DenseMultilinearExtension<F>,
    Wi_1_b: DenseMultilinearExtension<F>,
}

fn get_bits(k: usize, bit_len: usize) -> Vec<usize> {
    let mut result = vec![];
    for i in (0..bit_len).rev() {
        let kk = k >> i;
        let res = kk & 1;
        result.push(res);
    }

    result
}

pub fn get_reversed_vars_poly<F: Field>(poly: &DenseMultilinearExtension<F>) -> DenseMultilinearExtension<F> {
    let n = MultilinearExtension::num_vars(poly);
    let mut res_poly = poly.clone();


    if n == 0 || n == 1 {
        return res_poly;
    }

    for i in 0..n / 2 {
        res_poly.relabel_in_place(i, n - i - 1,1);
    }

    // res_poly.relabel_in_place(0, n / 2, n / 2);

    res_poly
}

fn to_f<F: Field>(arr: Vec<u64>) -> Vec<F> {
    arr.into_iter().map(F::from).collect()
}

fn interpolate_or_const<F: Field>(evals: &[F]) -> DensePolynomial<F> {
    if evals.len() == 1 {
        DensePolynomial::from_coefficients_slice(evals)
    } else if evals.len() == 2 {
        interpolate_univariate_on_evals(evals.try_into().unwrap())
    } else {
        panic!("interpolate_or_const")
    }
}

fn get_evaluations_by_mask<F: Field>(
    poly: &DenseMultilinearExtension<F>,
    mask: usize,
    nvars: usize,
) -> Vec<F> {
    if nvars == MultilinearExtension::num_vars(poly) {
        let bits = get_bits(mask, nvars).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
        return vec![Polynomial::evaluate(poly, &bits)];
    }

    let bit_set = 1 << nvars;
    // println!("nvars {} bit_set {:#08b}", nvars, bit_set);
    let mut mle_evals = poly.to_evaluations();
    remap_to_reverse_bits_indexing(&mut mle_evals, MultilinearExtension::num_vars(poly));
    // println!("mle remapped {:?}", mle_evals);

    // let index_1 = reverse_bits(bit_set | mask, MultilinearExtension::num_vars(poly));
    // let index_2 = reverse_bits(mask, MultilinearExtension::num_vars(poly));
    let index_1 = bit_set | mask;
    let index_2 = mask;

    // println!("index_1 {:#08b} {} index_2 {:#08b} {}", index_1, index_1,  index_2, index_2);

    let mut new_evals = vec![];

    new_evals.push(mle_evals[index_2]);
    new_evals.push(mle_evals[index_1]);

    new_evals
}

fn to_two_or_one_degree<F: Field>(
    poly: &DenseMultilinearExtension<F>,
    mask: usize,
    nvars: usize,
) -> DensePolynomial<F> {
    interpolate_or_const(&get_evaluations_by_mask(poly, mask, nvars))
}

fn generate_round_poly_eval_parameters(
    a_nvars: usize,
    b_nvars: usize,
) -> ((usize, usize), Vec<(usize, usize, usize)>){
    let (a_vars_fixing, b_vars_fixing) = {
        if a_nvars > 0 {
            (a_nvars - 1, b_nvars)
        } else {
            (a_nvars, b_nvars - 1)
        }
    };

    let combined_domain = 1 << (a_vars_fixing + b_vars_fixing); // for example 1111 where [1] = a and [111] = b
    let b_mask = (1 << b_vars_fixing) - 1; // 0111
    let a_mask = (1 << a_vars_fixing + b_vars_fixing) - 1 - b_mask; // 1000

    // println!("a_mask {:#08b} b_mask {:#08b}", a_mask, b_mask);

    let mut res = vec![];

    for i in 0..combined_domain {
        let combined_param = i;
        let b_param = i & b_mask;
        let a_param = (i & a_mask) >> b_vars_fixing;

        // println!("i {:#08b} a_param {:#08b} b_param {:#08b}", i, a_param, b_param);

        res.push((combined_param, a_param, b_param))
    }

    ((a_vars_fixing, b_vars_fixing), res)
}

impl<F: Field> Layer<'_, F> {
    fn get_gates_n(&self) -> usize {
        self.gates.len()
    }
}

fn test_intp<F: Field>(evals: &[F]) -> DensePolynomial<F> {
    // let mut current_sc_poly = self.clone();
    // let num_vars = current_sc_poly.num_vars();
    //
    // let evals = self.get_evaluations();
    // println!("evals {:?}", evals);

    let (zeroes, ones) = evals.iter().enumerate().fold((F::zero(), F::zero()), |(zeroes, ones), (i, eval)| {
        if i % 2 == 0 {
            (zeroes + eval, ones)
        } else {
            (zeroes, ones + eval)
        }
    });

    interpolate_univariate_on_evals(&[zeroes, ones])
}

impl<F: Field> SumCheckPoly<F> for RoundPoly<F> {
    fn get_evaluations(&self) -> Vec<F> {
        let mut res = vec![];
        // let num_vars = MultilinearExtension::num_vars(&self.add_i);
        let a_num_vars = MultilinearExtension::num_vars(&self.Wi_1_a);
        let b_num_vars = MultilinearExtension::num_vars(&self.Wi_1_b);
        let ab_num_vars = a_num_vars + b_num_vars;
        let a_domain = 1 << a_num_vars;
        let b_domain = 1 << b_num_vars;
        // println!("a_numvars {} b_num_vars {}", a_num_vars, b_num_vars);

        for a in 0..a_domain {
            for b in 0..b_domain {
                // let a_bits = get_bits(a, n / 2).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
                // let b_bits = get_bits(b, n / 2).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();;
                // let ab_combined_mask = (a << (n / 2)) | b;
                // let ab_bits = get_bits(ab_combined_mask, n).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
                //
                // let add_eval = Polynomial::evaluate(&self.add_i, &ab_bits)
                //     * (Polynomial::evaluate(&self.Wi_1, &a_bits) + Polynomial::evaluate(&self.Wi_1, &b_bits));
                // let mul_eval = Polynomial::evaluate(&self.mul_i, &ab_bits)
                //     * (Polynomial::evaluate(&self.Wi_1, &a_bits) * Polynomial::evaluate(&self.Wi_1, &b_bits));
                //
                // res.push(add_eval + mul_eval);
                let ab_combined_mask = (a << b_num_vars) | b;
                let ab_combined_mask_reversed = reverse_bits(ab_combined_mask, ab_num_vars);
                // println!("a {:#010b} b {:#010b} ab_combined_mask {:#010b}", a, b, ab_combined_mask_reversed);
                let ab_bits = get_bits(ab_combined_mask_reversed, ab_num_vars).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
                // let a_bits = get_bits(a, a_num_vars);
                // let b_bits = get_bits(b, b_num_vars);
                // let ab_bits = a_bits.into_iter()
                //     .chain(b_bits.into_iter())
                //     .map(|e| F::from(e as u64))
                //     .collect::<Vec<_>>();
                let eval = RoundPoly::evaluate(&self, &ab_bits);

                res.push(eval);
            }
        }

        res
    }

    fn num_vars(&self) -> usize {
        MultilinearExtension::num_vars(&self.add_i)
    }

    fn get_partial_sum_poly(&self) -> DensePolynomial<F> {
        let a_vars = MultilinearExtension::num_vars(&self.Wi_1_a);
        let b_vars = MultilinearExtension::num_vars(&self.Wi_1_b);
        // let (a_vars_fixing, b_vars_fixing) = {
        //     if a_vars > 0 {
        //         (a_vars - 1, b_vars)
        //     } else {
        //         (a_vars, b_vars - 1)
        //     }
        // };
        // let a_domain = 1 << a_vars_fixing;
        // let b_domain = 1 << b_vars_fixing;
        let Wi_1a_rev = get_reversed_vars_poly(&self.Wi_1_a);
        let Wi_1b_rev = get_reversed_vars_poly(&self.Wi_1_b);
        let add_rev = get_reversed_vars_poly(&self.add_i);
        let mul_rev = get_reversed_vars_poly(&self.mul_i);

        let mut res_poly = DensePolynomial::from_coefficients_vec(vec![F::zero()]);

        let ((a_vars_fixing, b_vars_fixing), eval_params) = generate_round_poly_eval_parameters(a_vars, b_vars);

        for (common_params, a_params, b_params) in eval_params {
            println!("common_params {:#08b}", common_params);
            let Wi_a_uni = to_two_or_one_degree(&self.Wi_1_a, a_params, a_vars_fixing);
            println!("a params {:#08b} Wi_a_uni {:?}", a_params, Wi_a_uni);
            let Wi_b_uni = to_two_or_one_degree(&self.Wi_1_b, b_params, b_vars_fixing);
            println!("b_params {:#08b} Wi_b_uni {:?}", b_params, Wi_b_uni);
            let add_i_uni = to_two_or_one_degree(&self.add_i, common_params, a_vars_fixing + b_vars_fixing);
            println!("add_i_uni {:?}", add_i_uni);
            let mul_uni = to_two_or_one_degree(&self.mul_i, common_params, a_vars_fixing + b_vars_fixing);
            println!("mul_i_uni {:?}", mul_uni);

            let add_poly = add_i_uni.naive_mul(&(&Wi_a_uni + &Wi_b_uni));
            let mul_poly = mul_uni.naive_mul(&(Wi_a_uni.naive_mul(&Wi_b_uni)));
            let local_res_poly = add_poly + mul_poly;

            println!("local_res_poly {:?}\n", local_res_poly);

            res_poly = res_poly + local_res_poly;
        }
        // println!("a_vars_fixing {}\nb_vars_fixing {}", a_vars_fixing, b_vars_fixing);
        // println!("Wi_a {:?}", self.Wi_1_a);
        // println!("Wi_b {:?}", self.Wi_1_b);
        // println!("add_i {:?}", self.add_i);
        // println!("mul_i {:?}", self.mul_i);
        //
        // let combined_domain = 1 << (a_vars_fixing + b_vars_fixing);
        //
        // let ab_domain = a_domain + b_domain;
        // println!("combined_domain {}", combined_domain);
        // let a_mask = (1 << a_vars_fixing) - 1;
        // let b_mask = (1 << b_vars_fixing) - 1;
        // println!("a_mask {:#08b} {:#08b}", a_mask, b_mask);
        //
        // for i in 0..combined_domain {
        //     let ab_mask = i;
        //     println!("\nab_mask {:#08b}", ab_mask);
        //     let ai = i & a_mask;
        //     let bi = i & b_mask;
        //
        //     println!("ai {:#08b}", ai);
        //     println!("bi {:#08b}", bi);
        //     let Wi_a_dense = to_two_or_one_degree(&self.Wi_1_a, ai, a_vars_fixing);
        //     println!("Wi_a_dense {:?}", Wi_a_dense);
        //     let Wi_b_dense = to_two_or_one_degree(&self.Wi_1_b, bi, b_vars_fixing);
        //     println!("Wi_b_dense {:?}", Wi_b_dense);
        //     let add_dense = to_two_or_one_degree(&self.add_i, ab_mask, a_vars_fixing + b_vars_fixing);
        //     println!("ab_mask {:#08b} add {:?}\nadd_dense {:?}", ab_mask, self.add_i.evaluations, add_dense);
        //     let mul_dense = to_two_or_one_degree(&self.mul_i, ab_mask, a_vars_fixing + b_vars_fixing);
        //     println!("mul_dense {:?}", mul_dense);
        //
        //     let add_poly = add_dense.naive_mul(&(&Wi_a_dense + &Wi_b_dense));
        //     let mul_poly = mul_dense.naive_mul(&(Wi_a_dense.naive_mul(&Wi_b_dense)));
        //     let local_res_poly = add_poly + mul_poly;
        //
        //     println!("\nlocal_res_poly {:?}", local_res_poly);
        //
        //     res_poly = res_poly + local_res_poly;
        // }
        //
        return res_poly;

        // println!("res poly new inter {:?}", res_poly);
        let res_poly_0 = res_poly.evaluate(&F::zero());
        let res_poly_1 = res_poly.evaluate(&F::one());
        // println!("res poly(0) {}", res_poly_0);
        // println!("res poly(1) {}", res_poly_1);
        // println!("res poly sum {}", res_poly_0 + res_poly_1);

        res_poly
    }

    // fn get_partial_sum_poly(&self) -> DensePolynomial<F> {
    //     // Choose 3 distinct field elements manually
    //     let t_values = vec![F::from(0u64), F::from(1u64), F::from(2u64)];
    //     let mut evals = vec![];
    //
    //     for t in &t_values {
    //         let fixed = self.fix_variables(&[*t]);
    //         let sum = fixed.get_evaluations().iter().copied().sum::<F>();
    //         evals.push(sum);
    //     }
    //
    //     // Interpolate over (t_values, evals)
    //     // let dense_poly: DensePolynomial<F> = lagrange_interpolate(&t_values, &evals);
    //     // dense_poly.into() // convert to sparse
    //     let mut res = vec![];
    //     let lagrange_polys = get_lagrange_polys(evals.len());
    //     for lp in lagrange_polys {
    //         res.push(lp.int)
    //     }
    // }

    fn fix_variables(&self, partial_point: &[F]) -> Self {
        println!("fixing {:?}", partial_point);
        if MultilinearExtension::num_vars(&self.Wi_1_a) > 0 {
            // println!("before fix mul {:?}", &self.mul_i);
            // println!("mul (0,0) {}", &Polynomial::evaluate(&self.mul_i, &vec![F::zero(), F::zero()]));
            // println!("mul (0,1) {}", &Polynomial::evaluate(&self.mul_i, &vec![F::zero(), F::one()]));
            // println!("mul (1,0) {}", &Polynomial::evaluate(&self.mul_i, &vec![F::one(), F::zero()]));
            // println!("mul (1,1) {}", &Polynomial::evaluate(&self.mul_i, &vec![F::one(), F::one()]));
            let res = RoundPoly {
                add_i: MultilinearExtension::fix_variables(&self.add_i, partial_point),
                mul_i: MultilinearExtension::fix_variables(&self.mul_i, partial_point),
                Wi_1_a: MultilinearExtension::fix_variables(&self.Wi_1_a, partial_point),
                Wi_1_b: self.Wi_1_b.clone(),
            };

            // println!("res.add {:?}", res.add_i);
            // println!("res.mul {:?}", res.mul_i);
            // println!("res.Wi_1a {:?}", res.Wi_1_a);
            // println!("res.Wi_1b {:?}", res.Wi_1_b);

            res
        } else {
            RoundPoly {
                add_i: MultilinearExtension::fix_variables(&self.add_i, partial_point),
                mul_i: MultilinearExtension::fix_variables(&self.mul_i, partial_point),
                Wi_1_a: self.Wi_1_a.clone(),
                Wi_1_b: MultilinearExtension::fix_variables(&self.Wi_1_b, partial_point)
            }
        }
    }

    fn evaluate(&self, point: &[F]) -> F {
        let binding = point.to_vec();
        let a_vars = MultilinearExtension::num_vars(&self.Wi_1_a);
        let (a_bits, b_bits) = binding.split_at(a_vars);
        let a_bits = a_bits.to_vec();
        let b_bits = b_bits.to_vec();
        // println!("evaluate a_bits {:?} b_bits {:?}", a_bits, b_bits);
        // let a_bits = get_bits(a, n / 2).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
        // let b_bits = get_bits(b, n / 2).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();;
        // let ab_combined_mask = (a << (n / 2)) | b;
        // let ab_bits = get_bits(ab_combined_mask, n).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();

        let add_eval = Polynomial::evaluate(&self.add_i, &point.to_vec())
            * (Polynomial::evaluate(&self.Wi_1_a, &a_bits) + Polynomial::evaluate(&self.Wi_1_b, &b_bits));
        let mul_eval = Polynomial::evaluate(&self.mul_i, &point.to_vec())
            * (Polynomial::evaluate(&self.Wi_1_a, &a_bits) * Polynomial::evaluate(&self.Wi_1_b, &b_bits));

        add_eval + mul_eval
    }
}

// fn sum_over_last_sc<F: Field>(poly: RoundPoly<F>) {
//     if
// }

impl<'a, F: Field> Circuit<'a, F> {
    fn solve(&self) -> Solution<F> {
        let mut solvedMap = HashMap::new();
        let mut evaluations = vec![];
        let mut inputs = vec![];
        let mut inputs_set = std::collections::HashSet::<&InputGate<F>>::new();
        
        for layer in &self.layers {
            let mut layer_evaluations = vec![];
            for gate in &layer.gates {
                match *gate.inputs {
                    GateInput::Gates((a, b)) => {
                        let a_res = solvedMap.remove(&a).unwrap();
                        let b_res = solvedMap.remove(&b).unwrap();
                        
                        let res = execute(&gate.executor, &a_res, &b_res);
                        layer_evaluations.push(res);
                        solvedMap.insert(*gate, res);
                    },
                    GateInput::Inputs((a, b)) => {
                        if !inputs_set.contains(a) {
                            inputs.push(a.value);
                            inputs_set.insert(a);
                        }
                        
                        if !inputs_set.contains(b) {
                            inputs.push(b.value);
                            inputs_set.insert(b);
                        }
                        
                        let res = execute(&gate.executor, &a.value, &b.value);
                        layer_evaluations.push(res);
                        solvedMap.insert(*gate, res);
                    }
                }
            }
            evaluations.push(layer_evaluations);
        }
        
        Solution {
            evaluations,
            inputs,
        }
    }
    
    fn in_1(&self, gate_i: usize, layer_i: usize) -> usize {
        self.bottom_left_count(gate_i, layer_i)
    }

    fn in_2(&self, gate_i: usize, layer_i: usize) -> usize {
        self.bottom_left_count(gate_i, layer_i) + 1
    }

    fn bottom_left_count(&self, gate_i: usize, layer_i: usize) -> usize {
        let layer = &self.layers[layer_i];
        // let gate = &layer.gates[gate_i];

        gate_i * 2
    }

    fn gates_n_at_layer_i_1(&self, layer_i: usize) -> usize {
        let layer = &self.layers[layer_i];

        match layer_i {
            0 => layer.gates.iter().fold(HashSet::new(),  |mut res, gate| {
                if let GateInput::Inputs((a, b)) = gate.inputs {
                    res.insert(a);
                    res.insert(b);
                    res
                } else {
                    panic!();
                }
            }).iter().count(),
            _ => layer.gates.len() * 2,
        }
    }

    fn get_gate_at_layer_i_1(&self, layer_i: usize, gate_i: usize) -> Option<&Gate<'a, F>> {
        if layer_i == 0 {
            panic!();
        }

        let layer = &self.layers[layer_i - 1];


        layer.gates.get(gate_i).map(|gate| *gate)
    }

    fn mul_i(&self, layer_i: usize) -> DenseMultilinearExtension<F> {
        let layer = &self.layers[layer_i];
        let mut gates_n = layer.gates.len();
        if !gates_n.is_power_of_two() {
            gates_n = gates_n.next_power_of_two();
        }
        // let
        let vars_num = f64::from(gates_n as u32).log2().ceil() as usize;
        // println!("gates_num {}", gates_n);
        // println!("vars {}", vars_num);
        // println!("trailing zeroes {}", (gates_n as u64).trailing_zeros());
        let domain_size = 1 << vars_num;
        // println!("domain_size {}", domain_size);
        // let bottom_layer = &self.layers[layer_i - 1];
        // let mut bottom_gates_n = bottom_layer.gates.len();
        let mut bottom_gates_n = self.gates_n_at_layer_i_1(layer_i);
        if !bottom_gates_n.is_power_of_two() {
            bottom_gates_n = bottom_gates_n.next_power_of_two();
        }
        let bottom_vars_num = f64::from(bottom_gates_n as u32).log2().ceil() as usize;
        let bottom_domain_size = 1 << bottom_vars_num;
        // println!("bottom_gates_n {} bottom_vars_num {} bottom_domain_size {}", bottom_gates_n, bottom_vars_num, bottom_domain_size);

        let mut evals = vec![];

        for z in 0..domain_size {
            for a in 0..bottom_domain_size {
                for b in 0..bottom_domain_size {
                    if let Some(current_gate) = layer.gates.get(z) {
                        let inputs = current_gate.inputs;
                        match inputs {
                            GateInput::Inputs((a_input, b_input)) => {
                                if a_input.index == a && b_input.index == b && matches!(current_gate.executor, ExecutorGateEnum::Mul(_)) {
                                    evals.push(F::one());
                                } else {
                                    evals.push(F::zero());
                                }
                            },
                            GateInput::Gates((a_gate, b_gate)) => {
                                if let (Some(bottom_a), Some(bottom_b)) = (self.get_gate_at_layer_i_1(layer_i, a), self.get_gate_at_layer_i_1(layer_i, b)) {
                                    if bottom_a == *a_gate && bottom_b == *b_gate && matches!(current_gate.executor, ExecutorGateEnum::Mul(_)) {
                                        evals.push(F::one());
                                    } else {
                                        evals.push(F::zero());
                                    }
                                } else {
                                    evals.push(F::zero());
                                }
                            }
                        };
                    } else {
                        evals.push(F::zero());
                    }
                }
            }
        }

        remap_to_reverse_bits_indexing(&mut evals, vars_num + 2 * bottom_vars_num);

        println!("mul evals {:?}", evals);

        DenseMultilinearExtension::from_evaluations_vec(vars_num + 2 * bottom_vars_num, evals)
    }

    fn add_i(&self, layer_i: usize) -> DenseMultilinearExtension<F> {
        let layer = &self.layers[layer_i];
        let mut gates_n = layer.gates.len();
        if !gates_n.is_power_of_two() {
            gates_n = gates_n.next_power_of_two();
        }
        // let
        let vars_num = f64::from(gates_n as u32).log2().ceil() as usize;
        // println!("gates_num {}", gates_n);
        // println!("vars {}", vars_num);
        // println!("trailing zeroes {}", (gates_n as u64).trailing_zeros());
        let domain_size = 1 << vars_num;
        // println!("domain_size {}", domain_size);
        // let bottom_layer = &self.layers[layer_i - 1];
        // let mut bottom_gates_n = bottom_layer.gates.len();
        let mut bottom_gates_n = self.gates_n_at_layer_i_1(layer_i);
        if !bottom_gates_n.is_power_of_two() {
            bottom_gates_n = bottom_gates_n.next_power_of_two();
        }
        let bottom_vars_num = f64::from(bottom_gates_n as u32).log2().ceil() as usize;
        let bottom_domain_size = 1 << bottom_vars_num;
        // println!("bottom_gates_n {} bottom_vars_num {} bottom_domain_size {}", bottom_gates_n, bottom_vars_num, bottom_domain_size);

        let mut evals = vec![];

        for z in 0..domain_size {
            for a in 0..bottom_domain_size {
                for b in 0..bottom_domain_size {
                    if let Some(current_gate) = layer.gates.get(z) {
                        let inputs = current_gate.inputs;
                        match inputs {
                            GateInput::Inputs((a_input, b_input)) => {
                                if a_input.index == a && b_input.index == b && matches!(current_gate.executor, ExecutorGateEnum::Add(_)) {
                                    evals.push(F::one());
                                } else {
                                    evals.push(F::zero());
                                }
                            },
                            GateInput::Gates((a_gate, b_gate)) => {
                                if let (Some(bottom_a), Some(bottom_b)) = (self.get_gate_at_layer_i_1(layer_i, a), self.get_gate_at_layer_i_1(layer_i, b)) {
                                    if bottom_a == *a_gate && bottom_b == *b_gate && matches!(current_gate.executor, ExecutorGateEnum::Add(_)) {
                                        evals.push(F::one());
                                    } else {
                                        evals.push(F::zero());
                                    }
                                } else {
                                    evals.push(F::zero());
                                }
                            }
                        };
                    } else {
                        evals.push(F::zero());
                    }
                }
            }
        }

        remap_to_reverse_bits_indexing(&mut evals, vars_num + 2 * bottom_vars_num);

        println!("evals {:?}", evals);

        DenseMultilinearExtension::from_evaluations_vec(vars_num + 2 * bottom_vars_num, evals)
    }

    // fn mul_i(&self, layer_i: usize) -> DenseMultilinearExtension<F> {
    //     let layer = &self.layers[layer_i];
    //     let mut gates_n = layer.gates.len();
    //     if !gates_n.is_power_of_two() {
    //         gates_n = gates_n.next_power_of_two();
    //     }
    //     // let
    //     let vars_num = f64::from(gates_n as u32).log2().ceil() as usize;
    //     let domain_size = 1 << vars_num;
    //     let bottom_layer = &self.layers[layer_i - 1];
    //     let mut bottom_gates_n = bottom_layer.gates.len();
    //     if !bottom_gates_n.is_power_of_two() {
    //         bottom_gates_n = bottom_gates_n.next_power_of_two();
    //     }
    //     let bottom_vars_num = f64::from(bottom_gates_n as u32).log2().ceil() as usize;
    //     let bottom_domain_size = 1 << bottom_vars_num;
    //
    //     let mut evals = vec![];
    //
    //     for z in 0..domain_size {
    //         for a in 0..bottom_domain_size {
    //             for b in 0..bottom_domain_size {
    //                 if let Some(current_gate) = layer.gates.get(z) {
    //                     let inputs = current_gate.inputs;
    //                     match inputs {
    //                         GateInput::Inputs((a, b)) => {
    //                             panic!("not now")
    //                         },
    //                         GateInput::Gates((a_gate, b_gate)) => {
    //                             if let (Some(bottom_a), Some(bottom_b)) = (bottom_layer.gates.get(a), bottom_layer.gates.get(b)) {
    //                                 if *bottom_a == *a_gate && *bottom_b == *b_gate && matches!(current_gate.executor, ExecutorGateEnum::Mul(_)) {
    //                                     evals.push(F::one());
    //                                 } else {
    //                                     evals.push(F::zero());
    //                                 }
    //                             } else {
    //                                 evals.push(F::zero());
    //                             }
    //                         }
    //                     };
    //                 } else {
    //                     evals.push(F::zero());
    //                 }
    //             }
    //         }
    //     }
    //
    //     println!("mul evals {:?}", evals);
    //
    //     DenseMultilinearExtension::from_evaluations_vec(vars_num + 2 * bottom_vars_num, evals)
    // }

    fn add(&self, layer_i: usize, gate_i: usize, inputs: GateInput<F>) -> bool {
        if let Some(gate) = self.get_gate_for_inputs(layer_i, gate_i, inputs) {
            return match gate.executor {
                ExecutorGateEnum::Add(_) => true,
                _ => false,
            }
        }

        false
    }

    fn mul(&self, layer_i: usize, gate_i: usize, inputs: GateInput<F>) -> bool {
        if let Some(gate) = self.get_gate_for_inputs(layer_i, gate_i, inputs) {
            return match gate.executor {
                ExecutorGateEnum::Mul(_) => true,
                _ => false,
            }
        }

        false
    }

    fn get_gate_for_inputs(&self, layer_i: usize, gate_i: usize, inputs: GateInput<F>) -> Option<&Gate<'a, F>> {
        let layer = &self.layers[layer_i];
        let gate = layer.gates[gate_i];

        let left_in = self.in_1(gate_i, layer_i);
        let right_in = self.in_2(gate_i, layer_i);

        // if left == left_in && right == right_in {
        //     return Some(gate);
        // }
        if *gate.inputs == inputs {
            return Some(gate);
        }

        None
    }
}

fn execute<F: Field>(executor: &ExecutorGateEnum, a: &F, b: &F) -> F {
    match executor {
        ExecutorGateEnum::Add(add) => add.execute(a, b),
        ExecutorGateEnum::Mul(mul) => mul.execute(a, b),
    }
}

pub fn test_gkr() {
    // (a + b) * (c * c)
    let input_a = InputGate {
        index: 0,
        value: Fr::from(10)
    };
    let input_b = InputGate {
        index: 1,
        value: Fr::from(200)
    };
    let input_c = InputGate {
        index: 2,
        value: Fr::from(20),
    };
    let input_d = InputGate {
        index: 3,
        value: Fr::from(300),
    };

    let add_gate_1 = Gate {
        inputs: &GateInput::Inputs((&input_a, &input_b)),
        executor: ExecutorGateEnum::Add(AddGate {})
    };
    let add_gate_2 = Gate {
        inputs: &GateInput::Inputs((&input_c, &input_d)),
        executor: ExecutorGateEnum::Add(AddGate {})
    };
    let mul_gate = Gate {
        inputs: &GateInput::Gates((&add_gate_1, &add_gate_2)),
        executor: ExecutorGateEnum::Mul(MulGate {})
    };

    println!("-12719300 {}\n-24282300 {}", Fr::from(-12719300), Fr::from(-24282300));

    let layer_1 = Layer {
        gates: vec![&add_gate_1, &add_gate_2],
    };
    let layer_2 = Layer {
        gates: vec![&mul_gate],
    };

    let circuit = Circuit {
        layers: vec![layer_1, layer_2],
    };

    let solution = circuit.solve();

    println!("{:?}", solution);

    let mut test_poly = DenseMultilinearExtension::from_evaluations_vec(
        3,
        vec![1,2,3,4,5,6,7,8].into_iter().map(Fr::from).collect(),
    );
    
    let random_points = vec![3,22,12,93,8181,12398,123]
        .into_iter()
        .map(Fr::from)
        .collect::<Vec<_>>();

    prove(circuit, solution, &random_points);
}

fn prove<F: Field>(
    circuit: Circuit<F>,
    solution: Solution<F>,
    random_points: &[F],
) {
    let mut spent_points = 0;
    for i in (0..solution.evaluations.len()).rev() {
        // println!("solution.evaluations[i] {:?}", solution.evaluations[i]);
        let mut padded = pad_with_zeroes(&solution.evaluations[i]);
        let vars_num = padded.len().ilog2();
        remap_to_reverse_bits_indexing(&mut padded, vars_num as usize);
        // println!("vars_num {}", vars_num);

        let circuit_layer = &circuit.layers[i];
        let Wi = interpolate(&padded);

        let r0 = &random_points[spent_points..MultilinearExtension::num_vars(&Wi) + spent_points];
        spent_points += r0.len();
        let m0 = Polynomial::evaluate(&Wi, &r0.to_vec());

        let add_poly = circuit.add_i(i);
        let mul_poly = circuit.mul_i(i);
        let Wi_1 = {
            if i == 0 {
                // println!("solution.inputs {:?}", solution.inputs);
                // println!("interpolate(&solution.inputs) {:?}", interpolate(&pad_with_zeroes(&solution.inputs)));
                let mut wi_1_evals = pad_with_zeroes(&solution.inputs);
                let len = wi_1_evals.len();
                remap_to_reverse_bits_indexing(&mut wi_1_evals, len.ilog2() as usize);
                let res = interpolate(&wi_1_evals);
                res
            } else {
                interpolate(&solution.evaluations[i - 1])
            }
        };

        println!("Wi_1 {:?}", Wi_1);

        if MultilinearExtension::num_vars(&Wi_1) == 2 {
            println!("Wi_1 (0,0) {:?}", Polynomial::evaluate(&Wi_1, &vec![F::zero(), F::zero()]));
            println!("Wi_1 (0,1) {:?}", Polynomial::evaluate(&Wi_1, &vec![F::zero(), F::one()]));
            println!("Wi_1 (1,0) {:?}", Polynomial::evaluate(&Wi_1, &vec![F::one(), F::zero()]));
            println!("Wi_1 (1,1) {:?}", Polynomial::evaluate(&Wi_1, &vec![F::one(), F::one()]));
        }

        if add_poly.num_vars == 5 {
            println!("add (0,0,0,0,1) {}", Polynomial::evaluate(&add_poly, &vec![F::zero(), F::zero(), F::zero(), F::zero(), F::one()]));
            println!("add (1,1,0,1,1) {}", Polynomial::evaluate(&add_poly, &vec![F::one(), F::one(), F::zero(), F::one(), F::one()]));
            println!("add (1,1,1,1,1) {}", Polynomial::evaluate(&add_poly, &vec![F::one(), F::one(), F::one(), F::one(), F::one()]));
        }

        let add_fixed = MultilinearExtension::fix_variables(&add_poly, &r0);
        let mul_fixed = MultilinearExtension::fix_variables(&mul_poly, &r0);
        println!("r0 {:?}", r0);
        println!("add fixed {:?}", add_fixed.get_evaluations());

        if add_fixed.num_vars == 4 {
            println!("add fixed (0,0,0,1) {}", Polynomial::evaluate(&add_fixed, &vec![ F::zero(), F::zero(), F::zero(), F::one()]));
            println!("add fixed (1,0,1,1) {}", Polynomial::evaluate(&add_fixed, &vec![ F::one(), F::zero(), F::one(), F::one()]));
        }
        // println!("Wi_1 {:?}", Wi_1);

        let interpolated = interpolate_round_for_sc(add_fixed.clone(), mul_fixed.clone(), Wi_1.clone());
        let sc_poly = RoundPoly {
            add_i: add_fixed.clone(),
            mul_i: mul_fixed.clone(),
            Wi_1_a: Wi_1.clone(),
            Wi_1_b: Wi_1.clone(),
        };
        // println!("sc eval {}", sc_poly.evaluate(&[F::zero(), F::one()]));

        println!("m0 {:?}", m0);
        // println!("interpolated {:?}", interpolated.evaluations);
        //
        let test_params = vec![(0,0), (0,1), (1,0), (1,1), (123123123, 456456456)]
            .into_iter()
            .map(|(a,b)| (F::from(a), F::from(b)))
            .collect::<Vec<_>>();

        // println!("sc evals {:?}", sc_poly.get_evaluations());
        // println!("interpolated evals {:?}", interpolated.evaluations);

        for (a,b) in test_params {
            // println!("interpolated ({},{}) {}", a, b, Polynomial::evaluate(&interpolated, &vec![a,b]));
            // println!(
            //     "base ({},{}) {}",a,b,
            //     Polynomial::evaluate(&add_fixed, &vec![a, b]) * (Polynomial::evaluate(&Wi_1, &vec![a]) + Polynomial::evaluate(&Wi_1, &vec![b]))
            //         + Polynomial::evaluate(&mul_fixed, &vec![a, b]) * (Polynomial::evaluate(&Wi_1, &vec![a]) * Polynomial::evaluate(&Wi_1, &vec![b]))
            // );
            // println!("sc ({}, {}) {}", a, b, sc_poly.evaluate(&vec![a, b]));
        }

        // let partial = sc_poly.get_partial_sum_poly();
        // let partial_0 = partial.evaluate(&F::zero());
        // let partial_1 = partial.evaluate(&F::one());
        // println!("partial(0) {:?}", partial_0);
        // println!("partial(1) {:?}", partial_1);
        // println!("partial {}", partial_0 + partial_1);
        // println!("sc sum evals {}", sc_poly.get_evaluations().iter().sum::<F>());
        //
        // let sc_poly = sc_poly.fix_variables(&[F::from(999888)]);
        // let partial = sc_poly.get_partial_sum_poly();
        // let partial_0 = partial.evaluate(&F::zero());
        // let partial_1 = partial.evaluate(&F::one());
        // println!("partial(0) {:?}", partial_0);
        // println!("partial(1) {:?}", partial_1);
        // println!("partial {}", partial_0 + partial_1);
        // println!("sc sum evals {}", sc_poly.get_evaluations().iter().sum::<F>());



        let (used_r, used_polys) = sumcheck::prove_2(sc_poly.clone());
        // let (used_r_inter, used_poly_inter) = sumcheck::prove(interpolated.clone());
        // println!("used_polys {:?}", used_polys);

        let (b, c) = used_r.split_at(used_r.len() / 2);
        let l = line(b, c);

        // println!("l {:?}", l[0]);
        // println!("l(0) {}", l[0].evaluate(&F::zero()));
        // println!("l(1) {}", l[0].evaluate(&F::one()));
        // println!("b {:?}", b[0]);
        // println!("c {:?}", c[0]);

        let q = restrict_poly(l, Wi_1.clone());
        let q_0 = q.evaluate(&F::zero());
        let q_1 = q.evaluate(&F::one());

        // println!("q_0 {:?}", q_0);
        // println!("q_1 {:?}", q_1);
        // // println!("q_sum {:?}", q_0 + q_1);
        // println!("Wi_1_b {}", Polynomial::evaluate(&Wi_1, &b.to_vec()));
        // println!("Wi_1_c {}", Polynomial::evaluate(&Wi_1, &c.to_vec()));

        let final_poly_eval_prover = Polynomial::evaluate(&add_fixed, &used_r) * (q_0 + q_1) + Polynomial::evaluate(&mul_fixed, &used_r) * (q_0 * q_1);
        let final_poly_eval_prover = Polynomial::evaluate(&add_fixed, &used_r) * (Polynomial::evaluate(&Wi_1, &b.to_vec()) + Polynomial::evaluate(&Wi_1, &c.to_vec()))
            + Polynomial::evaluate(&mul_fixed, &used_r) * (Polynomial::evaluate(&Wi_1, &b.to_vec()) * Polynomial::evaluate(&Wi_1, &c.to_vec()));

        println!("final_poly_eval_prover {:?}", final_poly_eval_prover);
        // println!("interpolared {:?}", Polynomial::evaluate(&interpolated, &used_r));
        let last_poly = used_polys.last().unwrap();
        // println!("last poly {:?}", last_poly);
        println!("last poly eval {:?}", last_poly.evaluate(&used_r.last().unwrap()));
        // println!("inter last poly eval {}", used_poly_inter.last().unwrap().evaluate(&used_r.last().unwrap()));
        println!("last poly sum 01 {}\n", last_poly.evaluate(&F::zero()) + last_poly.evaluate(&F::one()));
    }
    
    // println!("W0 {:?}", W0);
    // println!("W0_0 {}", W0.evaluate(&vec![F::from(0), F::from(0)]));
    // println!("W0_0 {}", W0.evaluate(&vec![F::from(0), F::from(1)]));
    // println!("W0_0 {}", W0.evaluate(&vec![F::from(1), F::from(0)]));
    // println!("W0_0 {}", W0.evaluate(&vec![F::from(1), F::from(1)]));
    // println!("m0 {}", m0);
}

// l(t) = (1 - t) * b + t * c = b + t * (c - b)
fn line<F: Field>(b: &[F], c: &[F]) -> Vec<DensePolynomial<F>> {
    let mut polys = vec![];
    for (b, c) in b.iter().zip(c.iter()) {
        polys.push(DensePolynomial::from_coefficients_slice(&[*b, *c - *b]));
    }

    polys
}

// f(x_1, x_2) | l(t) = f_l(x_1(t), x_2(t))
fn restrict_poly<F: Field>(
    line: Vec<DensePolynomial<F>>,
    w: DenseMultilinearExtension<F>,
) -> DensePolynomial<F> {
    let num_vars = MultilinearExtension::num_vars(&w);
    assert_eq!(line.len(), num_vars);

    let line_rev = line.clone().into_iter().rev().collect::<Vec<_>>();

    let mut res = DensePolynomial::from_coefficients_slice(&[F::zero()]);

    let mut evals_to_map = w.get_evaluations();
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

    // println!("res {:?}", res);

    res
}

fn pad_with_zeroes<F: Field>(evaluations: &[F]) -> Vec<F> {
    // return evaluations.clone().to_vec();
    let n = (evaluations.len() as f64).log2().ceil() as usize;
    let padded_len = 1 << n;
    let mut evaluations = Vec::from(evaluations);
    evaluations.resize(padded_len, F::zero());
    
    evaluations
}

fn interpolate<F: Field>(evaluations: &[F]) -> DenseMultilinearExtension<F> {
    // let lagrange_polys = get_lagrange_polys(evaluations.len() as u64);
    // F::zero();

    DenseMultilinearExtension::from_evaluations_slice(
        evaluations.len().ilog2() as usize,
        evaluations,
    )
    
    
    // DenseMultilinearExtension::from_evaluations_vec(0, vec![])
}

fn get_lagrange_polys<F: Field>(evaluations_count: usize) -> Vec<DenseMultilinearExtension<F>> {
    let mut res = vec![];
    
    // for x in 0..params_count {
    //     for w in 0..params_count {
    //         
    //     }
    // }
    
    // let evaluations_count = 2_usize.pow(params_count as u32);
    let params_count = evaluations_count.ilog2();
    
    for i in 0..evaluations_count {
        // println!("i {}", i);
        // 
        // (0..params_count).for_each(|j| {
        //     // DenseMultilinearExtension::from_evaluations_vec()
        //     // println!("{}", (i >> j) & 1)
        // })
        res.push(DenseMultilinearExtension::from_evaluations_vec(
            params_count as usize,
            (0..evaluations_count).into_iter().map(|e| if i == e { F::one() } else { F::zero() }).collect(),
        ));
    }
    
    res
}

fn reverse_bits(mut n: usize, k: usize) -> usize {
    let mut res = 0;
    for _ in 0..k {
        res <<= 1;
        res |= n & 1;
        n >>= 1;
    }

    res
}

fn remap_to_reverse_bits_indexing<T: Debug>(vec: &mut [T], k: usize) {
    // return;

    for i in 0..vec.len() {
        let i_reversed = reverse_bits(i, k);
        if i_reversed > i {
            vec.swap(i, i_reversed);
        }
    }

    // let n = vec.len();
    // let half = n / 2;
    // let quarter = half / 2;
    //
    // for i in 0..half {
    //
    // }

    // for (i, e) in vec.iter().enumerate() {
    //     println!("i {} {:#08b} vec[i] {:?}", i, i, e);
    // }
    // println!();
    // for i in 0..vec.len() / 2 {
    //     let new_i = reverse_bits(i, k);
    //     println!("i {} {:#08b} new_i {} {:#08b} vec[i] {:?} vec[new_i] {:?}", i, i, new_i, new_i, vec[i], vec[new_i]);
    //     vec.swap(i, new_i);
    // }
    //
    // println!();
    //
    // for i in vec.len() / 2..vec.len() {
    //     if i % 2 == 1 {
    //         let new_i = reverse_bits(i, k);
    //         println!("i {} {:#08b} new_i {} {:#08b} vec[i] {:?} vec[new_i] {:?}", i, i, new_i, new_i, vec[i], vec[new_i]);
    //         vec.swap(i, new_i);
    //     }
    // }
}

// fn add_poly<F: Field>(a: &[F], b: &[F], c: &[F]) -> DenseMultilinearExtension<F> {
//
// }

// fn w_sum_over_last_variable<F: Field>(poly: RoundPoly<F>) -> RoundPoly<F> {
//     if poly.Wi_1_b.num_vars > 0 {
//         RoundPoly {
//             add_i: sum_over_last_variable(&poly.add_i),
//             mul_i: sum_over_last_variable(&poly.mul_i),
//             Wi_1_a: poly.Wi_1_a,
//             Wi_1_b: sum_over_last_variable(&poly.Wi_1_b),
//         }
//     } else {
//         RoundPoly {
//             add_i: sum_over_last_variable(&poly.add_i),
//             mul_i: sum_over_last_variable(&poly.mul_i),
//             Wi_1_a: sum_over_last_variable(&poly.Wi_1_a),
//             Wi_1_b: poly.Wi_1_b,
//         }
//     }
// }

fn interpolate_round_for_sc<F: Field>(
    add_i: DenseMultilinearExtension<F>,
    mul_i: DenseMultilinearExtension<F>,
    Wi_1: DenseMultilinearExtension<F>,
) -> DenseMultilinearExtension<F> {
    let num_vars = MultilinearExtension::num_vars(&add_i);
    let w_num_vars = MultilinearExtension::num_vars(&Wi_1);
    let varians = 1 << w_num_vars;

    let mut evals = vec![];

    for a in 0..varians {
        for b in 0..varians {
            let a_bits = get_bits(a, w_num_vars).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
            let b_bits = get_bits(b, w_num_vars).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
            // println!("a_bits {:?}", a_bits);
            // println!("b_bits {:?}", b_bits);

            let ab_bitmask = (a << w_num_vars) | b;
            let ab_bits = get_bits(ab_bitmask, num_vars).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();;

            let add_eval = Polynomial::evaluate(&add_i, &ab_bits) * (Polynomial::evaluate(&Wi_1, &a_bits) + Polynomial::evaluate(&Wi_1, &b_bits));
            let mul_eval = Polynomial::evaluate(&mul_i, &ab_bits) * (Polynomial::evaluate(&Wi_1, &a_bits) * Polynomial::evaluate(&Wi_1, &b_bits));

            evals.push(add_eval + mul_eval);
        }
    }

    remap_to_reverse_bits_indexing(&mut evals, num_vars);

    DenseMultilinearExtension::from_evaluations_vec(num_vars, evals)
}


#[cfg(test)]
mod tests {
    use super::*;

    fn test_3_var_poly() -> DenseMultilinearExtension<Fr> {
        let nvars = 3;
        let mut evals = [1,2,3,4,5,6,7,8].into_iter().map(|e| Fr::from(e as u64)).collect::<Vec<_>>();
        remap_to_reverse_bits_indexing(&mut evals, nvars);
        DenseMultilinearExtension::from_evaluations_vec(nvars, evals)
    }

    fn get_test_round_poly_2_vars<F: Field>() -> RoundPoly<F> {
        let add_i = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![0,0,0,0].into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>(),
        );
        let mul_i = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![0,0,1,0].into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>(),
        );
        let Wi_1 = DenseMultilinearExtension::from_evaluations_vec(
            1,
            vec![210, 320].into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>(),
        );
        RoundPoly {
            add_i,
            mul_i,
            Wi_1_a: Wi_1.clone(),
            Wi_1_b: Wi_1.clone(),
        }
    }

    fn get_test_round_poly_4_vars<F: Field>() -> RoundPoly<F> {
        let add_i = DenseMultilinearExtension::from_evaluations_vec(
            5,
            (0..32).into_iter().map(|e| F::from((e == 16 || e == 27) as u64)).collect::<Vec<_>>(),
        );
        let mul_i = DenseMultilinearExtension::from_evaluations_vec(
            5,
            (0..32).into_iter().map(|e| F::from(0)).collect::<Vec<_>>(),
        );
        // let add_i = ark_poly::MultilinearExtension::fix_variables(&add_i, &[F::from(3)]);
        // let mul_i = ark_poly::MultilinearExtension::fix_variables(&mul_i, &[F::from(3)]);
        let Wi_1 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![10,20,200,300].into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>(),
        );

        for i in 0..1 << mul_i.num_vars {
            let bits = get_bits(i, mul_i.num_vars).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();

            let eval = ark_poly::Polynomial::evaluate(&mul_i, &bits);

            // println!("bits {:?} res {:?}", bits, eval);
        }

        RoundPoly {
            add_i: MultilinearExtension::fix_variables(&add_i, &[F::from(3)]),
            mul_i: MultilinearExtension::fix_variables(&mul_i, &[F::from(3)]),
            Wi_1_a: Wi_1.clone(),
            Wi_1_b: Wi_1.clone(),
        }
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

    // [
    // 1, // 000
    // 2, // 001
    // 3, // 010
    // 4, // 011
    // 5, // 100
    // 6, // 101
    // 7, // 110
    // 8] // 111

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

            println!("res {:?}", res);
            assert_eq!(res.coeffs.len(), 2);
            assert_eq!(res.coeffs, [Fr::from(coef_1), Fr::from(coef_2)]);
        }
    }

    #[test]
    fn get_partial_sum_poly_2_vars() {
        let round_poly = get_test_round_poly_2_vars::<Fr>();
        let partial_sum_poly_1 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_1.coeffs, [Fr::from(67200), Fr::from(-32000), Fr::from(-35200)]);

        let round_poly = round_poly.fix_variables(&[Fr::from(32)]);
        let partial_sum_poly_2 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_2.coeffs, [Fr::zero(), Fr::from(-24282300), Fr::from(-12719300)])
    }

    #[test]
    fn get_partial_sum_poly_4_vars() {
        let round_poly = get_test_round_poly_4_vars::<Fr>();
        let partial_sum_poly_1 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_1.coeffs, [Fr::from(-420), Fr::from(1330), Fr::from(50)]);
        
        let round_poly = round_poly.fix_variables(&[Fr::from(4)]);
        let partial_sum_poly_2 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_2.coeffs, [Fr::from(5700), Fr::from(4200), Fr::from(-9900)]);

        let round_poly = round_poly.fix_variables(&[Fr::from(5)]);
        let partial_sum_poly_3 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_3.coeffs, [Fr::from(-72000), Fr::from(-74400), Fr::from(-2400)]);

        let round_poly = round_poly.fix_variables(&[Fr::from(8)]);
        let partial_sum_poly_4 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_4.coeffs, [Fr::from(0), Fr::from(-624240), Fr::from(-196560)]);
    }

    #[test]
    fn generate_round_poly_params() {
        let a_vars = 2;
        let b_vars = 2;

        let ((a_vars_fixing, b_vars_fixing), res) = generate_round_poly_eval_parameters(a_vars, b_vars);

        assert_eq!(res, [
            (0b000, 0b000, 0b000),
            (0b001, 0b000, 0b001),
            (0b010, 0b000, 0b010),
            (0b011, 0b000, 0b011),
            (0b100, 0b001, 0b000),
            (0b101, 0b001, 0b001),
            (0b110, 0b001, 0b010),
            (0b111, 0b001, 0b011),
        ]);
        assert_eq!(a_vars_fixing, 1);
        assert_eq!(b_vars_fixing, 2);
    }

    #[test]
    fn generate_round_poly_params_one_var() {
        let a_vars = 1;
        let b_vars = 2;

        let ((a_vars_fixing, b_vars_fixing), res) = generate_round_poly_eval_parameters(a_vars, b_vars);

        assert_eq!(res, [
            (0b000, 0b000, 0b000),
            (0b001, 0b000, 0b001),
            (0b010, 0b000, 0b010),
            (0b011, 0b000, 0b011),
        ]);
        assert_eq!(a_vars_fixing, 0);
        assert_eq!(b_vars_fixing, 2);
    }

    #[test]
    fn generate_round_poly_params_zero_var() {
        let a_vars = 0;
        let b_vars = 2;

        let ((a_vars_fixing, b_vars_fixing), res) = generate_round_poly_eval_parameters(a_vars, b_vars);

        assert_eq!(res, [
            (0b000, 0b000, 0b000),
            (0b001, 0b000, 0b001),
        ]);
        assert_eq!(a_vars_fixing, 0);
        assert_eq!(b_vars_fixing, 1);
    }

    #[test]
    fn round_poly_fix_polys() {
        let round_poly = get_test_round_poly_4_vars::<Fr>();
        let mask = 0b0001;
        let a_mask = 0b0000;
        let b_mask = 0b0001;

        let add_uni = to_two_or_one_degree(&round_poly.add_i, mask, 3);
        let mul_uni = to_two_or_one_degree(&round_poly.mul_i, mask, 3);
        let Wi_1_a_uni = to_two_or_one_degree(&round_poly.Wi_1_a, a_mask, 1);
        let Wi_1_b_uni = to_two_or_one_degree(&round_poly.Wi_1_b, b_mask, 2);

        assert_eq!(add_uni.coeffs, [Fr::from(-2), Fr::from(2)]);
        assert_eq!(mul_uni.coeffs, []);
        assert_eq!(Wi_1_a_uni.coeffs, [Fr::from(10), Fr::from(10)]);
        assert_eq!(Wi_1_b_uni.coeffs, [Fr::from(200)]);
    }

    #[test]
    fn round_poly_fix_polys_other_mask() {
        let round_poly = get_test_round_poly_4_vars::<Fr>();
        let mask = 0b0101;
        let a_mask = 0b0001;
        let b_mask = 0b0001;

        println!("round ad {:?}", round_poly.add_i.evaluations);
        println!("round ad(0101) - {:?}", ark_poly::Polynomial::evaluate(&round_poly.add_i, &vec![Fr::zero(), Fr::one(), Fr::zero(), Fr::one()]));
        println!("round ad(1101) - {:?}", ark_poly::Polynomial::evaluate(&round_poly.add_i, &vec![Fr::one(), Fr::one(), Fr::zero(), Fr::one()]));
        let add_uni = to_two_or_one_degree(&round_poly.add_i, mask, 3);
        let mul_uni = to_two_or_one_degree(&round_poly.mul_i, mask, 3);
        let Wi_1_a_uni = to_two_or_one_degree(&round_poly.Wi_1_a, a_mask, 1);
        let Wi_1_b_uni = to_two_or_one_degree(&round_poly.Wi_1_b, b_mask, 2);

        assert_eq!(add_uni.coeffs, []);
        assert_eq!(mul_uni.coeffs, []);
        assert_eq!(Wi_1_a_uni.coeffs, [Fr::from(200), Fr::from(100)]);
        assert_eq!(Wi_1_b_uni.coeffs, [Fr::from(200)]);
    }

    #[test]
    fn round_poly_get_required_params() {
        let round_poly = get_test_round_poly_4_vars::<Fr>();
        let ((a_eval_params, b_eval_params), evals) = generate_round_poly_eval_parameters(
            round_poly.Wi_1_a.num_vars,
            round_poly.Wi_1_b.num_vars,
        );

        assert_eq!(a_eval_params, 1);
        assert_eq!(b_eval_params, 2);
        assert_eq!(evals[1], (0b0001usize, 0b0000usize, 0b0001usize));
    }
}