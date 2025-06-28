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

fn get_reversed_vars_poly<F: Field>(poly: &DenseMultilinearExtension<F>) -> DenseMultilinearExtension<F> {
    let n = MultilinearExtension::num_vars(poly);
    let mut res_poly = poly.clone();

    for i in 0..n / 2 {
        res_poly.relabel_in_place(i, n - i - 1,1);
    }

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

impl<F: Field> Layer<'_, F> {
    fn get_gates_n(&self) -> usize {
        self.gates.len()
    }
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
        let (a_vars_fixing, b_vars_fixing) = {
            if a_vars > 0 {
                (a_vars - 1, b_vars)
            } else {
                (a_vars, b_vars - 1)
            }
        };
        let a_domain = 1 << a_vars_fixing;
        let b_domain = 1 << b_vars_fixing;
        let Wi_1a_rev = get_reversed_vars_poly(&self.Wi_1_a);
        let Wi_1b_rev = get_reversed_vars_poly(&self.Wi_1_b);
        let add_rev = get_reversed_vars_poly(&self.add_i);
        let mul_rev = get_reversed_vars_poly(&self.mul_i);

        let mut res_poly = DensePolynomial::from_coefficients_vec(vec![F::zero()]);

        for ai in 0..a_domain {
            for bi in 0..b_domain {
                let a_bits = get_bits(ai, a_vars_fixing).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
                let b_bits = get_bits(bi, b_vars_fixing).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
                let mut ab_bits = a_bits.clone();
                ab_bits.extend(b_bits.clone());

                let Wi_a_fixed = MultilinearExtension::fix_variables(&Wi_1a_rev, a_bits.as_slice());
                let Wi_b_fixed = MultilinearExtension::fix_variables(&Wi_1b_rev, b_bits.as_slice());
                let add_fixed = MultilinearExtension::fix_variables(&add_rev, ab_bits.as_slice());
                let mul_fixed = MultilinearExtension::fix_variables(&mul_rev, ab_bits.as_slice());

                let Wi_a_dense = interpolate_or_const(&Wi_a_fixed.to_evaluations());
                let Wi_b_dense = interpolate_or_const(&Wi_b_fixed.to_evaluations());
                let add_dense = interpolate_or_const(&add_fixed.to_evaluations());
                let mul_dense = interpolate_or_const(&mul_fixed.to_evaluations());

                let add_poly = add_dense.naive_mul(&(&Wi_a_dense + &Wi_b_dense));
                let mul_poly = mul_dense.naive_mul(&(Wi_a_dense.naive_mul(&Wi_b_dense)));

                res_poly = res_poly + add_poly + mul_poly;
            }
        }

        // println!("res poly new inter {:?}", res_poly);
        let res_poly_0 = res_poly.evaluate(&F::zero());
        let res_poly_1 = res_poly.evaluate(&F::one());
        // println!("res poly(0) {}", res_poly_0);
        // println!("res poly(1) {}", res_poly_1);
        // println!("res poly sum {}", res_poly_0 + res_poly_1);

        return res_poly;


        let mut current_sc_poly = self.clone();
        let num_vars = current_sc_poly.num_vars();

        let evals = self.get_evaluations();
        println!("evals {:?}", evals);

        let (zeroes, ones) = evals.iter().enumerate().fold((F::zero(), F::zero()), |(zeroes, ones), (i, eval)| {
            if i % 2 == 0 {
                (zeroes + eval, ones)
            } else {
                (zeroes, ones + eval)
            }
        });

        interpolate_univariate_on_evals(&[zeroes, ones])
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
            RoundPoly {
                add_i: MultilinearExtension::fix_variables(&self.add_i, partial_point),
                mul_i: MultilinearExtension::fix_variables(&self.mul_i, partial_point),
                Wi_1_a: MultilinearExtension::fix_variables(&self.Wi_1_a, partial_point),
                Wi_1_b: self.Wi_1_b.clone(),
            }
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
        value: Fr::from(9)
    };
    let input_b = InputGate {
        index: 1,
        value: Fr::one()
    };
    let input_c = InputGate {
        index: 2,
        value: Fr::from(2),
    };

    let add_gate = Gate {
        inputs: &GateInput::Inputs((&input_a, &input_b)),
        executor: ExecutorGateEnum::Add(AddGate {})
    };
    let mul_gate = Gate {
        inputs: &GateInput::Inputs((&input_c, &input_c)),
        executor: ExecutorGateEnum::Mul(MulGate {})
    };
    let last_mul_gate = Gate {
        inputs: &GateInput::Gates((&add_gate, &mul_gate)),
        executor: ExecutorGateEnum::Mul(MulGate {})
    };

    let layer_1 = Layer {
        gates: vec![&add_gate, &mul_gate],
    };
    let layer_2 = Layer {
        gates: vec![&last_mul_gate],
    };

    let circuit = Circuit {
        layers: vec![layer_1, layer_2],
    };

    let solution = circuit.solve();
    
    // println!("sol {:?}", solution);
    
    let test_poly = DenseMultilinearExtension::from_evaluations_vec(
        3,
        vec![1,2,3,4,5,6,7,8].into_iter().map(Fr::from).collect(),
    );
    
    // let test_poly = get_lagrange_polys(8);
    
    // for i in 0..2_usize.pow(3) {
    //     let bits = (0..3).map(|b| Fr::from((i >> b & 1) as u64)).collect::<Vec<_>>();
    //     let res = test_poly[i].evaluate(bits.as_ref());
    //     let res_1 = test_poly[(i + 1) % 8].evaluate(bits.as_ref());
    //
    //     println!("i {} bits {:?} evals {:?}", i, bits, test_poly[i].evaluations);
    // }
    
    let random_points = vec![585,22,12,93,8181,12398,123]
        .into_iter()
        .map(Fr::from)
        .collect::<Vec<_>>();

    let test_w_poly = DenseMultilinearExtension::from_evaluations_vec(2, vec![
        Fr::from(0),
        Fr::from(0),
        Fr::from(2),
        Fr::from(5),
    ]);
    // println!("test_w_poly {:?}", test_w_poly);
    // println!("eval test_w_poly {}", Polynomial::evaluate(&test_w_poly, &vec![Fr::from(0), Fr::from(0)]));
    // println!("eval test_w_poly {}", Polynomial::evaluate(&test_w_poly, &vec![Fr::from(0), Fr::from(1)]));
    // println!("eval test_w_poly {}", Polynomial::evaluate(&test_w_poly, &vec![Fr::from(1), Fr::from(0)]));
    // println!("eval test_w_poly {}", Polynomial::evaluate(&test_w_poly, &vec![Fr::from(1), Fr::from(1)]));

    let b_test = vec![Fr::from(2), Fr::from(4)];
    let c_test = vec![Fr::from(3), Fr::from(2)];

    let line_test = line(&b_test, &c_test);
    // println!("line_test {:?}", line_test);

    let w_restricted = restrict_poly(line_test, test_w_poly.clone());
    // println!("w_restricted {:?}", w_restricted);
    // println!("W_test_b {}", Polynomial::evaluate(&test_w_poly, &b_test));
    // println!("W_test_c {}", Polynomial::evaluate(&test_w_poly, &c_test));
    // println!("test_q_0 {}", w_restricted.evaluate(&Fr::from(0)));
    // println!("test_q_1 {}", w_restricted.evaluate(&Fr::from(1)));

    prove(circuit, solution, &random_points);
    let bits = get_bits(3, 0b110);
    // println!("bits {:?}", bits);
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
                interpolate(&pad_with_zeroes(&solution.inputs))
            } else {
                interpolate(&solution.evaluations[i - 1])
            }
        };

        let add_fixed = MultilinearExtension::fix_variables(&add_poly, &r0);
        let mul_fixed = MultilinearExtension::fix_variables(&mul_poly, &r0);
        // println!("r0 {:?}", r0);
        // println!("Wi_1 {:?}", Wi_1);

        let interpolated = interpolate_round_for_sc(add_fixed.clone(), mul_fixed.clone(), Wi_1.clone());
        let sc_poly = RoundPoly {
            add_i: add_fixed.clone(),
            mul_i: mul_fixed.clone(),
            Wi_1_a: Wi_1.clone(),
            Wi_1_b: Wi_1.clone(),
        };
        // println!("sc eval {}", sc_poly.evaluate(&[F::zero(), F::one()]));

        // println!("m0 {:?}", m0);
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
    for i in 0..vec.len() / 2 {
        let new_i = reverse_bits(i, k);
        vec.swap(i, new_i);
    }
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
