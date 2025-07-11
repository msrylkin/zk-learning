use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::ops::{Add, Mul, MulAssign, Neg};
use ark_ff::{Field, One, Zero};
use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_test_curves::bls12_381::Fr;
use crate::sumcheck;
use crate::sumcheck::{interpolate_univariate_on_evals, OracleEvaluation, SumCheckPoly, SumCheckProof};

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

impl<F: Field> From<Solution<F>> for Vec<Vec<F>> {
    fn from(solution: Solution<F>) -> Self {
        let mut res = vec![solution.inputs];
        res.extend_from_slice(&solution.evaluations);

        res
    }
}

#[derive(Debug, Clone)]
struct LayerRoundPoly<F: Field> {
    add_i: DenseMultilinearExtension<F>,
    mul_i: DenseMultilinearExtension<F>,
    Wi_1_a: DenseMultilinearExtension<F>,
    Wi_1_b: DenseMultilinearExtension<F>,
}

#[derive(Debug, Clone)]
struct GKRFinalOracle<F: Field> {
    mul_i: DenseMultilinearExtension<F>,
    add_i: DenseMultilinearExtension<F>,
    q: DensePolynomial<F>
}

impl<F: Field> GKRFinalOracle<F> {
    fn new(add_i: DenseMultilinearExtension<F>, mul_i: DenseMultilinearExtension<F>, q: DensePolynomial<F>) -> Self {
        Self { mul_i, add_i, q }
    }
}

impl<F: Field> OracleEvaluation<F> for GKRFinalOracle<F> {
    fn final_eval(&self, r: &[F]) -> F {
        let q_0 = self.q.evaluate(&F::zero());
        let q_1 = self.q.evaluate(&F::one());
        let add_eval = ark_poly::Polynomial::evaluate(&self.add_i, &r.to_vec()) * (q_0 + q_1);
        let mul_eval = ark_poly::Polynomial::evaluate(&self.mul_i, &r.to_vec()) * (q_0 * q_1);
        
        add_eval + mul_eval
    }
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

fn to_f<F: Field>(arr: Vec<u64>) -> Vec<F> {
    arr.into_iter().map(F::from).collect()
}

fn interpolate_or_const<F: Field>(evals: &[F]) -> DensePolynomial<F> {
    match evals.len() {
        1 => DensePolynomial::from_coefficients_slice(evals),
        2 => interpolate_univariate_on_evals(evals.try_into().unwrap()),
        _ => panic!("wrong evaluations len"),
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

    if nvars != MultilinearExtension::num_vars(poly) - 1 {
        todo!("only single variable supported for now");
    }

    let bit_set = 1 << nvars;
    let bit_set_negated = !bit_set;
    let mut mle_evals = poly.to_evaluations();
    remap_to_reverse_bits_indexing(&mut mle_evals, MultilinearExtension::num_vars(poly));

    let index_1 = bit_set_negated & mask;
    let index_2 = bit_set | mask;

    let mut new_evals = vec![];

    new_evals.push(mle_evals[index_1]);
    new_evals.push(mle_evals[index_2]);

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

    let mut res = vec![];

    for i in 0..combined_domain {
        let combined_param = i;
        let b_param = i & b_mask;
        let a_param = (i & a_mask) >> b_vars_fixing;

        res.push((combined_param, a_param, b_param))
    }

    ((a_vars_fixing, b_vars_fixing), res)
}

impl<F: Field> Layer<'_, F> {
    fn get_gates_n(&self) -> usize {
        self.gates.len()
    }
}

impl<F: Field> SumCheckPoly<F> for LayerRoundPoly<F> {
    fn get_evaluations(&self) -> Vec<F> {
        let mut res = vec![];
        let a_num_vars = MultilinearExtension::num_vars(&self.Wi_1_a);
        let b_num_vars = MultilinearExtension::num_vars(&self.Wi_1_b);
        let ab_num_vars = a_num_vars + b_num_vars;
        let a_domain = 1 << a_num_vars;
        let b_domain = 1 << b_num_vars;

        for a in 0..a_domain {
            for b in 0..b_domain {
                let ab_combined_mask = (a << b_num_vars) | b;
                let ab_combined_mask_reversed = reverse_bits(ab_combined_mask, ab_num_vars);
                let ab_bits = get_bits(ab_combined_mask_reversed, ab_num_vars).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
                let eval = LayerRoundPoly::evaluate(&self, &ab_bits);

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

        let mut res_poly = DensePolynomial::from_coefficients_vec(vec![F::zero()]);

        let ((a_vars_fixing, b_vars_fixing), eval_params) = generate_round_poly_eval_parameters(a_vars, b_vars);

        for (common_params, a_params, b_params) in eval_params {
            let Wi_a_uni = to_two_or_one_degree(&self.Wi_1_a, a_params, a_vars_fixing);
            let Wi_b_uni = to_two_or_one_degree(&self.Wi_1_b, b_params, b_vars_fixing);
            let add_i_uni = to_two_or_one_degree(&self.add_i, common_params, a_vars_fixing + b_vars_fixing);
            let mul_uni = to_two_or_one_degree(&self.mul_i, common_params, a_vars_fixing + b_vars_fixing);

            let add_poly = add_i_uni.naive_mul(&(&Wi_a_uni + &Wi_b_uni));
            let mul_poly = mul_uni.naive_mul(&(Wi_a_uni.naive_mul(&Wi_b_uni)));
            let local_res_poly = add_poly + mul_poly;

            res_poly = res_poly + local_res_poly;
        }

        res_poly
    }

    fn fix_variables(&self, partial_point: &[F]) -> Self {
        if MultilinearExtension::num_vars(&self.Wi_1_a) > 0 {
            Self {
                add_i: MultilinearExtension::fix_variables(&self.add_i, partial_point),
                mul_i: MultilinearExtension::fix_variables(&self.mul_i, partial_point),
                Wi_1_a: MultilinearExtension::fix_variables(&self.Wi_1_a, partial_point),
                Wi_1_b: self.Wi_1_b.clone(),
            }
        } else {
            Self {
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

        let add_eval = Polynomial::evaluate(&self.add_i, &point.to_vec())
            * (Polynomial::evaluate(&self.Wi_1_a, &a_bits) + Polynomial::evaluate(&self.Wi_1_b, &b_bits));
        let mul_eval = Polynomial::evaluate(&self.mul_i, &point.to_vec())
            * (Polynomial::evaluate(&self.Wi_1_a, &a_bits) * Polynomial::evaluate(&self.Wi_1_b, &b_bits));

        add_eval + mul_eval
    }
}

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
        let vars_num = f64::from(gates_n as u32).log2().ceil() as usize;
        let domain_size = 1 << vars_num;
        let mut bottom_gates_n = self.gates_n_at_layer_i_1(layer_i);
        if !bottom_gates_n.is_power_of_two() {
            bottom_gates_n = bottom_gates_n.next_power_of_two();
        }
        let bottom_vars_num = f64::from(bottom_gates_n as u32).log2().ceil() as usize;
        let bottom_domain_size = 1 << bottom_vars_num;

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

        DenseMultilinearExtension::from_evaluations_vec(vars_num + 2 * bottom_vars_num, evals)
    }

    fn add_i(&self, layer_i: usize) -> DenseMultilinearExtension<F> {
        let layer = &self.layers[layer_i];
        let mut gates_n = layer.gates.len();
        if !gates_n.is_power_of_two() {
            gates_n = gates_n.next_power_of_two();
        }
        let vars_num = f64::from(gates_n as u32).log2().ceil() as usize;
        let domain_size = 1 << vars_num;
        let mut bottom_gates_n = self.gates_n_at_layer_i_1(layer_i);
        if !bottom_gates_n.is_power_of_two() {
            bottom_gates_n = bottom_gates_n.next_power_of_two();
        }
        let bottom_vars_num = f64::from(bottom_gates_n as u32).log2().ceil() as usize;
        let bottom_domain_size = 1 << bottom_vars_num;

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

        DenseMultilinearExtension::from_evaluations_vec(vars_num + 2 * bottom_vars_num, evals)
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
    
    let random_points = vec![3,22,12,93,8181,12398,123]
        .into_iter()
        .map(Fr::from)
        .collect::<Vec<_>>();

    let gkr_proof = prove(&circuit, &solution, &random_points);
    verify(circuit, &gkr_proof);
}

fn verify<F: Field>(
    circuit: Circuit<F>,
    proof: &GKRProof<F>,
) {
    let mut mi = Polynomial::evaluate(&proof.W0, &proof.r0);
    let mut ri = proof.r0.clone();

    for (i, GKRProofLayer { q, r_star, sumcheck_proof }) in proof.layers.iter().enumerate().rev() {
        let add_poly = circuit.add_i(i);
        let mul_poly = circuit.mul_i(i);
        let add_fixed = MultilinearExtension::fix_variables(&add_poly, &ri);
        let mul_fixed = MultilinearExtension::fix_variables(&mul_poly, &ri);
        let used_r = sumcheck_proof.get_used_randomness();
        let (b, c) = used_r.split_at(used_r.len() / 2);
        let l = line(b, c);

        let final_oracle = GKRFinalOracle::new(add_fixed, mul_fixed, q.clone());

        sumcheck::verify(&final_oracle, &sumcheck_proof, mi);

        ri = l.iter().map(|li| li.evaluate(r_star)).collect::<Vec<_>>();
        mi = q.evaluate(r_star);
    }

    let last_w = interpolate(&proof.inputs);
    let last_w_r = Polynomial::evaluate(&last_w, &ri.to_vec());

    assert_eq!(last_w_r, mi);
}

#[derive(Debug)]
struct GKRProofLayer<F: Field> {
    sumcheck_proof: SumCheckProof<F>,
    r_star: F,
    q: DensePolynomial<F>,
}

#[derive(Debug)]
struct GKRProof<F: Field> {
    layers: Vec<GKRProofLayer<F>>,
    outputs: Vec<F>,
    inputs: Vec<F>,
    r0: Vec<F>,
    W0: DenseMultilinearExtension<F>,
}

fn prove<F: Field>(
    circuit: &Circuit<F>,
    solution: &Solution<F>,
    random_points: &[F],
) -> GKRProof<F> {
    let outputs = solution.evaluations.last().unwrap();
    let W0 = interpolate(&outputs);
    let mut spent_points = MultilinearExtension::num_vars(&W0);
    let r0 = random_points[..spent_points].to_vec();
    let mut ri = r0.clone();

    let mut gkr_proof_layers = vec![];

    for i in (0..solution.evaluations.len()).rev() {
        let add_poly = circuit.add_i(i);
        let mul_poly = circuit.mul_i(i);
        let Wi_1 = {
            if i == 0 {
                interpolate(&solution.inputs)
            } else {
                interpolate(&solution.evaluations[i - 1])
            }
        };

        let add_fixed = MultilinearExtension::fix_variables(&add_poly, &ri);
        let mul_fixed = MultilinearExtension::fix_variables(&mul_poly, &ri);

        let sc_poly = LayerRoundPoly {
            add_i: add_fixed.clone(),
            mul_i: mul_fixed.clone(),
            Wi_1_a: Wi_1.clone(),
            Wi_1_b: Wi_1.clone(),
        };

        let sumcheck_proof = sumcheck::prove(&sc_poly);

        let used_r = sumcheck_proof.get_used_randomness();

        let (b, c) = used_r.split_at(used_r.len() / 2);
        let l = line(b, c);

        let q = restrict_poly(&l, Wi_1);
        let r_star = random_points[spent_points];
        spent_points += 1;

        ri = l.iter().map(|li| li.evaluate(&r_star)).collect();

        gkr_proof_layers.push(GKRProofLayer {
            q,
            sumcheck_proof,
            r_star,
        });
    }

    GKRProof {
        W0,
        outputs: outputs.clone(),
        inputs: solution.inputs.clone(),
        layers: gkr_proof_layers.into_iter().rev().collect(),
        r0,
    }
}

// l(t) = (1 - t) * b + t * c = b + t * (c - b)
fn line<F: Field>(b: &[F], c: &[F]) -> Vec<DensePolynomial<F>> {
    assert_eq!(b.len(), c.len());

    let mut polys = vec![];
    for (b, c) in b.iter().zip(c.iter()) {
        polys.push(DensePolynomial::from_coefficients_slice(&[*b, *c - *b]));
    }

    polys
}

// f(x_1, x_2) | l(t) = f_l(x_1(t), x_2(t))
fn restrict_poly<F: Field>(
    line: &Vec<DensePolynomial<F>>,
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

    res
}

fn pad_with_zeroes<F: Field>(evaluations: &[F]) -> Vec<F> {
    let n = (evaluations.len() as f64).log2().ceil() as usize;
    let padded_len = 1 << n;
    let mut evaluations = Vec::from(evaluations);
    evaluations.resize(padded_len, F::zero());
    
    evaluations
}

fn interpolate<F: Field>(evaluations: &[F]) -> DenseMultilinearExtension<F> {
    let mut evaluations = pad_with_zeroes(evaluations);
    let bitlen = evaluations.len().ilog2() as usize;
    // TODO: maybe get rid of this and create solution object directly with reversed bits order
    remap_to_reverse_bits_indexing(&mut evaluations, bitlen);

    DenseMultilinearExtension::from_evaluations_slice(
        bitlen,
        &evaluations,
    )
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
    for i in 0..vec.len() {
        let i_reversed = reverse_bits(i, k);
        if i_reversed > i {
            vec.swap(i, i_reversed);
        }
    }
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

    fn get_test_round_poly_2_vars<F: Field>() -> LayerRoundPoly<F> {
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
        LayerRoundPoly {
            add_i,
            mul_i,
            Wi_1_a: Wi_1.clone(),
            Wi_1_b: Wi_1.clone(),
        }
    }

    fn get_test_round_poly_4_vars<F: Field>() -> LayerRoundPoly<F> {
        let add_i = DenseMultilinearExtension::from_evaluations_vec(
            5,
            (0..32).into_iter().map(|e| F::from((e == 16 || e == 27) as u64)).collect::<Vec<_>>(),
        );
        let mul_i = DenseMultilinearExtension::from_evaluations_vec(
            5,
            (0..32).into_iter().map(|e| F::from(0)).collect::<Vec<_>>(),
        );
        let Wi_1 = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![10,20,200,300].into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>(),
        );

        LayerRoundPoly {
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

    #[test]
    fn line_test() {
        let b = &[Fr::from(3), Fr::from(7)];
        let c = &[Fr::from(4), Fr::from(5)];
        let l = line(b, c);

        assert_eq!(l.len(), 2);
        assert_eq!(l[0].coeffs, [Fr::from(3), Fr::one()]);
        assert_eq!(l[1].coeffs, [Fr::from(7), Fr::from(-2)]);
    }

    #[test]
    fn restrict_to_line_test() {
        let b = &[Fr::from(3), Fr::from(7)];
        let c = &[Fr::from(4), Fr::from(5)];
        let l = line(b, c);
        let w = DenseMultilinearExtension::from_evaluations_vec(
            2,
            vec![10, 20, 200, 300].into_iter().map(Fr::from).collect(),
        );

        let w_restricted = restrict_poly(&l, w.clone());
        let r_star = Fr::from(12);

        let l_evals = l.iter().map(|li| li.evaluate(&r_star)).collect::<Vec<_>>();

        assert_eq!(l_evals.len(), 2);
        assert_eq!(l_evals[0], Fr::from(15));
        assert_eq!(l_evals[1], Fr::from(-17));

        let w_r_star = ark_poly::Polynomial::evaluate(&w, &l_evals);
        assert_eq!(w_r_star, Fr::from(-26020));
        
        let w_restricted_r_star = w_restricted.evaluate(&r_star);
        assert_eq!(w_restricted_r_star, w_r_star);
    }
    
    #[test]
    fn add_i() {
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

        let layer_1 = Layer {
            gates: vec![&add_gate_1, &add_gate_2],
        };
        let layer_2 = Layer {
            gates: vec![&mul_gate],
        };

        let circuit = Circuit {
            layers: vec![layer_1, layer_2],
        };

        let add_i = circuit.add_i(0);
        println!("add_i evals {:?}", add_i.evaluations);
        println!("add_i varsnum {}", add_i.num_vars);
        println!("add_i(0,0,0,0,1) {}", ark_poly::Polynomial::evaluate(&add_i, &vec![Fr::zero(), Fr::zero(), Fr::zero(), Fr::zero(), Fr::one()]));
        println!("add_i(1,1,0,1,1) {}", ark_poly::Polynomial::evaluate(&add_i, &vec![Fr::one(), Fr::one(), Fr::zero(), Fr::one(), Fr::one()]));
        let add_i = circuit.add_i(1);
        println!("add_i evals {:?}", add_i.evaluations);
        println!("add_i varsnum {}", add_i.num_vars);

        let mul_i = circuit.mul_i(0);
        println!("mul_i evals {:?}", mul_i.evaluations);
        println!("mul_i varsnum {}", mul_i.num_vars);
        let mul_i = circuit.mul_i(1);
        println!("mul_i evals {:?}", mul_i.evaluations);
        println!("mul_i varsnum {}", mul_i.num_vars);
    }
}