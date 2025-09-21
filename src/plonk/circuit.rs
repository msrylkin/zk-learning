use std::collections::HashMap;
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::Radix2EvaluationDomain;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use ark_std::iterable::Iterable;
use ark_std::Zero;
use crate::evaluation_domain::MultiplicativeSubgroup;
use crate::plonk::prover::pick_coset_shifters;
// struct CircuitBuilder<F: FftField + PrimeField> {
//     circuit: Circuit<F>,
//     vars: HashMap<usize, Variable>,
//     last_var: usize,
// }

type DP<F> = DensePolynomial<F>;

enum Input<F: FftField + PrimeField> {
    Public(Variable),
    Witness(Variable),
    Constant(F),
    InnerWitness,
}

type Variable = (usize, usize); // (col, row)

// impl<F: PrimeField + FftField> CircuitBuilder<F> {
//     fn new() -> Self {
//         Self {
//             circuit: Circuit {
//                 gates: Vec::new(),
//             },
//             vars: HashMap::new(),
//             last_var: 0,
//         }
//     }
//
//     fn create_var() {}
//
//     // fn create_public_var(&mut self) -> usize {
//     //     // let var_id = self.last_var;
//     //     // let col = self.circuit.gates.len() - 1;
//     //     // self.circuit.gates.push(CircuitGate {
//     //     //     left: Input::Public((col, 0)),
//     //     //     right: Input::
//     //     // });
//     //     // self.last_var += 1;
//     //     //
//     //     // var_id
//     // }
//
//     // fn add_gate(
//     //     self,
//     //     left: Input<F>,
//     //     right: Input<F>,
//     //     out: Input<F>,
//     // ) -> Self {
//     //     let circuit = self.circuit;
//     //     let mut gates = circuit.gates;
//     //
//     //
//     // }
//
//     fn mul_gate() {
//
//     }
// }

enum CircuitOperation {
    Mul,
    Add,
}

struct CircuitGate<F: FftField + PrimeField> {
    // operation: CircuitOperation,
    left: Input<F>,
    right: Input<F>,
    output: Input<F>,
    ql: bool,
    qr: bool,
    qm: bool,
    qo: bool,
    qc: Input<F>,
}

struct CompiledGate<F: FftField + PrimeField> {
    // operation: CircuitOperation,
    left: F,
    right: F,
    output: F,
    ql: bool,
    qr: bool,
    qm: bool,
    qo: bool,
    qc: F,
}

pub struct CompiledCircuit<F: FftField + PrimeField> {
    gates: Vec<CompiledGate<F>>,
    sigma: Vec<usize>,
    public_inputs_count: usize,
}

#[derive(Clone, Debug)]
pub struct Solution<F: FftField + PrimeField> {
    pub a: DensePolynomial<F>,
    pub b: DensePolynomial<F>,
    pub c: DensePolynomial<F>,
    pub ql: DensePolynomial<F>,
    pub qr: DensePolynomial<F>,
    pub qm: DensePolynomial<F>,
    pub qo: DensePolynomial<F>,
    pub qc: DensePolynomial<F>,
    pub sid_1: DensePolynomial<F>,
    pub sid_2: DensePolynomial<F>,
    pub sid_3: DensePolynomial<F>,
    pub s_sigma_1: DensePolynomial<F>,
    pub s_sigma_2: DensePolynomial<F>,
    pub s_sigma_3: DensePolynomial<F>,
    pub pi: DensePolynomial<F>,
}

// pub struct Circuit<F: FftField + PrimeField> {
//
// }
//
// impl<F: FftField + PrimeField> CircuitBuilder<F> {
//
// }

fn format_bool<F: Field>(b: bool) -> F {
    match b {
        true => F::one(),
        false => F::zero(),
    }
}

// struct MakedCircuit<F: FftField + PrimeField> {
//     ql:
// }

fn compile_gates_to_circuit<F: FftField + PrimeField>(
    public_inputs_count: usize,
    gates: &[CompiledGate<F>],
) {

}

impl<F: FftField + PrimeField> CompiledCircuit<F> {
    pub fn get_abc_vectors(&self, domain: &[F]) -> (Vec<F>, Vec<F>, Vec<F>) {
        let mut a = vec![];
        let mut b = vec![];
        let mut c = vec![];

        for gate in &self.gates {
            a.push(gate.left);
            b.push(gate.right);
            c.push(gate.output);
        }

        for _ in self.gates.len()..domain.len() {
            a.push(F::zero());
            b.push(F::zero());
            c.push(F::zero());
        }

        // let padded_n = self.padded_len();
        // let padded_n =
        // for _ in 0..padded_n - self.gates.len() {
        //     a.push(F::zero());
        //     b.push(F::zero());
        //     c.push(F::zero());
        // }

        (a, b, c)
    }

    // pub fn padded_len(&self) -> usize {
    //     // self
    //     // let n = self.gates.len();
    //     //
    //     // match n.is_power_of_two() {
    //     //     true => n,
    //     //     false => n.next_power_of_two(),
    //     // }
    // }

    pub fn get_sigma(&self) -> Vec<usize> {
        self.sigma.clone()
    }

    pub fn get_selectors(&self, domain: &MultiplicativeSubgroup<F>) -> (
        DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>
    ) {
        let mut ql_values = vec![];
        let mut qr_values = vec![];
        let mut qm_values = vec![];
        let mut qo_values = vec![];
        let mut qc_values = vec![];

        for gate in &self.gates {
            ql_values.push(format_bool(gate.ql));
            qr_values.push(format_bool(gate.qr));
            qm_values.push(format_bool(gate.qm));
            qo_values.push(-format_bool::<F>(gate.qo));
            qc_values.push(gate.qc);
        }

        // let padded_n = self.padded_len();
        let padded_n = domain.len();
        for _ in 0..padded_n - self.gates.len() {
            ql_values.push(format_bool(false));
            qr_values.push(format_bool(false));
            qm_values.push(format_bool(false));
            qo_values.push(format_bool(false));
            qc_values.push(F::zero());
        }

        let ql = domain.interpolate_univariate(&ql_values);
        let qr = domain.interpolate_univariate(&qr_values);
        let qm = domain.interpolate_univariate(&qm_values);
        let qo = domain.interpolate_univariate(&qo_values);
        let qc = domain.interpolate_univariate(&qc_values);

        (ql, qr, qm , qo, qc)
    }

    fn sigma(&self, from: usize) -> usize {
        self.sigma[from]
    }

    pub fn get_sigma_polys(&self, domain: &MultiplicativeSubgroup<F>, k1: F, k2: F) -> (DP<F>, DP<F>, DP<F>) {
        (
            self.s_sigma_poly(1, domain, k1, k2),
            self.s_sigma_poly(2, domain, k1, k2),
            self.s_sigma_poly(3, domain, k1, k2),
        )
    }

    pub fn get_s_id_polys(&self, domain: &MultiplicativeSubgroup<F>, k1: F, k2: F) -> (DP<F>, DP<F>, DP<F>) {
        (
            self.s_id_poly(1, domain, k1, k2),
            self.s_id_poly(2, domain, k1, k2),
            self.s_id_poly(3, domain, k1, k2),
        )
    }

    fn s_sigma_poly(&self, col: usize, domain: &MultiplicativeSubgroup<F>, k1: F, k2: F) -> DensePolynomial<F> {
        let mut values = vec![];
        let n = self.gates.len();

        for i in 0..n {
            let permuted_index = self.sigma((col - 1) * n + i);
            let permuted_omega = self.map_index_to_coset_value(permuted_index, domain, k1, k2);
            values.push(permuted_omega);
        }

        // let padded_n = self.padded_len();
        let padded_n = domain.len();
        for i in self.gates.len()..padded_n {
            values.push(self.get_coset_shifter(col, k1, k2) * domain[i]);
        }

        domain.interpolate_univariate(&values)
    }

    fn s_id_poly(&self, col: usize, domain: &MultiplicativeSubgroup<F>, k1: F, k2: F) -> DP<F> {
        let mut values = vec![];
        // let n = self.padded_len();
        let n = domain.len();

        for i in 0..n {
            let permuted_omega = self.get_coset_shifter(col, k1, k2) * domain[i];
            values.push(permuted_omega);
        }

        domain.interpolate_univariate(&values)
    }

    fn map_index_to_coset_value(&self, index: usize, domain: &[F], k1: F, k2: F) -> F {
        let n = self.gates.len();

        match index {
            index if index < n => domain[index],
            index if index >= n && index < 2 * n => k1 * domain[index - n],
            _ => k2 * domain[index - 2 * n]
        }
    }

    fn get_coset_shifter(&self, col: usize, k1: F, k2: F) -> F {
        match col {
            1 => F::one(),
            2 => k1,
            3 => k2,
            _ => panic!("invalid col {}", col),
        }
    }

    pub fn get_public_input_poly(&self, domain: &MultiplicativeSubgroup<F>) -> DensePolynomial<F> {
        let values = self
            .get_abc_vectors(domain).0
            .iter()
            .enumerate()
            .map(|(i, e)| {
                if i < self.public_inputs_count {
                    return -*e;
                }

                F::zero()
            })
            .collect::<Vec<_>>();

        domain.interpolate_univariate(&values)
    }

    fn get_abc_polys(&self, domain: &MultiplicativeSubgroup<F>) -> (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>) {
        let (a, b, c) = self.get_abc_vectors(domain);

        (
            domain.interpolate_univariate(&a),
            domain.interpolate_univariate(&b),
            domain.interpolate_univariate(&c),
        )
    }
    
    pub fn get_solution(&self, domain: &MultiplicativeSubgroup<F>, k1: F, k2: F) -> Solution<F> {
        let (a, b, c) = self.get_abc_polys(domain);
        let (ql, qr, qm , qo, qc) = self.get_selectors(domain);
        let (sid_1, sid_2, sid_3) = self.get_s_id_polys(domain, k1, k2);
        let (s_sigma_1, s_sigma_2, s_sigma_3) = self.get_sigma_polys(domain, k1, k2);
        let pi = self.get_public_input_poly(domain);
        
        Solution {
            a,
            b,
            c,
            ql,
            qr,
            qm,
            qo,
            qc,
            sid_1,
            sid_2,
            sid_3,
            s_sigma_1,
            s_sigma_2,
            s_sigma_3,
            pi,
        }
    }
}

// struct Gate<F: FftField + PrimeField> {
//     left: usize,
//     right: usize,
//     output: usize,
//     // out_result: F,
//     constant: F,
// }

#[derive(Clone, Debug)]
struct ArithmeticGate {
        left: usize,
        right: usize,
        output: usize,
        // out_result: F,
        // constant: F,
}

#[derive(Clone, Debug)]
enum Gate {
    Addition (ArithmeticGate),
    Multiplication (ArithmeticGate),
    // Constant {
    //     var_i: usize,
    //     e: F,
    // }
}

impl Gate {
    pub fn execute<F: FftField + PrimeField>(&self, left: F, right: F) -> F {
        match self {
            Gate::Addition(_) => left + right,
            Gate::Multiplication(_) => left * right,
            // _ => panic!()
        }
    }
}

// struct CircuitBuilder<F: FftField + PrimeField> {
//     // gates: Vec<CircuitGate<F>>,
//     // variables: Vec<F>,
//     variables_count: usize,
//     witness: Vec<usize>,
//     public_inputs: Vec<usize>,
//     gates: Vec<Gate<F>>,
// }

#[derive(Clone, Debug)]
struct CircuitDescription<F: FftField + PrimeField> {
    variables_count: usize,
    witness: Vec<usize>,
    public_inputs: Vec<usize>,
    gates: Vec<Gate>,
    constants: Vec<(usize, F)>
}

enum VarOrConst<F> {
    Var(usize),
    Const(F)
}

#[derive(Debug, Clone)]
struct CompiledGate2<F: PrimeField + FftField> {
    // left: VarOrConst<F>,
    left: usize,
    right: usize,
    out: usize,
    ql: bool,
    qr: bool,
    qm: bool,
    qo: bool,
    constant: F,
    operation: Option<Operation>
}

#[derive(Debug, Clone)]
enum Operation {
    Multiplication,
    Addition,
}

struct CircuitCompiler<'a, F: FftField + PrimeField> {
    public_input_count: usize,
    public_inputs: Vec<usize>,
    constants: Vec<(usize, F)>,
    sigma: Vec<usize>,
    gates: Vec<CompiledGate2<F>>,
    domain: &'a MultiplicativeSubgroup<F>,
    k1: F,
    k2: F,
}

impl<'a, F: FftField + PrimeField> CircuitCompiler<'a, F> {
    fn new(
        public_input_count: usize,
        sigma: Vec<usize>,
        gates: Vec<CompiledGate2<F>>,
        domain: &'a MultiplicativeSubgroup<F>,
        k1: F,
        k2: F,
        public_inputs: Vec<usize>,
        constants: Vec<(usize, F)>,
    ) -> Self {
        Self { public_input_count, gates, domain, k1, k2, sigma, public_inputs, constants }
    }

    pub fn compile(self) -> CompiledCircuit2<'a, F>  {
        let (ql, qr, qm , qo, qc) = self.get_selectors();
        let (sid_1, sid_2, sid_3) = self.get_s_id_polys();
        let (s_sigma_1, s_sigma_2, s_sigma_3) = self.get_sigma_polys();

        CompiledCircuit2 {
            // circuit_description,
            domain: self.domain,
            public_inputs: self.public_inputs,
            constants: self.constants,
            // constants: circuit_description.constants.clone(),
            ql,
            qr,
            qm,
            qo,
            qc,
            sid_1,
            sid_2,
            sid_3,
            s_sigma_1,
            s_sigma_2,
            s_sigma_3,
            gates: self.gates,
            public_inputs_count: self.public_input_count,
            sigma: self.sigma,
        }
    }

    pub fn get_selectors(&self) -> (
        DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>
    ) {
        let mut ql_values = vec![];
        let mut qr_values = vec![];
        let mut qm_values = vec![];
        let mut qo_values = vec![];
        let mut qc_values = vec![];

        for gate in &self.gates {
            ql_values.push(format_bool(gate.ql));
            qr_values.push(format_bool(gate.qr));
            qm_values.push(format_bool(gate.qm));
            qo_values.push(-format_bool::<F>(gate.qo));
            qc_values.push(gate.constant);
        }

        // let padded_n = self.padded_len();
        let padded_n = self.domain.len();
        for _ in 0..padded_n - self.gates.len() {
            ql_values.push(format_bool(false));
            qr_values.push(format_bool(false));
            qm_values.push(format_bool(false));
            qo_values.push(format_bool(false));
            qc_values.push(F::zero());
        }

        let ql = self.domain.interpolate_univariate(&ql_values);
        let qr = self.domain.interpolate_univariate(&qr_values);
        let qm = self.domain.interpolate_univariate(&qm_values);
        let qo = self.domain.interpolate_univariate(&qo_values);
        let qc = self.domain.interpolate_univariate(&qc_values);

        (ql, qr, qm , qo, qc)
    }

    fn sigma(&self, from: usize) -> usize {
        self.sigma[from]
    }

    pub fn get_sigma_polys(&self) -> (DP<F>, DP<F>, DP<F>) {
        (
            self.s_sigma_poly(1),
            self.s_sigma_poly(2),
            self.s_sigma_poly(3),
        )
    }

    pub fn get_s_id_polys(&self) -> (DP<F>, DP<F>, DP<F>) {
        (
            self.s_id_poly(1),
            self.s_id_poly(2),
            self.s_id_poly(3),
        )
    }

    fn s_sigma_poly(&self, col: usize) -> DensePolynomial<F> {
        let mut values = vec![];
        let n = self.gates.len();

        for i in 0..n {
            let permuted_index = self.sigma((col - 1) * n + i);
            let permuted_omega = self.map_index_to_coset_value(permuted_index);
            values.push(permuted_omega);
        }

        // let padded_n = self.padded_len();
        let padded_n = self.domain.len();
        for i in self.gates.len()..padded_n {
            values.push(self.get_coset_shifter(col) * self.domain[i]);
        }

        self.domain.interpolate_univariate(&values)
    }

    fn s_id_poly(&self, col: usize) -> DP<F> {
        let mut values = vec![];
        // let n = self.padded_len();
        let n = self.domain.len();

        for i in 0..n {
            let permuted_omega = self.get_coset_shifter(col) * self.domain[i];
            values.push(permuted_omega);
        }

        self.domain.interpolate_univariate(&values)
    }

    fn map_index_to_coset_value(&self, index: usize) -> F {
        let n = self.gates.len();

        match index {
            index if index < n => self.domain[index],
            index if index >= n && index < 2 * n => self.k1 * self.domain[index - n],
            _ => self.k2 * self.domain[index - 2 * n]
        }
    }

    fn get_coset_shifter(&self, col: usize) -> F {
        match col {
            1 => F::one(),
            2 => self.k1,
            3 => self.k2,
            _ => panic!("invalid col {}", col),
        }
    }
}

impl<F: FftField + PrimeField> CircuitDescription<F> {
    pub fn compile<'a>(&self, domain: &'a MultiplicativeSubgroup<F>) -> CompiledCircuit2<'a, F> {
        let mut compiled_gates = vec![];
        let mut variables_map = HashMap::new();

        for pi_var_i in &self.public_inputs {
            compiled_gates.push(CompiledGate2 {
                // left: VarOrConst::Var(*pi_var_i),
                left: *pi_var_i,
                right: 0,
                out: 0,
                ql: true,
                qr: false,
                qm: false,
                qo: false,
                constant: F::zero(),
                operation: None,
            });
        }

        for (const_i, e) in &self.constants {
            variables_map.insert(*const_i, *e);
            compiled_gates.push(CompiledGate2 {
                // left: VarOrConst::Const(*e),
                left: *const_i,
                right: 0,
                out: 0,
                ql: true,
                qr: false,
                qm: false,
                qo: false,
                constant: -*e,
                operation: None,
            });
        }

        for gate in &self.gates {
            match gate {
                Gate::Addition (gate) => {
                    compiled_gates.push(CompiledGate2 {
                        // left: VarOrConst::Var(gate.left),
                        left: gate.left,
                        right: gate.right,
                        out: gate.output,
                        ql: true,
                        qr: true,
                        qo: true,
                        qm: false,
                        constant: F::zero(),
                        operation: Some(Operation::Addition),
                    });
                },
                Gate::Multiplication (gate) => {
                    compiled_gates.push(CompiledGate2 {
                        // left: VarOrConst::Var(gate.left),
                        left: gate.left,
                        right: gate.right,
                        out: gate.output,
                        ql: false,
                        qr: false,
                        qo: true,
                        qm: true,
                        constant: F::zero(),
                        operation: Some(Operation::Multiplication),
                    });
                }
            }
        }

        let (k1, k2) = pick_coset_shifters(domain);

        CircuitCompiler::new(
            self.public_inputs.len(),
            Self::build_sigma_permutations(&compiled_gates),
            compiled_gates,
            domain,
            k1,
            k2,
            self.public_inputs.clone(),
            self.constants.clone(),
        ).compile()

        // CompiledCircuit2 {
        //     public_inputs_count: self.public_inputs.len(),
        //     ql: DensePolynomial::zero(),
        //     qc: DensePolynomial::zero(),
        //     qm: DensePolynomial::zero(),
        //     qr: DensePolynomial::zero(),
        //     qo: DensePolynomial::zero(),
        //     sid_1: DensePolynomial::zero(),
        //     sid_2: DensePolynomial::zero(),
        //     sid_3: DensePolynomial::zero(),
        //     s_sigma_1: DensePolynomial::zero(),
        //     s_sigma_3: DensePolynomial::zero(),
        //     s_sigma_2: DensePolynomial::zero(),
        //     sigma: Self::build_sigma_permutations(&compiled_gates),
        //     gates: compiled_gates,
        // }
    }

    fn build_sigma_permutations(gates: &[CompiledGate2<F>]) -> Vec<usize> {
        // let mut left_i = vec![];
        // let mut right_i = vec![];
        // let mut out_i = vec![];

        let mut items = vec![];

        for gate in gates {
            items.push((gate.left, 0));
        }

        for gate in gates {
            items.push((gate.right, 1));
        }

        for gate in gates {
            items.push((gate.out, 2));
        }

        // let items: Vec<(usize, usize)> = (0..=2)
        //     .flat_map(|col| {
        //         (0..gates.len()).map(move |index| (index, col))
        //     })
        //     .collect();

        println!("items {:?}", items);

        let mut prev_index_col_map = HashMap::new();
        let mut first_index_map = HashMap::new();
        let mut sigma = vec![0; gates.len() * 3];

        for (i, (var_id, col)) in items.into_iter().enumerate() {
            let current = prev_index_col_map.remove(&var_id);

            match current {
                None => {
                    prev_index_col_map.insert(var_id, i);
                    first_index_map.insert(var_id, i);
                    sigma[i] = i;
                },
                Some(last_i) => {
                    prev_index_col_map.insert(var_id, i);
                    sigma[i] = *first_index_map.get(&var_id).unwrap();
                    sigma[last_i] = i;
                }
            };
        }

        sigma
    }
    pub fn solve(&self, public_inputs: &[F], witness: &[F]) -> CompiledCircuit<F> {
        assert_eq!(public_inputs.len(), self.public_inputs.len());
        assert_eq!(witness.len(), self.witness.len());

        let mut variables_map = HashMap::<usize, F>::new();
        let mut compiled_gates = vec![];

        for (public_input_index, value) in self.public_inputs.iter().zip(public_inputs) {
            variables_map.insert(*public_input_index, *value);

            // compiled_gates.push(CompiledGate {
            //     left: *value,
            // });
        }

        for (witness_index, value) in self.witness.iter().zip(witness) {
            let insert_res = variables_map.insert(*witness_index, *value);

            if insert_res.is_some() {
                panic!("index already present");
            }
        }

        for gate in &self.gates {
            match gate {
                Gate::Multiplication(gate) => {
                    println!("gate {:?}", gate);
                    let left = *variables_map.get(&gate.left).unwrap();
                    let right = *variables_map.get(&gate.right).unwrap();

                    let res = left * right;

                    variables_map.insert(gate.output, res);
                    compiled_gates.push(CompiledGate {
                        left,
                        right,
                        output: res,
                        ql: false,
                        qr: false,
                        qm: true,
                        qo: true,
                        qc: F::zero(),
                    });
                },
                Gate::Addition(gate) => {
                    let left = *variables_map.get(&gate.left).unwrap();
                    let right = *variables_map.get(&gate.right).unwrap();

                    let res = left + right;

                    variables_map.insert(gate.output, res);
                    compiled_gates.push(CompiledGate {
                        left,
                        right,
                        output: res,
                        ql: true,
                        qr: true,
                        qm: false,
                        qo: true,
                        qc: F::zero(),
                    });
                },
                _ => panic!(),
            }
        }

        panic!()
    }
}

struct CircuitBuilder<F: FftField + PrimeField> {
    circuit_description: CircuitDescription<F>
}

impl<F: FftField + PrimeField> CircuitBuilder<F> {
    pub fn new() -> Self {
        CircuitBuilder {
            // variables: vec![F::zero()],
            // variables_count: 1, // reserved for zero
            // witness: vec![],
            // public_inputs: vec![],
            // gates: vec![],
            circuit_description: CircuitDescription {
                variables_count: 1, // reserved for zero
                witness: vec![],
                public_inputs: vec![],
                gates: vec![],
                constants: vec![],
            }
        }
    }

    fn get_zero_var_index(&self) -> usize {
        0
    }

    pub fn add_witness(&mut self) -> usize {
        let index = self.add_variable();
        self.circuit_description.witness.push(index);
        index
    }

    pub fn add_public_input(&mut self) -> usize {
        let index = self.add_variable();
        self.circuit_description.public_inputs.push(index);
        index
    }

    pub fn make_public(&mut self, var_i: usize) {
        self.circuit_description.public_inputs.push(var_i);
    }

    pub fn multiplication_gate(&mut self, left: usize, right: usize) -> usize {
        // let left_var = self.variables[left];
        // let right_var = self.variables[right];

        // let out_result = left_var * right_var;
        let output = self.add_variable();

        // self.gates.push(Gate {
        //     left,
        //     right,
        //     output,
        //     // out_result,
        //     constant: F::zero(),
        // });
        self.circuit_description.gates.push(Gate::Multiplication(ArithmeticGate {
            left,
            right,
            output,
        }));

        output
    }

    pub fn addition_gate(&mut self, left: usize, right: usize) -> usize {
        // let left_var = self.variables[left];
        // let right_var = self.variables[right];
        //
        // let out_result = left_var + right_var;
        let output = self.add_variable();

        // self.gates.push(Gate {
        //     left,
        //     right,
        //     output,
        //     // out_result,
        //     constant: F::zero(),
        // });
        self.circuit_description.gates.push(Gate::Addition(ArithmeticGate {
            left,
            right,
            output,
        }));

        output
    }

    pub fn constant_var(&mut self, e: F) -> usize {
        let const_i = self.add_variable();

        self.circuit_description.constants.push((const_i, e));

        // self.gates.push(Gate {
        //     left: const_i,
        //     right: self.get_zero_var_index(),
        //     output: self.get_zero_var_index(),
        //     out_result: F::zero(),
        //     constant: -e,
        // });
        // self.circuit_description.gates.push(Gate::Constant{ var_i: const_i, e });

        const_i
    }

    fn add_variable(&mut self) -> usize {
        let new_var = self.circuit_description.variables_count;

        self.circuit_description.variables_count += 1;

        new_var
    }

    fn make(self) -> CircuitDescription<F> {
        self.circuit_description
    }
}

// struct CircuitDescription<F: FftField + PrimeField> {
//     public_inputs: Vec<usize>,
//     constants: Vec<(usize, F)>,
//     gates: Vec<Gate<F>>,
// }

fn build_test_circuit<F: FftField + PrimeField>() -> CircuitDescription<F> {
    let mut circuit_builder = CircuitBuilder::new();

    // let a = circuit_builder.constant_var(F::from(9));
    let a = circuit_builder.add_variable();
    let b = circuit_builder.constant_var(F::from(82));
    // let b = circuit_builder.add_variable();
    circuit_builder.make_public(a);

    let mul_result_1 = circuit_builder.multiplication_gate(a, b);
    let mul_result_2 = circuit_builder.multiplication_gate(mul_result_1, mul_result_1);
    let add_result = circuit_builder.addition_gate(mul_result_2, b);

    circuit_builder.make_public(add_result);

    circuit_builder.make()
}

struct CompiledCircuit2<'a, F: FftField + PrimeField> {
    // circuit_description: CircuitDescription<F>,
    domain: &'a MultiplicativeSubgroup<F>,
    public_inputs_count: usize,
    public_inputs: Vec<usize>,
    constants: Vec<(usize, F)>,
    pub ql: DensePolynomial<F>,
    pub qr: DensePolynomial<F>,
    pub qm: DensePolynomial<F>,
    pub qo: DensePolynomial<F>,
    pub qc: DensePolynomial<F>,
    pub sid_1: DensePolynomial<F>,
    pub sid_2: DensePolynomial<F>,
    pub sid_3: DensePolynomial<F>,
    pub s_sigma_1: DensePolynomial<F>,
    pub s_sigma_2: DensePolynomial<F>,
    pub s_sigma_3: DensePolynomial<F>,
    gates: Vec<CompiledGate2<F>>,
    sigma: Vec<usize>,
}

struct Solution2<'a, F: FftField + PrimeField> {
    domain: &'a MultiplicativeSubgroup<F>,
    solution_gates: Vec<GateSolution<F>>
}

impl<F: FftField + PrimeField> Solution2<'_, F> {
    pub fn get_abc_polys(&self) -> (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>) {
        let mut a = vec![];
        let mut b = vec![];
        let mut c = vec![];

        for solution_gate in &self.solution_gates {
            a.push(solution_gate.left);
            b.push(solution_gate.right);
            c.push(solution_gate.out);
        }

        let a = self.domain.interpolate_univariate(&a);
        let b = self.domain.interpolate_univariate(&b);
        let c = self.domain.interpolate_univariate(&c);

        (a, b, c)
    }
}

#[derive(Debug)]
struct GateSolution<F: FftField + PrimeField> {
    left: F,
    right: F,
    out: F,
}

impl<'a, F: FftField + PrimeField> CompiledCircuit2<'a, F> {
    pub fn solve(&self, public_input: &[F], witness: &[F]) -> Solution2<F> {
        assert_eq!(public_input.len(), self.public_inputs_count);
        // assert_eq!()

        let mut solution_gates = vec![];
        let mut variables_map = HashMap::new();
        variables_map.insert(0usize, F::zero());

        for (pi_index, value) in self.public_inputs.iter().zip(public_input) {
            variables_map.insert(*pi_index, *value);
        }

        for (const_i, value) in &self.constants {
            variables_map.insert(*const_i, *value);
        }

        println!("vars map {:?}", variables_map);

        for gate in &self.gates {
            println!("gate {:?}", gate);
            let left = *variables_map.get(&gate.left).unwrap();
            let right = *variables_map.get(&gate.right).unwrap();

            let out = match gate.operation {
                Some(Operation::Multiplication) => {
                    let res = left * right;
                    variables_map.insert(gate.out, res);

                    res
                },
                Some(Operation::Addition) => {
                    let res = left + right;
                    variables_map.insert(gate.out, res);

                    res
                },
                None => *variables_map.get(&gate.out).unwrap(),
            };

            solution_gates.push(GateSolution {
                left,
                right,
                out,
            });
        }

        Solution2 {
            domain: self.domain,
            solution_gates
        }
    }
}


pub fn get_test_circuit<F: FftField + PrimeField>() -> CompiledCircuit<F> {
    let circuit = CompiledCircuit {
        public_inputs_count: 2,
        sigma: vec![
            3, // 0 -> 3,
            17, // 1 -> 17
            9, // 2 -> 9
            0, // 3 -> 0
            10, // 4 -> 10,
            16, // 5 -> 16,
            7, // 6 -> 7 // 0
            8, // 7 -> 8 // 0
            12, // 8 -> 12 // 0
            11, // 9 -> 11
            15, // 10 -> 15
            2, // 11 -> 2
            13, // 12 -> 13 // 0
            14, // 13 -> 14 // 0
            6, // 14 -> 6 // 0
            4, // 15 -> 4
            5, // 16 -> 5
            1, // 17 -> 1
        ],
        gates: vec![
            CompiledGate { // PI1
                left: F::from(9),
                right: F::from(0),
                output: F::from(0),
                ql: true,
                qr: false,
                qm: false,
                qo: false,
                // qc: F::from(-9),
                qc: F::from(0),
            },
            CompiledGate { // solution
                left: F::from(544726),
                right: F::from(0),
                output: F::from(0),
                ql: true,
                qr: false,
                qm: false,
                qo: false,
                // qc: F::from(-544726),
                qc: F::from(0),
            },
            CompiledGate { // const
                left: F::from(82),
                right: F::from(0),
                output: F::from(0),
                ql: true,
                qr: false,
                qm: false,
                qo: false,
                qc: F::from(-82),
            },
            CompiledGate { // mult gate
                left: F::from(9),
                right: F::from(82),
                output: F::from(738),
                ql: false,
                qr: false,
                qm: true,
                qo: true,
                qc: F::from(0),
            },
            CompiledGate { // mult gate
                left: F::from(738),
                right: F::from(738),
                output: F::from(544644),
                ql: false,
                qr: false,
                qm: true,
                qo: true,
                qc: F::from(0),
            },
            CompiledGate { // add gate
                left: F::from(544644),
                right: F::from(82),
                output: F::from(544726),
                ql: true,
                qr: true,
                qm: false,
                qo: true,
                qc: F::from(0),
            },
        ],
    };

    circuit
}


#[cfg(test)]
mod tests {
    use ark_poly::Polynomial;
    use ark_std::iterable::Iterable;
    use ark_std::Zero;
    use ark_test_curves::bls12_381::Fr;
    use crate::plonk::circuit::{build_test_circuit, get_test_circuit, ArithmeticGate, CircuitGate, CompiledCircuit, CompiledGate, Gate, Input, Variable};
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::plonk::prover::pick_coset_shifters;

    #[test]
    pub fn test_circuit_description() {
        let test_circuit = build_test_circuit::<Fr>();

        assert_eq!(test_circuit.variables_count, 6);
        assert_eq!(test_circuit.public_inputs, vec![1, 5]);
        assert_eq!(test_circuit.constants, vec![(2, Fr::from(82))]);
        assert_eq!(test_circuit.gates.len(), 3);
        assert!(matches!(test_circuit.gates[0], Gate::Multiplication(ArithmeticGate { left: 1, right: 2, output: 3 })));
        assert!(matches!(test_circuit.gates[1], Gate::Multiplication(ArithmeticGate { left: 3, right: 3, output: 4 })));
        assert!(matches!(test_circuit.gates[2], Gate::Addition(ArithmeticGate { left: 4, right: 2, output: 5 })));
    }

    #[test]
    pub fn test_circuit_compile() {
        let test_circuit = build_test_circuit::<Fr>();
        let domain = generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>();
        let compiled_circuit = test_circuit.compile(&domain);

        assert_eq!(compiled_circuit.sigma, vec![
            3, // 0 -> 3,
            17, // 1 -> 17
            9, // 2 -> 9
            0, // 3 -> 0
            10, // 4 -> 10,
            16, // 5 -> 16,
            7, // 6 -> 6 // 0
            8, // 7 -> 7 // 0
            12, // 8 -> 8 // 0
            11, // 9 -> 11
            15, // 10 -> 15
            2, // 11 -> 2
            13, // 12 -> 12 // 0
            14, // 13 -> 13 // 0
            6, // 14 -> 14 // 0
            4, // 15 -> 4
            5, // 16 -> 5
            1, // 17 -> 1
        ]);

        assert_eq!(compiled_circuit.public_inputs_count, test_circuit.public_inputs.len());

        assert_eq!(compiled_circuit.gates[0].ql, true);
        assert_eq!(compiled_circuit.gates[0].qr, false);
        assert_eq!(compiled_circuit.gates[0].qm, false);
        assert_eq!(compiled_circuit.gates[0].qo, false);
        assert_eq!(compiled_circuit.gates[0].constant, Fr::zero());
        assert_eq!(compiled_circuit.gates[0].left, 1);
        assert_eq!(compiled_circuit.gates[0].right, 0);
        assert_eq!(compiled_circuit.gates[0].out, 0);

        assert_eq!(compiled_circuit.gates[1].ql, true);
        assert_eq!(compiled_circuit.gates[1].qr, false);
        assert_eq!(compiled_circuit.gates[1].qm, false);
        assert_eq!(compiled_circuit.gates[1].qo, false);
        assert_eq!(compiled_circuit.gates[1].constant, Fr::zero());
        assert_eq!(compiled_circuit.gates[1].left, 5);
        assert_eq!(compiled_circuit.gates[1].right, 0);
        assert_eq!(compiled_circuit.gates[1].out, 0);

        assert_eq!(compiled_circuit.gates[2].ql, true);
        assert_eq!(compiled_circuit.gates[2].qr, false);
        assert_eq!(compiled_circuit.gates[2].qm, false);
        assert_eq!(compiled_circuit.gates[2].qo, false);
        assert_eq!(compiled_circuit.gates[2].constant, -Fr::from(82));
        assert_eq!(compiled_circuit.gates[2].left, 2);
        assert_eq!(compiled_circuit.gates[2].right, 0);
        assert_eq!(compiled_circuit.gates[2].out, 0);

        assert_eq!(compiled_circuit.gates[3].ql, false);
        assert_eq!(compiled_circuit.gates[3].qr, false);
        assert_eq!(compiled_circuit.gates[3].qm, true);
        assert_eq!(compiled_circuit.gates[3].qo, true);
        assert_eq!(compiled_circuit.gates[3].constant, Fr::zero());
        assert_eq!(compiled_circuit.gates[3].left, 1);
        assert_eq!(compiled_circuit.gates[3].right, 2);
        assert_eq!(compiled_circuit.gates[3].out, 3);

        assert_eq!(compiled_circuit.gates[4].ql, false);
        assert_eq!(compiled_circuit.gates[4].qr, false);
        assert_eq!(compiled_circuit.gates[4].qm, true);
        assert_eq!(compiled_circuit.gates[4].qo, true);
        assert_eq!(compiled_circuit.gates[4].constant, Fr::zero());
        assert_eq!(compiled_circuit.gates[4].left, 3);
        assert_eq!(compiled_circuit.gates[4].right, 3);
        assert_eq!(compiled_circuit.gates[4].out, 4);

        assert_eq!(compiled_circuit.gates[5].ql, true);
        assert_eq!(compiled_circuit.gates[5].qr, true);
        assert_eq!(compiled_circuit.gates[5].qm, false);
        assert_eq!(compiled_circuit.gates[5].qo, true);
        assert_eq!(compiled_circuit.gates[5].constant, Fr::zero());
        assert_eq!(compiled_circuit.gates[5].left, 4);
        assert_eq!(compiled_circuit.gates[5].right, 2);
        assert_eq!(compiled_circuit.gates[5].out, 5);
    }

    #[test]
    pub fn test_circuit_solve() {
        let test_circuit = build_test_circuit::<Fr>();
        let domain = generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>();
        let compiled_circuit = test_circuit.compile(&domain);
        let solution = compiled_circuit.solve(&[Fr::from(9), Fr::from(544726)], &[]);

        let expected = vec![(9, 0, 0), (544726, 0, 0), (82, 0, 0), (9, 82, 738), (738, 738, 544644), (544644, 82, 544726)];

        assert_eq!(expected.len(), solution.solution_gates.len());

        for ((left, right, out), gate) in expected.into_iter().zip(solution.solution_gates) {
            assert_eq!(Fr::from(left), gate.left);
            assert_eq!(Fr::from(right), gate.right);
            assert_eq!(Fr::from(out), gate.out);
        }
    }

    #[test]
    pub fn test_abc() {
        let domain = generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>();
        let (k1, k2) = pick_coset_shifters(&domain);
        let test_circuit = get_test_circuit::<Fr>();

        let (a,b,c) = test_circuit.get_abc_vectors(&domain);
        println!("a {:?}", a);
        // let (a_poly, b_poly, c_poly) = test_circuit.get_

        // for
    }

    #[test]
    pub fn test_s_id() {
        let domain = generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>();
        let (k1, k2) = pick_coset_shifters(&domain);
        let k1_coset = domain.iter().map(|e| k1 * *e).collect::<Vec<_>>();
        let k2_coset = domain.iter().map(|e| k2 * *e).collect::<Vec<_>>();
        let test_circuit = get_test_circuit();

        let (sid_1, sid_2, sid_3) = test_circuit.get_s_id_polys(&domain, k1, k2);

        for e in &domain {
            assert_eq!(*e, sid_1.evaluate(e));
        }

        for (ke, e) in k1_coset.iter().zip(domain.iter()) {
            assert_eq!(*ke, sid_2.evaluate(e));
        }

        for (ke, e) in k2_coset.iter().zip(domain.iter()) {
            assert_eq!(*ke, sid_3.evaluate(e));
        }
    }

    #[test]
    pub fn test_sigma_poly() {
        return;
        let domain = generate_multiplicative_subgroup::<{ 1<< 3 }, Fr>();
        let (k1, k2) = pick_coset_shifters(&domain);
        let k1_coset = domain.iter().map(|e| k1 * *e).collect::<Vec<_>>();
        let k2_coset = domain.iter().map(|e| k2 * *e).collect::<Vec<_>>();
        let test_circuit = get_test_circuit();

        let (sigma_1, sigma_2, sigma_3) = test_circuit.get_sigma_polys(&domain, k1, k2);

        let sigma_1_permutations = [
            domain[3],
            k2_coset[5],
            k1_coset[3],
            domain[0],
            k1_coset[4],
            k2_coset[4],
            domain[6],
            domain[7]
        ];
        let sigma_2_permutations = [
            k1_coset[0],
            k1_coset[1],
            k1_coset[2],
            k1_coset[5],
            k2_coset[3],
            domain[2],
            k1_coset[6],
            k1_coset[7],
        ];
        let sigma_3_permutations = [
            k2_coset[0],
            k2_coset[1],
            k2_coset[2],
            domain[4],
            domain[5],
            domain[1],
            k2_coset[6],
            k2_coset[7],
        ];

        for ((i, x), y) in domain.iter().enumerate().zip(sigma_1_permutations.iter()) {
            assert_eq!(sigma_1.evaluate(x), y);
        }

        for (x, y) in domain.iter().zip(sigma_2_permutations.iter()) {
            assert_eq!(sigma_2.evaluate(x), y);
        }

        for (x, y) in domain.iter().zip(sigma_3_permutations.iter()) {
            assert_eq!(sigma_3.evaluate(x), y);
        }
    }
}