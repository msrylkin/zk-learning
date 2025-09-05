use std::collections::HashMap;
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::univariate::DensePolynomial;
use ark_std::iterable::Iterable;
use ark_test_curves::bls12_381::Fr;
use crate::poly_utils::interpolate_univariate;

struct CircuitBuilder<F: FftField + PrimeField> {
    circuit: Circuit<F>,
    vars: HashMap<usize, Variable>,
    last_var: usize,
}

type DP<F> = DensePolynomial<F>;

enum Input<F: FftField + PrimeField> {
    Public(Variable),
    Witness(Variable),
    Constant(F),
    InnerWitness,
}

type Variable = (usize, usize); // (col, row)

impl<F: PrimeField + FftField> CircuitBuilder<F> {
    fn new() -> Self {
        Self {
            circuit: Circuit {
                gates: Vec::new(),
            },
            vars: HashMap::new(),
            last_var: 0,
        }
    }

    fn create_var() {

    }

    // fn create_public_var(&mut self) -> usize {
    //     // let var_id = self.last_var;
    //     // let col = self.circuit.gates.len() - 1;
    //     // self.circuit.gates.push(CircuitGate {
    //     //     left: Input::Public((col, 0)),
    //     //     right: Input::
    //     // });
    //     // self.last_var += 1;
    //     //
    //     // var_id
    // }

    // fn add_gate(
    //     self,
    //     left: Input<F>,
    //     right: Input<F>,
    //     out: Input<F>,
    // ) -> Self {
    //     let circuit = self.circuit;
    //     let mut gates = circuit.gates;
    //
    //
    // }

    fn mul_gate() {

    }
}

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

fn format_bool<F: Field>(b: bool) -> F {
    match b {
        true => F::one(),
        false => F::zero(),
    }
}

impl<F: FftField + PrimeField> CompiledCircuit<F> {
    pub fn get_abc_vectors(&self) -> (Vec<F>, Vec<F>, Vec<F>) {
        let mut a = vec![];
        let mut b = vec![];
        let mut c = vec![];

        for gate in &self.gates {
            a.push(gate.left);
            b.push(gate.right);
            c.push(gate.output);
        }

        let padded_n = self.padded_len();
        for _ in 0..padded_n - self.gates.len() {
            a.push(F::zero());
            b.push(F::zero());
            c.push(F::zero());
        }

        (a, b, c)
    }

    pub fn padded_len(&self) -> usize {
        let n = self.gates.len();

        match n.is_power_of_two() {
            true => n,
            false => n.next_power_of_two(),
        }
    }

    pub fn get_sigma(&self) -> Vec<usize> {
        self.sigma.clone()
    }

    pub fn get_selectors(&self, domain: &[F]) -> (
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

        let padded_n = self.padded_len();
        for _ in 0..padded_n - self.gates.len() {
            ql_values.push(format_bool(false));
            qr_values.push(format_bool(false));
            qm_values.push(format_bool(false));
            qo_values.push(format_bool(false));
            qc_values.push(F::zero());
        }

        let ql = interpolate_univariate(domain, &ql_values);
        let qr = interpolate_univariate(domain, &qr_values);
        let qm = interpolate_univariate(domain, &qm_values);
        let qo = interpolate_univariate(domain, &qo_values);
        let qc = interpolate_univariate(domain, &qc_values);

        (ql, qr, qm , qo, qc)
    }

    fn sigma(&self, from: usize) -> usize {
        self.sigma[from]
    }

    pub fn get_sigma_polys(&self, domain: &[F], k1: F, k2: F) -> (DP<F>, DP<F>, DP<F>) {
        (
            self.s_sigma_poly(1, domain, k1, k2),
            self.s_sigma_poly(2, domain, k1, k2),
            self.s_sigma_poly(3, domain, k1, k2),
        )
    }

    pub fn get_s_id_polys(&self, domain: &[F], k1: F, k2: F) -> (DP<F>, DP<F>, DP<F>) {
        (
            self.s_id_poly(1, domain, k1, k2),
            self.s_id_poly(2, domain, k1, k2),
            self.s_id_poly(3, domain, k1, k2),
        )
    }

    fn s_sigma_poly(&self, col: usize, domain: &[F], k1: F, k2: F) -> DensePolynomial<F> {
        let mut values = vec![];
        let n = self.gates.len();

        for i in 0..n {
            let permuted_index = self.sigma((col - 1) * n + i);
            let permuted_omega = self.map_index_to_coset_value(permuted_index, domain, k1, k2);
            values.push(permuted_omega);
        }

        let padded_n = self.padded_len();
        for i in self.gates.len()..padded_n {
            values.push(self.get_coset_shifter(col, k1, k2) * domain[i]);
        }

        interpolate_univariate(domain, &values)
    }

    fn s_id_poly(&self, col: usize, domain: &[F], k1: F, k2: F) -> DP<F> {
        let mut values = vec![];
        let n = self.padded_len();

        for i in 0..n {
            let permuted_omega = self.get_coset_shifter(col, k1, k2) * domain[i];
            values.push(permuted_omega);
        }

        interpolate_univariate(domain, &values)
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

    pub fn get_public_input_poly(&self, domain: &[F]) -> DensePolynomial<F> {
        let values = self
            .get_abc_vectors().0
            .iter()
            .enumerate()
            .map(|(i, e)| {
                if i < self.public_inputs_count {
                    return -*e;
                }

                F::zero()
            })
            .collect::<Vec<_>>();

        interpolate_univariate(&domain, &values)
    }
    
    pub fn get_solution(&self, domain: &[F], k1: F, k2: F) -> Solution<F> {
        let (a, b, c) = self.get_abc_vectors();
        let (a, b, c) = (
            interpolate_univariate(domain, &a),
            interpolate_univariate(domain, &b),
            interpolate_univariate(domain, &c),
        );
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

struct Circuit<F: FftField + PrimeField> {
    gates: Vec<CircuitGate<F>>
}

impl<F: FftField + PrimeField> Circuit<F> {
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
            6, // 6 -> 6
            7, // 7 -> 7
            8, // 8 -> 8
            11, // 9 -> 11
            15, // 10 -> 15
            2, // 11 -> 2
            12, // 12 -> 12
            13, // 13 -> 13
            14, // 14 -> 14
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
    use ark_test_curves::bls12_381::Fr;
    use crate::plonk::circuit::{get_test_circuit, Circuit, CircuitBuilder, CircuitGate, CompiledCircuit, CompiledGate, Input, Variable};
    use crate::plonk::prover::pick_coset_shifters;
    use crate::poly_utils::generate_multiplicative_subgroup;

    #[test]
    pub fn test_s_id() {
        let domain = generate_multiplicative_subgroup::<{ 1<< 3 }, Fr>();
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