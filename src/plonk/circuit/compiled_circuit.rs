use std::collections::HashMap;
use ark_ff::{FftField, PrimeField};
use ark_poly::univariate::DensePolynomial;
use ark_std::iterable::Iterable;
use crate::plonk::circuit::circuit_description::{CircuitDescription, Gate};
use crate::plonk::domain::PlonkDomain;
use crate::poly_utils::format_bool;

type Selectors<F> = (
    DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>
);

#[derive(Debug, Clone)]
struct CompiledGate<F: PrimeField + FftField> {
    left: usize,
    right: usize,
    out: usize,
    ql: bool,
    qr: bool,
    qm: bool,
    qo: bool,
    constant: F,
}

/// Compiled circuit with interpolated public polynomials
pub struct CompiledCircuit<F: FftField + PrimeField> {
    public_inputs: Vec<usize>,
    outputs: Vec<usize>,
    constants: Vec<(usize, F)>,
    /// `ql` poly (left witness poly selector)
    pub ql: DensePolynomial<F>,
    /// `qr` poly (right witness poly selector)
    pub qr: DensePolynomial<F>,
    /// `qm` poly (multiplication selector)
    pub qm: DensePolynomial<F>,
    /// `qo` poly (output column selector)
    pub qo: DensePolynomial<F>,
    /// `qc` poly (constant selector)
    pub qc: DensePolynomial<F>,
    /// `sid_1` poly (identity permutation for left column)
    pub sid_1: DensePolynomial<F>,
    /// `sid_2` poly (identity permutation for right column)
    pub sid_2: DensePolynomial<F>,
    /// `sid_3` poly (identity permutation for output column)
    pub sid_3: DensePolynomial<F>,
    /// `sigma_1` poly (hash permutation poly for left column)
    pub s_sigma_1: DensePolynomial<F>,
    /// `sigma_2` poly (hash permutation poly for right column)
    pub s_sigma_2: DensePolynomial<F>,
    /// `sigma_3` poly (hash permutation poly for output column)
    pub s_sigma_3: DensePolynomial<F>,
    gates: Vec<CompiledGate<F>>,
    sigma: Vec<usize>,
}

pub struct CircuitCompiler<'a, F: FftField + PrimeField> {
    public_inputs: Vec<usize>,
    outputs: Vec<usize>,
    constants: Vec<(usize, F)>,
    gates: Vec<CompiledGate<F>>,
    domain: &'a PlonkDomain<F>,
    sigma: Vec<usize>,
}

impl<'a, 'b, F: FftField + PrimeField> CircuitCompiler<'a, F> {
    pub fn new(
        circuit_description: &'b CircuitDescription<F>,
        domain: &'a PlonkDomain<F>,
    ) -> Self {
        let mut compiled_gates = vec![];
        let mut variables_map = HashMap::new();

        for pi_var_i in [circuit_description.public_inputs(), circuit_description.outputs()].concat() {
            compiled_gates.push(CompiledGate {
                left: pi_var_i,
                right: 0,
                out: 0,
                ql: true,
                qr: false,
                qm: false,
                qo: false,
                constant: F::zero(),
            });
        }

        for (const_i, e) in circuit_description.constants() {
            variables_map.insert(*const_i, *e);
            compiled_gates.push(CompiledGate {
                left: *const_i,
                right: 0,
                out: 0,
                ql: true,
                qr: false,
                qm: false,
                qo: false,
                constant: -*e,
            });
        }

        for gate in circuit_description.gates() {
            match gate {
                Gate::Addition (gate) => {
                    compiled_gates.push(CompiledGate {
                        left: gate.left,
                        right: gate.right,
                        out: gate.output,
                        ql: true,
                        qr: true,
                        qo: true,
                        qm: false,
                        constant: F::zero(),
                    });
                },
                Gate::Multiplication (gate) => {
                    compiled_gates.push(CompiledGate {
                        left: gate.left,
                        right: gate.right,
                        out: gate.output,
                        ql: false,
                        qr: false,
                        qo: true,
                        qm: true,
                        constant: F::zero(),
                    });
                }
            }
        }

        assert!(compiled_gates.len() <= domain.len());

        Self {
            outputs: circuit_description.outputs().to_vec(),
            domain,
            sigma: Self::build_sigma_permutations(&compiled_gates),
            gates: compiled_gates,
            public_inputs: circuit_description.public_inputs().to_vec(),
            constants: circuit_description.constants().to_vec(),
        }
    }

    fn build_sigma_permutations(gates: &[CompiledGate<F>]) -> Vec<usize> {
        let mut items = vec![];

        for gate in gates {
            items.push(gate.left);
        }

        for gate in gates {
            items.push(gate.right);
        }

        for gate in gates {
            items.push(gate.out);
        }

        let mut prev_index_col_map = HashMap::new();
        let mut first_index_map = HashMap::new();
        let mut sigma = vec![0; gates.len() * 3];

        for (i, var_id) in items.into_iter().enumerate() {
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

    pub fn compile(self) -> CompiledCircuit<F>  {
        let (ql, qr, qm , qo, qc) = self.get_selectors();
        let (sid_1, sid_2, sid_3) = self.get_s_id_polys();
        let (s_sigma_1, s_sigma_2, s_sigma_3) = self.get_sigma_polys();

        CompiledCircuit {
            outputs: self.outputs,
            public_inputs: self.public_inputs,
            constants: self.constants,
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
            sigma: self.sigma,
        }
    }

    pub fn get_selectors(&self) -> Selectors<F> {
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

    pub fn get_sigma_polys(&self) -> (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>) {
        (
            self.s_sigma_poly(1),
            self.s_sigma_poly(2),
            self.s_sigma_poly(3),
        )
    }

    pub fn get_s_id_polys(&self) -> (DensePolynomial<F>, DensePolynomial<F>, DensePolynomial<F>) {
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

        let padded_n = self.domain.len();
        for i in self.gates.len()..padded_n {
            values.push(self.get_coset_shifter(col) * self.domain[i]);
        }

        self.domain.interpolate_univariate(&values)
    }

    fn s_id_poly(&self, col: usize) -> DensePolynomial<F> {
        let mut values = vec![];
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
            index if index >= n && index < 2 * n => self.domain.k1() * self.domain[index - n],
            _ => self.domain.k2() * self.domain[index - 2 * n]
        }
    }

    fn get_coset_shifter(&self, col: usize) -> F {
        match col {
            1 => F::one(),
            2 => self.domain.k1(),
            3 => self.domain.k2(),
            _ => panic!("invalid col {}", col),
        }
    }
}

#[cfg(test)]
mod tests {
    use ark_poly::Polynomial;
    use ark_std::Zero;
    use ark_test_curves::bls12_381::Fr;
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::plonk::domain::PlonkDomain;
    use crate::plonk::test_utils::test_utils::build_test_circuit;
    use crate::poly_utils::format_bool;

    #[test]
    pub fn test_circuit_compile() {
        let test_circuit = build_test_circuit::<Fr>();
        let domain = generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);
        let compiled_circuit = test_circuit.compile(&domain);

        assert_eq!(compiled_circuit.sigma, vec![
            5, // 0 -> 5
            7, // 1 -> 7
            29, // 2 -> 29
            15, // 3 -> 15
            18, // 4 -> 18
            0, // 5 -> 0
            16, // 6 -> 16
            9, // 7 -> 9
            27, // 8 -> 27
            26, // 9 -> 26

            11, // 10 -> 11
            12, // 11 -> 12
            13, // 12 -> 13
            14, // 13 -> 14
            20, // 14 -> 20
            3, // 15 -> 3
            25, // 16 -> 25
            17, // 17 -> 17
            4, // 18 -> 4
            28, // 19 -> 28

            21, // 20 -> 21
            22, // 21 -> 22
            23, // 22 -> 23
            24, // 23 -> 24
            10, // 24 -> 10
            6, // 25 -> 6
            1, // 26 -> 1
            8, // 27 -> 8
            19, // 28 -> 19
            2, // 29 -> 2
        ]);

        assert_eq!(compiled_circuit.public_inputs.len(), test_circuit.public_inputs_count());

        // ql, qr, qm, qo, const (qc), left, right, out
        let rows = vec![
            (true, false, false, false, 0, 1, 0, 0),
            (true, false, false, false, 0, 6, 0, 0),
            (true, false, false, false, 0, 9, 0, 0),
            (true, false, false, false, -82, 2, 0, 0),
            (true, false, false, false, 1, 3, 0, 0),
            (false, false, true, true, 0, 1, 2, 5),
            (false, false, true, true, 0, 5, 5, 6),
            (true, true, false, true, 0, 6, 4, 7),
            (false, false, true, true, 0, 7, 3, 8),
            (false, false, true, true, 0, 6, 8, 9),
        ];

        assert_eq!(rows.len(), compiled_circuit.gates.len());

        for (i, (ql, qr, qm, qo, constant, left, right, out)) in rows.into_iter().enumerate() {
            assert_eq!(compiled_circuit.gates[i].ql, ql);
            assert_eq!(compiled_circuit.gates[i].qr, qr);
            assert_eq!(compiled_circuit.gates[i].qm, qm);
            assert_eq!(compiled_circuit.gates[i].qo, qo);
            assert_eq!(compiled_circuit.gates[i].constant, Fr::from(constant));
            assert_eq!(compiled_circuit.gates[i].left, left);
            assert_eq!(compiled_circuit.gates[i].right, right);
            assert_eq!(compiled_circuit.gates[i].out, out);
        }

        for (gate, w) in compiled_circuit.gates.iter().zip(&domain) {
            assert_eq!(compiled_circuit.ql.evaluate(w), format_bool(gate.ql));
            assert_eq!(compiled_circuit.qr.evaluate(w), format_bool(gate.qr));
            assert_eq!(compiled_circuit.qm.evaluate(w), format_bool(gate.qm));
            assert_eq!(compiled_circuit.qo.evaluate(w), -format_bool::<Fr>(gate.qo));
            assert_eq!(compiled_circuit.qc.evaluate(w), gate.constant);
        }
    }

    #[test]
    pub fn test_circuit_solve() {
        let test_circuit = build_test_circuit::<Fr>();
        let domain = generate_multiplicative_subgroup::<{ 1 << 5 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);
        let compiled_circuit = test_circuit.compile(&domain);
        let public_vector = vec![Fr::from(9)];
        let private_vector = vec![Fr::from(-544000)];
        let solution = test_circuit.solve(&public_vector, &private_vector, &domain);

        assert_eq!(solution.public_witness.inputs_vector, vec![Fr::from(9)]);
        assert_eq!(solution.public_witness.output_vector, vec![Fr::from(544644), Fr::from(-350750736)]);

        // a, b, c
        let expected = vec![
            (9, 0, 0),
            (544644, 0, 0),
            (-350750736, 0, 0),
            (82, 0, 0),
            (-1, 0, 0),
            (9, 82, 738),
            (738, 738, 544644),
            (544644, -544000, 644),
            (644, -1, -644),
            (544644, -644, -350750736),
        ];

        assert_eq!(expected.len(), solution.solution_gates.len());

        for (((left, right, out), gate), w) in expected.into_iter().zip(&solution.solution_gates).zip(&domain) {
            assert_eq!(Fr::from(left), gate.left);
            assert_eq!(Fr::from(right), gate.right);
            assert_eq!(Fr::from(out), gate.out);

            assert_eq!(solution.a.evaluate(w), Fr::from(left));
            assert_eq!(solution.b.evaluate(w), Fr::from(right));
            assert_eq!(solution.c.evaluate(w), Fr::from(out));
        }

        for i in solution.solution_gates.len()..domain.len() {
            assert_eq!(solution.a.evaluate(&domain[i]), Fr::zero());
            assert_eq!(solution.b.evaluate(&domain[i]), Fr::zero());
            assert_eq!(solution.c.evaluate(&domain[i]), Fr::zero());
        }

        for i in 0..public_vector.len() {
            assert_eq!(
                solution.public_witness.pi_combined.evaluate(&domain[i]),
                -public_vector[i],
            );
            assert_eq!(compiled_circuit.qc.evaluate(&domain[i]), Fr::zero());
        }

        for i in public_vector.len()..solution.public_witness.output_vector.len() + public_vector.len() {
            assert_eq!(
                solution.public_witness.pi_combined.evaluate(&domain[i]),
                -solution.public_witness.output_vector[i - public_vector.len()],
            );
            assert_eq!(compiled_circuit.qc.evaluate(&domain[i]), Fr::zero());
        }

        for i in solution.public_witness.output_vector.len() + public_vector.len()..compiled_circuit.gates.len() {
            assert_eq!(solution.public_witness.pi_combined.evaluate(&domain[i]), Fr::zero());
            assert_eq!(compiled_circuit.qc.evaluate(&domain[i]), compiled_circuit.gates[i].constant);
        }

        for i in compiled_circuit.gates.len()..domain.len() {
            assert_eq!(solution.public_witness.pi_combined.evaluate(&domain[i]), Fr::zero());
            assert_eq!(compiled_circuit.qc.evaluate(&domain[i]), Fr::zero());
        }
    }
}