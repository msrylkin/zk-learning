use std::collections::HashMap;
use ark_ff::{FftField, Field, PrimeField};
use ark_poly::univariate::DensePolynomial;
use crate::plonk::circuit::circuit_description::{CircuitDescription, Gate};
use crate::plonk::domain::PlonkDomain;
use crate::poly_utils::format_bool;

fn pad_up_to_len<F: Field>(mut vec: Vec<F>, len: usize) -> Vec<F> {
    vec.resize(len, F::zero());

    vec
}

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

pub struct CompiledCircuit<F: FftField + PrimeField> {
    public_inputs: Vec<usize>,
    outputs: Vec<usize>,
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
    gates: Vec<CompiledGate<F>>,
    sigma: Vec<usize>,
}

impl<'a, F: FftField + PrimeField> CompiledCircuit<F> {}

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
                Gate::Addition (gate, is_output) => {
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
                Gate::Multiplication (gate, is_output) => {
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
            items.push((gate.left, 0));
        }

        for gate in gates {
            items.push((gate.right, 1));
        }

        for gate in gates {
            items.push((gate.out, 2));
        }

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
    use crate::plonk::test_utils::{build_test_circuit, build_test_circuit_private_witness};
    use crate::poly_utils::format_bool;

    #[test]
    pub fn test_circuit_compile() {
        let test_circuit = build_test_circuit::<Fr>();
        let domain = generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);
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

        assert_eq!(compiled_circuit.public_inputs.len(), test_circuit.public_inputs_count());

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

        for (gate, w) in compiled_circuit.gates.iter().zip(&domain) {
            assert_eq!(compiled_circuit.ql.evaluate(w), format_bool(gate.ql));
            assert_eq!(compiled_circuit.qr.evaluate(w), format_bool(gate.qr));
            assert_eq!(compiled_circuit.qm.evaluate(w), format_bool(gate.qm));
            assert_eq!(compiled_circuit.qo.evaluate(w), -format_bool::<Fr>(gate.qo));
            assert_eq!(compiled_circuit.qc.evaluate(w), gate.constant);
        }
    }

    #[test]
    pub fn test_circuit_compile_private_input() {
        let test_circuit = build_test_circuit_private_witness::<Fr>();
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
            1, // 26 -> 7
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

        // assert_eq!(compiled_circuit.gates[0].ql, true);
        // assert_eq!(compiled_circuit.gates[0].qr, false);
        // assert_eq!(compiled_circuit.gates[0].qm, false);
        // assert_eq!(compiled_circuit.gates[0].qo, false);
        // assert_eq!(compiled_circuit.gates[0].constant, Fr::zero());
        // assert_eq!(compiled_circuit.gates[0].left, 1);
        // assert_eq!(compiled_circuit.gates[0].right, 0);
        // assert_eq!(compiled_circuit.gates[0].out, 0);
        //
        // assert_eq!(compiled_circuit.gates[1].ql, true);
        // assert_eq!(compiled_circuit.gates[1].qr, false);
        // assert_eq!(compiled_circuit.gates[1].qm, false);
        // assert_eq!(compiled_circuit.gates[1].qo, false);
        // assert_eq!(compiled_circuit.gates[1].constant, Fr::zero());
        // assert_eq!(compiled_circuit.gates[1].left, 5);
        // assert_eq!(compiled_circuit.gates[1].right, 0);
        // assert_eq!(compiled_circuit.gates[1].out, 0);
        //
        // assert_eq!(compiled_circuit.gates[2].ql, true);
        // assert_eq!(compiled_circuit.gates[2].qr, false);
        // assert_eq!(compiled_circuit.gates[2].qm, false);
        // assert_eq!(compiled_circuit.gates[2].qo, false);
        // assert_eq!(compiled_circuit.gates[2].constant, -Fr::from(82));
        // assert_eq!(compiled_circuit.gates[2].left, 2);
        // assert_eq!(compiled_circuit.gates[2].right, 0);
        // assert_eq!(compiled_circuit.gates[2].out, 0);
        //
        // assert_eq!(compiled_circuit.gates[3].ql, false);
        // assert_eq!(compiled_circuit.gates[3].qr, false);
        // assert_eq!(compiled_circuit.gates[3].qm, true);
        // assert_eq!(compiled_circuit.gates[3].qo, true);
        // assert_eq!(compiled_circuit.gates[3].constant, Fr::zero());
        // assert_eq!(compiled_circuit.gates[3].left, 1);
        // assert_eq!(compiled_circuit.gates[3].right, 2);
        // assert_eq!(compiled_circuit.gates[3].out, 3);
        //
        // assert_eq!(compiled_circuit.gates[4].ql, false);
        // assert_eq!(compiled_circuit.gates[4].qr, false);
        // assert_eq!(compiled_circuit.gates[4].qm, true);
        // assert_eq!(compiled_circuit.gates[4].qo, true);
        // assert_eq!(compiled_circuit.gates[4].constant, Fr::zero());
        // assert_eq!(compiled_circuit.gates[4].left, 3);
        // assert_eq!(compiled_circuit.gates[4].right, 3);
        // assert_eq!(compiled_circuit.gates[4].out, 4);
        //
        // assert_eq!(compiled_circuit.gates[5].ql, true);
        // assert_eq!(compiled_circuit.gates[5].qr, true);
        // assert_eq!(compiled_circuit.gates[5].qm, false);
        // assert_eq!(compiled_circuit.gates[5].qo, true);
        // assert_eq!(compiled_circuit.gates[5].constant, Fr::zero());
        // assert_eq!(compiled_circuit.gates[5].left, 4);
        // assert_eq!(compiled_circuit.gates[5].right, 2);
        // assert_eq!(compiled_circuit.gates[5].out, 5);

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
        let domain = generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);
        let compiled_circuit = test_circuit.compile(&domain);
        let pi_vector = vec![Fr::from(9)];
        let solution = test_circuit.solve(&pi_vector, &[], &domain);

        assert_eq!(solution.public_witness.pi_vector, vec![Fr::from(9)]);
        assert_eq!(solution.public_witness.output_vector, vec![Fr::from(544726)]);

        let expected = vec![(9, 0, 0), (544726, 0, 0), (82, 0, 0), (9, 82, 738), (738, 738, 544644), (544644, 82, 544726)];

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

        for i in 0..pi_vector.len() {
            assert_eq!(
                solution.public_witness.pi_combined.evaluate(&domain[i]),
                -pi_vector[i],
            );
            assert_eq!(compiled_circuit.qc.evaluate(&domain[i]), Fr::zero());
        }

        for i in pi_vector.len()..solution.public_witness.output_vector.len() + pi_vector.len() {
            assert_eq!(
                solution.public_witness.pi_combined.evaluate(&domain[i]),
                -solution.public_witness.output_vector[i - pi_vector.len()],
            );
            assert_eq!(compiled_circuit.qc.evaluate(&domain[i]), Fr::zero());
        }

        for i in solution.public_witness.output_vector.len() + pi_vector.len()..compiled_circuit.gates.len() {
            assert_eq!(solution.public_witness.pi_combined.evaluate(&domain[i]), Fr::zero());
            assert_eq!(compiled_circuit.qc.evaluate(&domain[i]), compiled_circuit.gates[i].constant);
        }

        for i in compiled_circuit.gates.len()..domain.len() {
            assert_eq!(solution.public_witness.pi_combined.evaluate(&domain[i]), Fr::zero());
            assert_eq!(compiled_circuit.qc.evaluate(&domain[i]), Fr::zero());
        }
    }
}