use std::collections::HashMap;
use ark_ff::{FftField, Field, PrimeField};
use ark_std::iterable::Iterable;
use crate::plonk::circuit::{CompiledCircuit, GateSolution, PublicWitness};
use crate::plonk::circuit::compiled_circuit::CircuitCompiler;
use crate::plonk::circuit::solution::Solution;
use crate::plonk::domain::PlonkDomain;

#[derive(Clone, Debug)]
pub struct ArithmeticGate {
    pub left: usize,
    pub right: usize,
    pub output: usize,
}

#[derive(Clone, Debug)]
pub enum Gate {
    Addition (ArithmeticGate, bool),
    Multiplication (ArithmeticGate, bool),
}

impl Gate {
    pub fn make_output(&mut self) {
        match self {
            Gate::Addition(_, out) => *out = true,
            Gate::Multiplication(_, out) => *out = true,
        }
    }
}

pub struct GateInfo {
    pub gate_i: usize,
    pub out_var_i: usize,
}

#[derive(Clone, Debug)]
pub struct CircuitDescription<F: FftField + PrimeField> {
    variables_count: usize,
    witness: Vec<usize>,
    public_inputs: Vec<usize>,
    outputs: Vec<usize>,
    gates: Vec<Gate>,
    constants: Vec<(usize, F)>
}

impl<F: FftField + PrimeField> CircuitDescription<F> {
    pub fn new() -> Self {
        Self {
            variables_count: 1, // reserved for zero
            witness: vec![],
            public_inputs: vec![],
            gates: vec![],
            constants: vec![],
            outputs: vec![],
        }
    }

    pub fn public_inputs_count(&self) -> usize {
        self.public_inputs.len()
    }

    pub fn public_inputs(&self) -> &[usize] {
        self.public_inputs.as_slice()
    }

    pub fn outputs(&self) -> &[usize] {
        self.outputs.as_slice()
    }

    pub fn constants(&self) -> &[(usize, F)] {
        self.constants.as_slice()
    }

    pub fn gates(&self) -> &[Gate] {
        self.gates.as_slice()
    }

    pub fn make_public(&mut self, var_i: usize) {
        self.add_public_input(var_i);
    }

    pub fn add_public_input(&mut self, var_id: usize)  {
        self.public_inputs.push(var_id);
    }

    pub fn addition_gate(&mut self, left: usize, right: usize) -> GateInfo {
        let output = self.add_variable();

        self.gates.push(Gate::Addition(ArithmeticGate {
            left,
            right,
            output,
        }, false));

        GateInfo {
            out_var_i: output,
            gate_i: self.gates.len() - 1,
        }
    }

    pub fn multiplication_gate(&mut self, left: usize, right: usize) -> GateInfo {
        let output = self.add_variable();

        self.gates.push(Gate::Multiplication(ArithmeticGate {
            left,
            right,
            output,
        }, false));

        GateInfo {
            out_var_i: output,
            gate_i: self.gates.len() - 1,
        }
    }

    pub fn constant_var(&mut self, e: F) -> usize {
        let const_i = self.add_variable();

        self.constants.push((const_i, e));

        const_i
    }

    pub fn add_variable(&mut self) -> usize {
        let new_var = self.variables_count;

        self.variables_count += 1;

        new_var
    }

    pub fn make_output(&mut self, gate_i: usize) {
        let gate_out_var = self.gate_out_var(gate_i);
        self.outputs.push(gate_out_var);
        self.gates[gate_i].make_output();
    }

    pub fn gate_out_var(&self, gate_i: usize) -> usize {
        match self.gates[gate_i] {
            Gate::Addition(ArithmeticGate {output, ..}, _) => output,
            Gate::Multiplication(ArithmeticGate {output, ..}, _) => output,
        }
    }

    pub fn compile(&self, domain: &PlonkDomain<F>) -> CompiledCircuit<F> {
        CircuitCompiler::new(self, domain).compile()
    }

    pub fn solve(
        &self,
        public_input: &[F],
        private_input: &[F],
        domain: &PlonkDomain<F>,
    ) -> Solution<F> {
        assert_eq!(public_input.len(), self.public_inputs.len());

        let mut variables_map = HashMap::new();
        variables_map.insert(0usize, F::zero());

        let mut pi_values = vec![];
        let mut out_values = vec![];

        let mut constraint_gates = vec![];
        let mut constants_gates = vec![];
        let mut input_gates = vec![];

        for (pi_index, value) in self.public_inputs.iter().zip(public_input) {
            variables_map.insert(*pi_index, *value);
            pi_values.push(*value);
            input_gates.push(GateSolution {
                left: *value,
                right: F::zero(),
                out: F::zero(),
            });
        }

        for (const_i, value) in &self.constants {
            variables_map.insert(*const_i, *value);
            constants_gates.push(GateSolution {
                left: *value,
                right: F::zero(),
                out: F::zero(),
            });
        }

        for gate in &self.gates {
            match gate {
                Gate::Addition(gate, is_out) => {
                    let left = variables_map[&gate.left];
                    let right = variables_map[&gate.right];
                    let res = left + right;

                    constraint_gates.push(GateSolution {
                        left,
                        right,
                        out: res,
                    });

                    variables_map.insert(gate.output, res);
                },
                Gate::Multiplication(gate, is_out) => {
                    let left = variables_map[&gate.left];
                    let right = variables_map[&gate.right];
                    let res = left * right;

                    constraint_gates.push(GateSolution {
                        left,
                        right,
                        out: res,
                    });

                    variables_map.insert(gate.output, res);
                }
            }
        }

        let mut out_gates = vec![];

        for out_i in &self.outputs {
            let out_value = variables_map[out_i];
            out_values.push(out_value);
            out_gates.push(GateSolution {
                left: out_value,
                right: F::zero(),
                out: F::zero(),
            });
        }

        let solution_gates = [input_gates, out_gates, constants_gates, constraint_gates].concat();

        let pi_combined_values = pi_values
            .iter()
            .chain(out_values.iter())
            .map(|e| -*e)
            .collect::<Vec<_>>();
        let pi_combined_values = pad_up_to_len(pi_combined_values, domain.len());

        Solution::new(solution_gates, domain, PublicWitness {
            pi_combined: domain.interpolate_univariate(&pi_combined_values),
            pi_vector: pi_values,
            output_vector: out_values,
        })
    }
}

fn pad_up_to_len<F: Field>(mut vec: Vec<F>, len: usize) -> Vec<F> {
    vec.resize(len, F::zero());

    vec
}

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::Fr;
    use crate::plonk::circuit::circuit_description::{ArithmeticGate, Gate};
    use crate::plonk::test_utils::build_test_circuit;

    #[test]
    pub fn test_circuit_description() {
        let test_circuit = build_test_circuit::<Fr>();

        assert_eq!(test_circuit.variables_count, 6);
        assert_eq!(test_circuit.public_inputs, vec![1]);
        assert_eq!(test_circuit.outputs, vec![5]);
        assert_eq!(test_circuit.constants, vec![(2, Fr::from(82))]);
        assert_eq!(test_circuit.gates.len(), 3);
        assert!(matches!(test_circuit.gates[0], Gate::Multiplication(ArithmeticGate { left: 1, right: 2, output: 3 }, false)));
        assert!(matches!(test_circuit.gates[1], Gate::Multiplication(ArithmeticGate { left: 3, right: 3, output: 4 }, false)));
        assert!(matches!(test_circuit.gates[2], Gate::Addition(ArithmeticGate { left: 4, right: 2, output: 5 }, true)));
    }
}