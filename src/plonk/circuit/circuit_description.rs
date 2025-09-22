use ark_ff::{FftField, PrimeField};
use crate::evaluation_domain::MultiplicativeSubgroup;
use crate::plonk::circuit::{CompiledCircuit};
use crate::plonk::circuit::compiled_circuit::CircuitCompiler;
use crate::plonk::domain::PlonkDomain;

#[derive(Clone, Debug)]
pub struct ArithmeticGate {
    pub left: usize,
    pub right: usize,
    pub output: usize,
}

#[derive(Clone, Debug)]
pub enum Gate {
    Addition (ArithmeticGate),
    Multiplication (ArithmeticGate),
}

#[derive(Clone, Debug)]
pub struct CircuitDescription<F: FftField + PrimeField> {
    variables_count: usize,
    witness: Vec<usize>,
    public_inputs: Vec<usize>,
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
        }
    }

    pub fn public_inputs_count(&self) -> usize {
        self.public_inputs.len()
    }

    pub fn public_inputs(&self) -> &[usize] {
        self.public_inputs.as_slice()
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

    pub fn addition_gate(&mut self, left: usize, right: usize) -> usize {
        let output = self.add_variable();

        self.gates.push(Gate::Addition(ArithmeticGate {
            left,
            right,
            output,
        }));

        output
    }

    pub fn multiplication_gate(&mut self, left: usize, right: usize) -> usize {
        let output = self.add_variable();

        self.gates.push(Gate::Multiplication(ArithmeticGate {
            left,
            right,
            output,
        }));

        output
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

    pub fn compile<'a, 'b>(&'b self, domain: &'a PlonkDomain<'a, F>) -> CompiledCircuit<'a, F> {
        CircuitCompiler::new(self, domain).compile()
    }
}

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::Fr;
    use crate::plonk::circuit::build_test_circuit;
    use crate::plonk::circuit::circuit_description::{ArithmeticGate, Gate};

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
}