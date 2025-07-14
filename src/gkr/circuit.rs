use ark_ff::Field;
use ark_poly::DenseMultilinearExtension;
use crate::poly_utils::{reverse_bits};

#[derive(Debug)]
struct Gate {
    inputs: (usize, usize), // indexes
}

#[derive(Debug)]
struct InputGate<F: Field> {
    value: F,
}

#[derive(Debug)]
enum GateType {
    AddGate(Gate),
    MulGate(Gate),
}

impl GateType {
    fn execute<F: Field>(&self, previous_layer: &[F]) -> F {
        match self {
            GateType::AddGate(Gate { inputs: (i_a, i_b)}) => previous_layer[*i_a] + previous_layer[*i_b],
            GateType::MulGate(Gate { inputs: (i_a, i_b )}) => previous_layer[*i_a] * previous_layer[*i_b],
        }
    }
    
    fn get_inputs(&self) -> (usize, usize) {
        match self {
            GateType::AddGate(Gate { inputs: (i_a, i_b)}) => (*i_a, *i_b),
            GateType::MulGate(Gate { inputs: (i_a, i_b )}) => (*i_a, *i_b),
        }
    }
    
    fn is_add(&self) -> bool {
        if let GateType::AddGate(_) = self {
            return true;
        }
        
        false
    }

    fn is_mul(&self) -> bool {
        if let GateType::MulGate(_) = self {
            return true;
        }

        false
    }
}

#[derive(Debug)]
pub struct Layer {
    gates: Vec<GateType>,
}

impl Layer {
    pub fn new() -> Self {
        Self { gates: vec![] }
    }

    pub fn add_addition_gate(self, inputs: (usize, usize)) -> Self {
        let mut gates = self.gates;
        gates.push(GateType::AddGate(Gate { inputs }));

        Self {
            gates,
            ..self
        }
    }

    pub fn add_multiplication_gate(self, inputs: (usize, usize)) -> Self {
        let mut gates = self.gates;
        gates.push(GateType::MulGate(Gate { inputs }));
        
        Self {
            gates,
            ..self
        }
    }
}

#[derive(Debug)]
pub struct Circuit<F: Field> {
    inputs: Vec<F>,
    layers: Vec<Layer>
}


pub struct CircuitBuilder<F: Field> {
    circuit: Circuit<F>,
}

impl<F: Field> CircuitBuilder<F> {
    pub fn new(inputs: Vec<F>) -> Self {
        assert_ne!(inputs.len(), 0);

        Self {
            circuit: Circuit {
                inputs,
                layers: vec![],
            }
        }
    }

    pub fn add_layer(self, layer: Layer) -> Self {
        let mut circuit = self.circuit;
        circuit.layers.push(layer);
        
        Self {
            circuit
        }
    }
    
    pub fn build(self) -> Circuit<F> {
        self.circuit
    }
}

impl<F: Field> Solution<F>  {
    pub fn to_evaluations(&self) -> Vec<Vec<F>> {
        self.evaluations.clone()
    }
    
    pub fn inputs(&self) -> Vec<F> {
        self.inputs.clone()
    }
}

#[derive(Debug)]
pub struct Solution<F: Field> {
    evaluations: Vec<Vec<F>>,
    inputs: Vec<F>,
}

impl<F: Field> Circuit<F> {
    pub fn solve(&self) -> Solution<F> {
        let mut evaluations = vec![];
        for layer in &self.layers {
            let previous_layer_values = evaluations.last().unwrap_or(&self.inputs);
            let mut current_layer_values = vec![];
            
            for gate in &layer.gates {
                let calculated_result = gate.execute(&previous_layer_values);
                current_layer_values.push(calculated_result);
            }
            
            evaluations.push(current_layer_values);
        }
        
        Solution {
            inputs: self.inputs.clone(),
            evaluations,
        }
    }
    
    fn eq<P: Fn(&GateType) -> bool>(
        &self,
        layer_i: usize,
        check_gate: P,
    ) -> DenseMultilinearExtension<F> {
        let layer = &self.layers[layer_i];
        let gates_n = layer.gates.len();
        let vars_num = f64::from(gates_n as u32).log2().ceil() as usize;

        let bottom_gates_n = self.get_bottom_layer_gates_count(layer_i);
        let bottom_vars_num = f64::from(bottom_gates_n as u32).log2().ceil() as usize;

        let total_vars_num = vars_num + bottom_vars_num * 2;

        let mut evals = vec![F::zero(); 1 << total_vars_num];

        for (i_z, gate) in layer.gates.iter().enumerate() {
            if check_gate(gate) {
                let (i_b, i_c) = gate.get_inputs();
                let z_mask = i_z << (bottom_vars_num * 2);
                let b_mask = i_b << bottom_vars_num;
                let c_mask = i_c;
                let zbc_mask = z_mask | b_mask | c_mask;

                evals[reverse_bits(zbc_mask, total_vars_num)] = F::one();
            }
        }

        DenseMultilinearExtension::from_evaluations_vec(total_vars_num, evals)
    }

    pub fn get_bottom_layer_gates_count(&self, layer_i: usize) -> usize {
        match layer_i {
            0 => self.inputs.len(),
            _ => self.layers[layer_i - 1].gates.len(),
        }
    }
    
    pub fn add_i(&self, layer_i: usize) -> DenseMultilinearExtension<F> {
        self.eq(layer_i, |gate| matches!(gate, GateType::AddGate(_)))
    }

    pub fn mul_i(&self, layer_i: usize) -> DenseMultilinearExtension<F> {
        self.eq(layer_i, |gate| matches!(gate, GateType::MulGate(_)))
    }
}

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::Fr;
    use crate::gkr::test_utils::get_test_circuit;
    
    #[test]
    fn circuit_test() {
        let circuit = get_test_circuit();
        
        let solution = circuit.solve();
        assert_eq!(solution.inputs, vec![10, 200, 20, 300].into_iter().map(Fr::from).collect::<Vec<_>>());
        assert_eq!(solution.evaluations, vec![
            vec![210, 320].into_iter().map(Fr::from).collect::<Vec<_>>(),
            vec![67200].into_iter().map(Fr::from).collect::<Vec<_>>(),
        ]);
        
        let add_i = circuit.add_i(0);
        assert_eq!(
            add_i.evaluations,
            vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0].into_iter().map(Fr::from).collect::<Vec<_>>(),
        );
        assert_eq!(add_i.num_vars, 5);

        let add_i = circuit.add_i(1);
        assert_eq!(
            add_i.evaluations,
            vec![0, 0, 0, 0].into_iter().map(Fr::from).collect::<Vec<_>>(),
        );
        assert_eq!(add_i.num_vars, 2);
        
        let mul_i = circuit.mul_i(0);
        assert_eq!(
            mul_i.evaluations,
            vec![0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0].into_iter().map(Fr::from).collect::<Vec<_>>(),
        );
        assert_eq!(mul_i.num_vars, 5);

        let mul_i = circuit.mul_i(1);
        assert_eq!(
            mul_i.evaluations,
            vec![0, 0, 1, 0].into_iter().map(Fr::from).collect::<Vec<_>>(),
        );
        assert_eq!(mul_i.num_vars, 2);
    }
}