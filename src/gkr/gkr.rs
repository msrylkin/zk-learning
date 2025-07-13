use std::collections::{HashMap, HashSet};
use std::fmt::Debug;
use std::ops::{Add, Mul, MulAssign, Neg};
use ark_ff::{Field, One, Zero};
use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use ark_test_curves::bls12_381::Fr;
use crate::sumcheck;
// use crate::sumcheck::{interpolate_univariate_on_evals, OracleEvaluation, SumCheckPoly, SumCheckProof};


// pub fn test_gkr() {
//     // (a + b) * (c * c)
//     let input_a = InputGate {
//         index: 0,
//         value: Fr::from(10)
//     };
//     let input_b = InputGate {
//         index: 1,
//         value: Fr::from(200)
//     };
//     let input_c = InputGate {
//         index: 2,
//         value: Fr::from(20),
//     };
//     let input_d = InputGate {
//         index: 3,
//         value: Fr::from(300),
//     };
// 
//     let add_gate_1 = Gate {
//         inputs: &GateInput::Inputs((&input_a, &input_b)),
//         executor: ExecutorGateEnum::Add(AddGate {})
//     };
//     let add_gate_2 = Gate {
//         inputs: &GateInput::Inputs((&input_c, &input_d)),
//         executor: ExecutorGateEnum::Add(AddGate {})
//     };
//     let mul_gate = Gate {
//         inputs: &GateInput::Gates((&add_gate_1, &add_gate_2)),
//         executor: ExecutorGateEnum::Mul(MulGate {})
//     };
// 
//     let layer_1 = Layer {
//         gates: vec![&add_gate_1, &add_gate_2],
//     };
//     let layer_2 = Layer {
//         gates: vec![&mul_gate],
//     };
// 
//     let circuit = Circuit {
//         layers: vec![layer_1, layer_2],
//     };
// 
//     let solution = circuit.solve();
// 
//     println!("{:?}", solution);
//     
//     let random_points = vec![3,22,12,93,8181,12398,123]
//         .into_iter()
//         .map(Fr::from)
//         .collect::<Vec<_>>();
// 
//     let gkr_proof = prove(&circuit, &solution, &random_points);
//     verify(circuit, &gkr_proof);
// }
