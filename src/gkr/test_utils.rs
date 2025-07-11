use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_test_curves::bls12_381::Fr;
use crate::gkr::circuit::{Circuit, CircuitBuilder, Layer};
use crate::gkr::prover::LayerRoundPoly;
use crate::gkr::random_oracle::FixedRandomOracle;
use crate::poly_utils::remap_to_reverse_bits_indexing;

#[cfg(test)]
pub fn get_test_round_poly_2_vars<F: Field>() -> LayerRoundPoly<F> {
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
    // LayerRoundPoly {
    //     add_i,
    //     mul_i,
    //     Wi_1_a: Wi_1.clone(),
    //     Wi_1_b: Wi_1.clone(),
    // }

    LayerRoundPoly::new(add_i, mul_i, Wi_1.clone(), Wi_1.clone())
}

#[cfg(test)]
pub fn get_test_round_poly_4_vars<F: Field>() -> LayerRoundPoly<F> {
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

    // LayerRoundPoly {
    //     add_i: MultilinearExtension::fix_variables(&add_i, &[F::from(3)]),
    //     mul_i: MultilinearExtension::fix_variables(&mul_i, &[F::from(3)]),
    //     Wi_1_a: Wi_1.clone(),
    //     Wi_1_b: Wi_1.clone(),
    // }
    LayerRoundPoly::new(
        MultilinearExtension::fix_variables(&add_i, &[F::from(3)]),
        MultilinearExtension::fix_variables(&mul_i, &[F::from(3)]),
        Wi_1.clone(),
        Wi_1.clone(),
    )
}

#[cfg(test)]
pub fn get_test_circuit() -> Circuit<Fr> {
    // let circuit = Circuit {
    //     inputs: vec![10, 200, 20, 300].into_iter().map(Fr::from).collect(),
    //     layers: vec![
    //         Layer {
    //             gates: vec![
    //                 crate::gkr::circuit::GateType::AddGate(crate::gkr::circuit::Gate {
    //                     inputs: (0, 1),
    //                 }),
    //                 crate::gkr::circuit::GateType::AddGate(crate::gkr::circuit::Gate {
    //                     inputs: (2, 3),
    //                 }),
    //             ]
    //         },
    //         crate::gkr::circuit::Layer {
    //             gates: vec![
    //                 crate::gkr::circuit::GateType::MulGate(crate::gkr::circuit::Gate {
    //                     inputs: (0, 1),
    //                 }),
    //             ]
    //         },
    //     ],
    // };
    
    let circuit = CircuitBuilder::new(
        vec![10, 200, 20, 300].into_iter().map(Fr::from).collect()
    )
        .add_layer(
            Layer::new()
                .add_addition_gate((0, 1))
                .add_addition_gate((2, 3))
        )
        .add_layer(
            Layer::new()
                .add_multiplication_gate((0, 1))
        )
        .build();

    circuit
}

pub fn test_fixed_random_oracle<F: Field>() -> FixedRandomOracle<F> {
    let test_values = vec![3,22,12,93,8181,12398,123].into_iter().map(F::from).collect::<Vec<_>>();
    
    FixedRandomOracle::new(test_values)
}