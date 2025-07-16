#![allow(unused_imports)]

use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_test_curves::bls12_381::Fr;
use crate::gkr::circuit::{Circuit, CircuitBuilder, Layer};
use crate::gkr::prover::LayerRoundPoly;

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

    LayerRoundPoly::new(add_i, mul_i, Wi_1)
}

#[cfg(test)]
pub fn get_test_round_poly_4_vars<F: Field>() -> LayerRoundPoly<F> {
    let add_i = DenseMultilinearExtension::from_evaluations_vec(
        5,
        (0..32).into_iter().map(|e| F::from((e == 16 || e == 27) as u64)).collect::<Vec<_>>(),
    );
    let mul_i = DenseMultilinearExtension::from_evaluations_vec(
        5,
        vec![F::from(0); 32],
    );
    let Wi_1 = DenseMultilinearExtension::from_evaluations_vec(
        2,
        vec![10,20,200,300].into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>(),
    );

    LayerRoundPoly::new(
        add_i.fix_variables(&[F::from(3)]),
        mul_i.fix_variables(&[F::from(3)]),
        Wi_1,
    )
}

#[cfg(test)]
pub fn get_test_circuit() -> Circuit<Fr> {
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