use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, MultilinearExtension};
use ark_test_curves::bls12_381::Fr;
use crate::gkr::prover::LayerRoundPoly;
use crate::poly_utils::remap_to_reverse_bits_indexing;

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