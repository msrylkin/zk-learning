pub(crate) mod gkr;
mod circuit;
mod prover;
mod verifier;
mod common;
mod test_utils;

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::Fr;
    use ark_ff::{Field, Zero};
    use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension};
    use ark_poly::univariate::DensePolynomial;
    use crate::gkr::common::GKRProofLayer;
    use crate::gkr::prover::{prove, LayerRoundPoly};
    use crate::gkr::test_utils::get_test_circuit;
    use crate::poly_utils::{get_evaluations_by_mask, remap_to_reverse_bits_indexing, to_two_or_one_degree};
    use crate::sumcheck::{prove as sc_prove, SumCheckPoly, SumCheckProof};
    use super::*;

    #[test]
    fn test_prove() {
        let circuit = get_test_circuit();
        let solution = circuit.solve();
        let random_points = vec![3,22,12,93,8181,12398,123]
            .into_iter()
            .map(Fr::from)
            .collect::<Vec<_>>();

        let gkr_proof = prove(&circuit, &solution, &random_points);
        
        println!("{:?}", gkr_proof);
        
        assert_eq!(gkr_proof.inputs, vec![10, 200, 20, 300].into_iter().map(Fr::from).collect::<Vec<_>>());
        assert_eq!(gkr_proof.outputs, vec![Fr::from(67200)]);
        assert_eq!(gkr_proof.r0, []);
        assert_eq!(gkr_proof.W0.num_vars, 0);
        assert_eq!(gkr_proof.W0.evaluations, vec![Fr::from(67200)]);
        
        assert_eq!(gkr_proof.layers.len(), 2);
        // layer 1
        // assert_eq!(gkr_proof.layers[0].r_star, Fr::from(22));
        // assert_eq!(gkr_proof.layers[0], GKRProofLayer {
        //     sumcheck_proof: SumCheckProof {
        //         first_round_poly: DensePolynomial::from_coefficients_vec(vec![]),
        //         last_round_r: Fr::from(123),
        //         steps: vec![],
        //     },
        //     r_star: Fr::from(123),
        //     q: DensePolynomial::from_coefficients_vec(vec![])
        // })
    }
}