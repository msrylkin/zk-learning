mod circuit;
mod prover;
mod verifier;
mod proof;
mod test_utils;
mod protocol;

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::Fr;
    use crate::gkr::protocol::GKRProtocol;
    use crate::gkr::test_utils::{get_test_circuit};
    use crate::poly_utils::to_f;
    use crate::random_oracle::{FixedRandomOracle};
    
    fn get_test_random_oracle() -> FixedRandomOracle<Fr> {
        let random_points = to_f(vec![3,4,5,9,11,12,13, 14]);
        FixedRandomOracle::new(random_points)
    }

    #[test]
    fn test_prove() {
        let circuit = get_test_circuit();
        let solution = circuit.solve();

        let random_oracle = get_test_random_oracle();
        let gkr_protocol = GKRProtocol::new(&random_oracle);

        let gkr_proof = gkr_protocol.prove(&circuit, &solution);

        assert_eq!(gkr_proof.inputs, to_f(vec![10, 200, 20, 300]));
        assert_eq!(gkr_proof.r0, []);
        assert_eq!(gkr_proof.W0.num_vars, 0);
        assert_eq!(gkr_proof.W0.evaluations, vec![Fr::from(67200)]);

        assert_eq!(gkr_proof.layers.len(), 2);
        
        assert_eq!(gkr_proof.layers[1].q.coeffs, [Fr::from(540), Fr::from(110)]);
        assert_eq!(gkr_proof.layers[1].sumcheck_proof.first_round_poly.coeffs, to_f(vec![67200, -32000, -35200]));
        assert_eq!(gkr_proof.layers[1].sumcheck_proof.steps[0].r, Fr::from(3));
        assert_eq!(gkr_proof.layers[1].sumcheck_proof.steps[0].poly.coeffs, to_f(vec![0, -226800, -118800]));
        assert_eq!(gkr_proof.layers[1].sumcheck_proof.last_round_r, Fr::from(4));
        assert_eq!(gkr_proof.layers[1].r_star, Fr::from(5));

        assert_eq!(gkr_proof.layers[0].q.coeffs, [Fr::from(11100), Fr::from(5000), Fr::from(540)]);
        assert_eq!(gkr_proof.layers[0].sumcheck_proof.first_round_poly.coeffs, [Fr::from(-1470), Fr::from(3880), Fr::from(150)]);
        assert_eq!(gkr_proof.layers[0].sumcheck_proof.steps[0].r, Fr::from(9));
        assert_eq!(gkr_proof.layers[0].sumcheck_proof.steps[0].poly.coeffs, to_f(vec![45600, 82400, -128000]));
        assert_eq!(gkr_proof.layers[0].sumcheck_proof.steps[1].r, Fr::from(11));
        assert_eq!(gkr_proof.layers[0].sumcheck_proof.steps[1].poly.coeffs, to_f(vec![-6328000, -1864000, -16000]));
        assert_eq!(gkr_proof.layers[0].sumcheck_proof.steps[2].r, Fr::from(12));
        assert_eq!(gkr_proof.layers[0].sumcheck_proof.steps[2].poly.coeffs, to_f(vec![0, -27850400, -3149600]));
        assert_eq!(gkr_proof.layers[0].sumcheck_proof.last_round_r, Fr::from(13));
        assert_eq!(gkr_proof.layers[0].r_star, Fr::from(14));
    }
    
    #[test]
    fn test_verify() {
        let circuit = get_test_circuit();
        let solution = circuit.solve();

        let random_oracle = get_test_random_oracle();
        let gkr_protocol = GKRProtocol::new(&random_oracle);

        let gkr_proof = gkr_protocol.prove(&circuit, &solution);
        
        // must not panic
        gkr_protocol.verify(&circuit, &gkr_proof);
    }

    #[test]
    #[should_panic]
    fn test_invalid_proof_inputs() {
        let circuit = get_test_circuit();
        let solution = circuit.solve();

        let random_oracle = get_test_random_oracle();
        let gkr_protocol = GKRProtocol::new(&random_oracle);

        let mut gkr_proof = gkr_protocol.prove(&circuit, &solution);
        
        // mess the inputs
        gkr_proof.inputs[0] = Fr::from(999);

        // must panic
        gkr_protocol.verify(&circuit, &gkr_proof);
    }

    #[test]
    #[should_panic]
    fn test_invalid_proof_output_poly() {
        let circuit = get_test_circuit();
        let solution = circuit.solve();

        let random_oracle = get_test_random_oracle();
        let gkr_protocol = GKRProtocol::new(&random_oracle);

        let mut gkr_proof = gkr_protocol.prove(&circuit, &solution);

        // mess the outputs
        gkr_proof.W0.evaluations[0] = Fr::from(999);

        // must panic
        gkr_protocol.verify(&circuit, &gkr_proof);
    }

    #[test]
    #[should_panic]
    fn test_invalid_proof_sc_poly() {
        let circuit = get_test_circuit();
        let solution = circuit.solve();

        let random_oracle = get_test_random_oracle();
        let gkr_protocol = GKRProtocol::new(&random_oracle);

        let mut gkr_proof = gkr_protocol.prove(&circuit, &solution);

        // mess the some poly coeffs
        gkr_proof.layers[0].sumcheck_proof.steps[1].poly.coeffs[1] = Fr::from(999);

        // must panic
        gkr_protocol.verify(&circuit, &gkr_proof);
    }

    #[test]
    #[should_panic]
    fn test_invalid_q_degree() {
        let circuit = get_test_circuit();
        let solution = circuit.solve();

        let random_oracle = get_test_random_oracle();
        let gkr_protocol = GKRProtocol::new(&random_oracle);

        let mut gkr_proof = gkr_protocol.prove(&circuit, &solution);

        // mess the some poly coeffs
        gkr_proof.layers[0].q.coeffs.push(Fr::from(999));

        // must panic
        gkr_protocol.verify(&circuit, &gkr_proof);
    }
}