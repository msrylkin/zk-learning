mod protocol;
mod prover;
mod verifier;
mod sumcheck_poly;
mod multilin_sumcheck_poly;
mod proof;

pub use proof::*;
pub use sumcheck_poly::*;
pub use protocol::*;

#[cfg(test)]
mod tests {
    use ark_poly::{DenseUVPolynomial};
    use ark_poly::univariate::DensePolynomial;
    use ark_test_curves::bls12_381::Fr;
    use crate::poly_utils::to_f;
    use crate::random_oracle::FixedRandomOracle;
    use crate::sumcheck::multilin_sumcheck_poly::{DenseMultilinFinalEvaluationOracle, DenseMultilinSumcheckPoly};
    use crate::sumcheck::{SumCheckProof, SumCheckProtocol};

    fn get_test_random_oracle() -> FixedRandomOracle<Fr> {
        let random_points = to_f(vec![3,4,5,9,11,12,13, 14]);
        FixedRandomOracle::new(random_points)
    }

    #[test]
    pub fn test_sumcheck_prove() {
        let random_oracle = get_test_random_oracle();
        let sumcheck_protocol = SumCheckProtocol::new(1, &random_oracle);
        let sumcheck_multilin_poly = DenseMultilinSumcheckPoly::from_evaluations_vec(
            3,
            to_f::<Fr>(vec![10, 55, 300, 77, 20, 66, 400, 88])
        );
        
        let sumcheck_proof = sumcheck_protocol.prove(&sumcheck_multilin_poly);
        
        assert_eq!(sumcheck_proof.first_round_poly.coeffs, to_f::<Fr>(vec![730, -444]));
        assert_eq!(sumcheck_proof.steps[0].poly.coeffs, to_f::<Fr>(vec![303, -1208]));
        assert_eq!(sumcheck_proof.steps[0].r, Fr::from(3));
        assert_eq!(sumcheck_proof.steps[1].poly.coeffs, to_f::<Fr>(vec![-1911, -707]));
        assert_eq!(sumcheck_proof.steps[1].r, Fr::from(4));
        assert_eq!(sumcheck_proof.last_round_r, Fr::from(5));
    }
    
    #[test]
    pub fn test_verify() {
        let random_oracle = get_test_random_oracle();
        let sumcheck_protocol = SumCheckProtocol::new(1, &random_oracle);
        let sumcheck_multilin_poly = DenseMultilinSumcheckPoly::from_evaluations_vec(
            3,
            to_f::<Fr>(vec![10, 55, 300, 77, 20, 66, 400, 88])
        );
        
        let claimed_sum = sumcheck_multilin_poly.evaluations.iter().sum::<Fr>();

        let sumcheck_proof = sumcheck_protocol.prove(&sumcheck_multilin_poly);
        let final_evaluation_oracle = DenseMultilinFinalEvaluationOracle::new(sumcheck_multilin_poly);
        
        sumcheck_protocol.verify(
            &final_evaluation_oracle,
            &sumcheck_proof,
            claimed_sum,
        );
    }

    #[test]
    #[should_panic]
    pub fn test_verify_wrong_sum() {
        let random_oracle = get_test_random_oracle();
        let sumcheck_protocol = SumCheckProtocol::new(1, &random_oracle);
        let sumcheck_multilin_poly = DenseMultilinSumcheckPoly::from_evaluations_vec(
            3,
            to_f::<Fr>(vec![10, 55, 300, 77, 20, 66, 400, 88])
        );

        let sumcheck_proof = sumcheck_protocol.prove(&sumcheck_multilin_poly);
        let final_evaluation_oracle = DenseMultilinFinalEvaluationOracle::new(sumcheck_multilin_poly);
        
        sumcheck_protocol.verify(
            &final_evaluation_oracle,
            &sumcheck_proof,
            Fr::from(123),
        );
    }

    #[test]
    #[should_panic]
    pub fn test_verify_wrong_degree() {
        let random_oracle = get_test_random_oracle();
        let sumcheck_protocol = SumCheckProtocol::new(1, &random_oracle);
        let sumcheck_multilin_poly = DenseMultilinSumcheckPoly::from_evaluations_vec(
            3,
            to_f::<Fr>(vec![10, 55, 300, 77, 20, 66, 400, 88])
        );

        let claimed_sum = sumcheck_multilin_poly.evaluations.iter().sum::<Fr>();

        let sumcheck_proof = sumcheck_protocol.prove(&sumcheck_multilin_poly);
        let fake_first_round_poly = DensePolynomial::from_coefficients_vec(
            vec![Fr::from(286), Fr::from(-444), Fr::from(888)]
        );
        
        assert_ne!(sumcheck_proof.first_round_poly.coeffs.iter().sum::<Fr>(), claimed_sum);

        let sumcheck_proof = SumCheckProof {
            first_round_poly: fake_first_round_poly,
            ..sumcheck_proof
        };
        let final_evaluation_oracle = DenseMultilinFinalEvaluationOracle::new(sumcheck_multilin_poly);
        
        sumcheck_protocol.verify(
            &final_evaluation_oracle,
            &sumcheck_proof,
            claimed_sum,
        );
    }
}