use ark_ff::Field;
use ark_poly::DenseMultilinearExtension;
use ark_poly::univariate::DensePolynomial;
use crate::sumcheck::SumCheckProof;

#[derive(Debug)]
pub struct GKRProofLayer<F: Field> {
    pub sumcheck_proof: SumCheckProof<F>,
    pub r_star: F,
    pub q: DensePolynomial<F>,
}

#[derive(Debug)]
pub struct GKRProof<F: Field> {
    pub layers: Vec<GKRProofLayer<F>>,
    pub outputs: Vec<F>,
    pub inputs: Vec<F>,
    pub r0: Vec<F>,
    pub W0: DenseMultilinearExtension<F>,
}