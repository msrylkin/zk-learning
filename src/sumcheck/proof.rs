use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

#[derive(Debug)]
pub struct SumCheckStep<F: Field> {
    pub poly: DensePolynomial<F>,
    pub r: F,
}

impl<F: Field> SumCheckStep<F> {
    pub fn new(poly: DensePolynomial<F>, r: F) -> Self {
        SumCheckStep {
            poly,
            r
        }
    }
}

#[derive(Debug)]
pub struct SumCheckProof<F: Field> {
    pub first_round_poly: DensePolynomial<F>,
    pub last_round_r: F,
    pub steps: Vec<SumCheckStep<F>>
}

impl<F: Field> SumCheckProof<F> {
    pub fn new(first_round_poly: DensePolynomial<F>, last_round_r: F, steps: Vec<SumCheckStep<F>>) -> Self {
        SumCheckProof {
            first_round_poly,
            last_round_r,
            steps,
        }
    }
    
    pub fn get_used_randomness(&self) -> Vec<F> {
        let mut used_r = self.steps
            .iter()
            .map(|step| step.r)
            .collect::<Vec<_>>();
        used_r.push(self.last_round_r);

        used_r
    }
}
