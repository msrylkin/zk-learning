use ark_ff::Field;
use ark_poly::univariate::DensePolynomial;

/// A single round in the sumcheck protocol.
#[derive(Debug)]
pub struct SumCheckStep<F: Field> {
    /// Partial sum polynomial, summed over the Boolean domain except for the first variable.
    pub poly: DensePolynomial<F>,
    /// Randomness used to fix the first variable.
    pub r: F,
}

impl<F: Field> SumCheckStep<F> {
    /// Creates a new sumcheck step from the given `poly` and `r`.
    pub fn new(poly: DensePolynomial<F>, r: F) -> Self {
        SumCheckStep {
            poly,
            r
        }
    }
}

/// The main sumcheck proof representation.
#[derive(Debug)]
pub struct SumCheckProof<F: Field> {
    /// First-round partial sum polynomial with no variables fixed.
    pub first_round_poly: DensePolynomial<F>,
    /// Randomness used in the final round.
    pub last_round_r: F,
    /// All intermediate steps in the sumcheck protocol.
    ///
    /// In each round, the prover samples randomness, fixes the first (leftmost) variable,
    /// and sends the resulting partial sum polynomial to the verifier.
    pub steps: Vec<SumCheckStep<F>>
}

impl<F: Field> SumCheckProof<F> {
    /// Creates a new sumcheck proof.
    pub fn new(first_round_poly: DensePolynomial<F>, last_round_r: F, steps: Vec<SumCheckStep<F>>) -> Self {
        SumCheckProof {
            first_round_poly,
            last_round_r,
            steps,
        }
    }

    /// Returns all randomness values `r` used in the protocol:
    /// the per-round randomness values, followed by the final one.
    pub fn get_used_randomness(&self) -> Vec<F> {
        let mut used_r = self.steps
            .iter()
            .map(|step| step.r)
            .collect::<Vec<_>>();
        used_r.push(self.last_round_r);

        used_r
    }
}
