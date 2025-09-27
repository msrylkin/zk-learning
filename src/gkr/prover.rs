use ark_ff::Field;
use ark_poly::{DenseMultilinearExtension, DenseUVPolynomial, MultilinearExtension, Polynomial};
use ark_poly::univariate::DensePolynomial;
use crate::gkr::circuit::{Circuit, Solution};
use crate::gkr::proof::{GKRProof, GKRProofLayer};
use crate::random_oracle::RandomOracle;
use crate::poly_utils::{get_bits, interpolate, line, restrict_poly, reverse_bits, to_two_or_one_degree};
use crate::sumcheck::{SumCheckPoly, SumCheckProtocol};

#[derive(Debug, Clone)]
pub struct LayerRoundPoly<F: Field> {
    add_i: DenseMultilinearExtension<F>,
    mul_i: DenseMultilinearExtension<F>,
    Wi_1_a: DenseMultilinearExtension<F>,
    Wi_1_b: DenseMultilinearExtension<F>,
}

type FixingCounts = (usize, usize);
type EvaluationMasks = (usize, usize, usize);

impl<F: Field> LayerRoundPoly<F> {
    pub fn new(
        add_i: DenseMultilinearExtension<F>,
        mul_i: DenseMultilinearExtension<F>,
        Wi_1: DenseMultilinearExtension<F>,
    ) -> Self {
        Self {
            add_i,
            mul_i,
            Wi_1_a: Wi_1.clone(),
            Wi_1_b: Wi_1,
        }
    }
}

fn generate_round_poly_eval_parameters(
    a_nvars: usize,
    b_nvars: usize,
) -> (FixingCounts, Vec<EvaluationMasks>){
    let (a_vars_fixing, b_vars_fixing) = {
        if a_nvars > 0 {
            (a_nvars - 1, b_nvars)
        } else {
            (a_nvars, b_nvars - 1)
        }
    };

    let combined_domain = 1 << (a_vars_fixing + b_vars_fixing); // for example 1111 where [1] = a and [111] = b
    let b_mask = (1 << b_vars_fixing) - 1; // 0111
    let a_mask = (1 << (a_vars_fixing + b_vars_fixing)) - 1 - b_mask; // 1000

    let mut res = vec![];

    for i in 0..combined_domain {
        let combined_param = i;
        let b_param = i & b_mask;
        let a_param = (i & a_mask) >> b_vars_fixing;

        res.push((combined_param, a_param, b_param))
    }

    ((a_vars_fixing, b_vars_fixing), res)
}

impl<F: Field> SumCheckPoly<F> for LayerRoundPoly<F> {
    fn get_evaluations(&self) -> Vec<F> {
        let mut res = vec![];
        let a_num_vars = self.Wi_1_a.num_vars();
        let b_num_vars = self.Wi_1_b.num_vars();
        let ab_num_vars = a_num_vars + b_num_vars;
        let a_domain = 1 << a_num_vars;
        let b_domain = 1 << b_num_vars;

        for a in 0..a_domain {
            for b in 0..b_domain {
                let ab_combined_mask = (a << b_num_vars) | b;
                let ab_combined_mask_reversed = reverse_bits(ab_combined_mask, ab_num_vars);
                let ab_bits = get_bits(ab_combined_mask_reversed, ab_num_vars).into_iter().map(|e| F::from(e as u64)).collect::<Vec<_>>();
                let eval = self.evaluate(&ab_bits);

                res.push(eval);
            }
        }

        res
    }

    fn num_vars(&self) -> usize {
        self.add_i.num_vars()
    }

    fn get_partial_sum_poly(&self) -> DensePolynomial<F> {
        let a_vars = self.Wi_1_a.num_vars();
        let b_vars = self.Wi_1_b.num_vars();

        let mut res_poly = DensePolynomial::from_coefficients_vec(vec![F::zero()]);

        let ((a_vars_fixing, b_vars_fixing), eval_params) = generate_round_poly_eval_parameters(a_vars, b_vars);

        for (common_params, a_params, b_params) in eval_params {
            let Wi_a_uni = to_two_or_one_degree(&self.Wi_1_a, a_params, a_vars_fixing);
            let Wi_b_uni = to_two_or_one_degree(&self.Wi_1_b, b_params, b_vars_fixing);
            let add_i_uni = to_two_or_one_degree(&self.add_i, common_params, a_vars_fixing + b_vars_fixing);
            let mul_uni = to_two_or_one_degree(&self.mul_i, common_params, a_vars_fixing + b_vars_fixing);

            let add_poly = add_i_uni.naive_mul(&(&Wi_a_uni + &Wi_b_uni));
            let mul_poly = mul_uni.naive_mul(&(Wi_a_uni.naive_mul(&Wi_b_uni)));
            let local_res_poly = add_poly + mul_poly;

            res_poly = res_poly + local_res_poly;
        }

        res_poly
    }

    fn fix_variable(&self, e: F) -> Self {
        if self.Wi_1_a.num_vars() > 0 {
            Self {
                add_i: self.add_i.fix_variables(&[e]),
                mul_i: self.mul_i.fix_variables(&[e]),
                Wi_1_a: self.Wi_1_a.fix_variables(&[e]),
                Wi_1_b: self.Wi_1_b.clone(),
            }
        } else {
            Self {
                add_i: self.add_i.fix_variables(&[e]),
                mul_i: self.mul_i.fix_variables(&[e]),
                Wi_1_a: self.Wi_1_a.clone(),
                Wi_1_b: self.Wi_1_b.fix_variables(&[e]),
            }
        }
    }

    fn evaluate(&self, point: &[F]) -> F {
        let binding = point.to_vec();
        let a_vars = self.Wi_1_a.num_vars();
        let (a_bits, b_bits) = binding.split_at(a_vars);
        let a_bits = a_bits.to_vec();
        let b_bits = b_bits.to_vec();
        
        let add_eval = self.add_i.evaluate(&point.to_vec()) * (self.Wi_1_a.evaluate(&a_bits) + self.Wi_1_b.evaluate(&b_bits));
        let mul_eval = self.mul_i.evaluate(&point.to_vec()) * (self.Wi_1_a.evaluate(&a_bits) * self.Wi_1_b.evaluate(&b_bits));

        add_eval + mul_eval
    }
}

pub fn prove<F: Field, O: RandomOracle<Item = F>>(
    circuit: &Circuit<F>,
    solution: &Solution<F>,
    random_oracle: &O,
    sum_check_protocol: &SumCheckProtocol<O>,
) -> GKRProof<F> {
    let solution_evaluations = solution.to_evaluations();
    let solution_inputs = solution.inputs();
    let outputs = solution_evaluations.last().unwrap();

    let W0 = interpolate(outputs);
    let r0 = random_oracle.get_randomness(W0.num_vars());
    let mut ri = r0.clone();

    let mut gkr_proof_layers = vec![];

    for i in (0..solution_evaluations.len()).rev() {
        let add_poly = circuit.add_i(i);
        let mul_poly = circuit.mul_i(i);
        let Wi_1 = {
            if i == 0 {
                interpolate(&solution_inputs)
            } else {
                interpolate(&solution_evaluations[i - 1])
            }
        };

        let add_fixed = add_poly.fix_variables(&ri);
        let mul_fixed = mul_poly.fix_variables(&ri);

        let sc_poly = LayerRoundPoly::new(
            add_fixed,
            mul_fixed,
            Wi_1.clone()
        );

        let sumcheck_proof = sum_check_protocol.prove(&sc_poly);

        let used_r = sumcheck_proof.get_used_randomness();

        let (b, c) = used_r.split_at(used_r.len() / 2);
        let l = line(b, c);

        let q = restrict_poly(&l, Wi_1);
        let r_star = random_oracle.get_randomness(1)[0];

        ri = l.iter().map(|li| li.evaluate(&r_star)).collect();

        gkr_proof_layers.push(GKRProofLayer {
            q,
            sumcheck_proof,
            r_star,
        });
    }

    GKRProof {
        W0,
        inputs: solution_inputs,
        layers: gkr_proof_layers.into_iter().rev().collect(),
        r0,
    }
}

#[cfg(test)]
mod tests {
    use ark_ff::Zero;
    use ark_test_curves::bls12_381::Fr;
    use crate::gkr::test_utils::{get_test_round_poly_2_vars, get_test_round_poly_4_vars};
    use super::*;

    #[test]
    fn generate_round_poly_params() {
        let a_vars = 2;
        let b_vars = 2;

        let ((a_vars_fixing, b_vars_fixing), res) = generate_round_poly_eval_parameters(a_vars, b_vars);

        assert_eq!(res, [
            (0b000, 0b000, 0b000),
            (0b001, 0b000, 0b001),
            (0b010, 0b000, 0b010),
            (0b011, 0b000, 0b011),
            (0b100, 0b001, 0b000),
            (0b101, 0b001, 0b001),
            (0b110, 0b001, 0b010),
            (0b111, 0b001, 0b011),
        ]);
        assert_eq!(a_vars_fixing, 1);
        assert_eq!(b_vars_fixing, 2);
    }

    #[test]
    fn generate_round_poly_params_one_var() {
        let a_vars = 1;
        let b_vars = 2;

        let ((a_vars_fixing, b_vars_fixing), res) = generate_round_poly_eval_parameters(a_vars, b_vars);

        assert_eq!(res, [
            (0b000, 0b000, 0b000),
            (0b001, 0b000, 0b001),
            (0b010, 0b000, 0b010),
            (0b011, 0b000, 0b011),
        ]);
        assert_eq!(a_vars_fixing, 0);
        assert_eq!(b_vars_fixing, 2);
    }

    #[test]
    fn generate_round_poly_params_zero_var() {
        let a_vars = 0;
        let b_vars = 2;

        let ((a_vars_fixing, b_vars_fixing), res) = generate_round_poly_eval_parameters(a_vars, b_vars);

        assert_eq!(res, [
            (0b000, 0b000, 0b000),
            (0b001, 0b000, 0b001),
        ]);
        assert_eq!(a_vars_fixing, 0);
        assert_eq!(b_vars_fixing, 1);
    }

    #[test]
    fn round_poly_fix_polys() {
        let round_poly = get_test_round_poly_4_vars::<Fr>();
        let mask = 0b0001;
        let a_mask = 0b0000;
        let b_mask = 0b0001;
    
        let add_uni = to_two_or_one_degree(&round_poly.add_i, mask, 3);
        let mul_uni = to_two_or_one_degree(&round_poly.mul_i, mask, 3);
        let Wi_1_a_uni = to_two_or_one_degree(&round_poly.Wi_1_a, a_mask, 1);
        let Wi_1_b_uni = to_two_or_one_degree(&round_poly.Wi_1_b, b_mask, 2);
    
        assert_eq!(add_uni.coeffs, [Fr::from(-2), Fr::from(2)]);
        assert_eq!(mul_uni.coeffs, []);
        assert_eq!(Wi_1_a_uni.coeffs, [Fr::from(10), Fr::from(10)]);
        assert_eq!(Wi_1_b_uni.coeffs, [Fr::from(200)]);
    }
    
    #[test]
    fn round_poly_fix_polys_other_mask() {
        let round_poly = get_test_round_poly_4_vars::<Fr>();
        let mask = 0b0101;
        let a_mask = 0b0001;
        let b_mask = 0b0001;
    
        let add_uni = to_two_or_one_degree(&round_poly.add_i, mask, 3);
        let mul_uni = to_two_or_one_degree(&round_poly.mul_i, mask, 3);
        let Wi_1_a_uni = to_two_or_one_degree(&round_poly.Wi_1_a, a_mask, 1);
        let Wi_1_b_uni = to_two_or_one_degree(&round_poly.Wi_1_b, b_mask, 2);
    
        assert_eq!(add_uni.coeffs, []);
        assert_eq!(mul_uni.coeffs, []);
        assert_eq!(Wi_1_a_uni.coeffs, [Fr::from(200), Fr::from(100)]);
        assert_eq!(Wi_1_b_uni.coeffs, [Fr::from(200)]);
    }

    #[test]
    fn get_partial_sum_poly_2_vars() {
        let round_poly = get_test_round_poly_2_vars::<Fr>();
        let partial_sum_poly_1 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_1.coeffs, [Fr::from(67200), Fr::from(-32000), Fr::from(-35200)]);

        let round_poly = round_poly.fix_variable(Fr::from(32));
        let partial_sum_poly_2 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_2.coeffs, [Fr::zero(), Fr::from(-24282300), Fr::from(-12719300)])
    }

    #[test]
    fn get_partial_sum_poly_4_vars() {
        let round_poly = get_test_round_poly_4_vars::<Fr>();
        let partial_sum_poly_1 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_1.coeffs, [Fr::from(-420), Fr::from(1330), Fr::from(50)]);

        let round_poly = round_poly.fix_variable(Fr::from(4));
        let partial_sum_poly_2 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_2.coeffs, [Fr::from(5700), Fr::from(4200), Fr::from(-9900)]);

        let round_poly = round_poly.fix_variable(Fr::from(5));
        let partial_sum_poly_3 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_3.coeffs, [Fr::from(-72000), Fr::from(-74400), Fr::from(-2400)]);

        let round_poly = round_poly.fix_variable(Fr::from(8));
        let partial_sum_poly_4 = round_poly.get_partial_sum_poly();

        assert_eq!(partial_sum_poly_4.coeffs, [Fr::from(0), Fr::from(-624240), Fr::from(-196560)]);
    }

    #[test]
    fn round_poly_get_required_params() {
        let round_poly = get_test_round_poly_4_vars::<Fr>();
        let ((a_eval_params, b_eval_params), evals) = generate_round_poly_eval_parameters(
            round_poly.Wi_1_a.num_vars,
            round_poly.Wi_1_b.num_vars,
        );

        assert_eq!(a_eval_params, 1);
        assert_eq!(b_eval_params, 2);
        assert_eq!(evals[1], (0b0001usize, 0b0000usize, 0b0001usize));
    }
}