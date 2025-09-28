use std::ops::Deref;
use ark_ec::CurveGroup;
use ark_ec::pairing::Pairing;
use ark_ff::{Field};
use ark_poly::{DenseUVPolynomial, Polynomial};
use ark_poly::univariate::{DenseOrSparsePolynomial, DensePolynomial};
use ark_std::iterable::Iterable;
use ark_std::{One, Zero};
use crate::plonk::big_quotient_poly::BigQuotientPoly;
use crate::plonk::blinder::{blind_big_quotient, blind_solution, blind_z_poly};
use crate::plonk::circuit::PlonkSolution;
use crate::plonk::permutation::PermutationArgument;
use crate::plonk::proof::{Commitments, OpeningProofs, Openings, Proof};
use crate::plonk::protocol::Party;
use crate::plonk::transcript_protocol::TranscriptProtocol;
use crate::poly_utils::{const_poly};

pub struct PlonkProver<'a, P: Pairing> {
    party_data: Party<'a, P>,
}

impl<'a, P: Pairing> Deref for PlonkProver<'a, P> {
    type Target = Party<'a, P>;

    fn deref(&self) -> &Self::Target {
        &self.party_data
    }
}

impl<'a, P: Pairing> PlonkProver<'a, P> {
    pub fn new(
        party_data: Party<'a, P>,
    ) -> Self {
        Self {
            party_data
        }
    }

    pub fn prove(
        &self,
        solution: PlonkSolution<P::ScalarField>,
    ) -> Proof<P> {
        let omega = self.domain.generator();

        let mut transcript = TranscriptProtocol::<P>::new(
            omega,
            [
                solution.public_witness.inputs_vector.clone(),
                solution.public_witness.output_vector.clone(),
            ].concat().as_slice()
        );
        let Zh = &self.Zh;

        let solution = blind_solution(solution, Zh);

        let (a_comm, b_comm, c_comm) = (
            self.kzg.commit(&solution.a),
            self.kzg.commit(&solution.b),
            self.kzg.commit(&solution.c),
        );

        transcript.append_abc_commitments(a_comm.into_affine(), b_comm.into_affine(), c_comm.into_affine());
        let (beta, gamma ) = transcript.get_beta_gamma();

        let permutation_argument = PermutationArgument::new(self.domain, beta, gamma, self.circuit, &solution);

        let z_poly = permutation_argument.z_poly();
        let z_poly = blind_z_poly(&z_poly, Zh);
        let z_comm = self.kzg.commit(&z_poly);
        let z_shifted = shift_poly(&z_poly, &omega);

        transcript.append_z_commitment(z_comm.into_affine());

        let big_t = self.compute_big_quotient(&solution, &z_poly, &z_shifted, &permutation_argument, &transcript);

        let big_t= blind_big_quotient(big_t, self.domain.len());
        let (lo_t, mid_t, hi_t) = big_t.get_splitted_polys();
        let (t_lo_comm, t_mid_comm, t_hi_comm) = (
            self.kzg.commit(lo_t),
            self.kzg.commit(mid_t),
            self.kzg.commit(hi_t),
        );

        transcript.append_t(t_lo_comm.into_affine(), t_mid_comm.into_affine(), t_hi_comm.into_affine());
        let zeta = transcript.get_zeta();

        let openings = self.compute_openings(&solution, &z_shifted, &transcript);

        let r_poly = self.linearization_poly(
            &openings,
            &solution,
            &permutation_argument,
            &z_poly,
            &big_t,
            &transcript,
        );

        assert_eq!(r_poly.evaluate(&zeta), P::ScalarField::zero());

        transcript.append_openings(&openings);
        let vi = transcript.get_vi();

        let zeta_openings = self.kzg.batch_open(
            &[&r_poly, &solution.a, &solution.b, &solution.c, &self.circuit.s_sigma_1, &self.circuit.s_sigma_2],
            &zeta,
            &vi,
        );
        let zeta_omega_opening = self.kzg.open(
            &z_poly,
            &(zeta * omega),
        );

        Proof::new(
            openings,
            Commitments {
                a: a_comm,
                b: b_comm,
                c: c_comm,
                z: z_comm,
                t_lo: t_lo_comm,
                t_mid: t_mid_comm,
                t_hi: t_hi_comm,
            },
            OpeningProofs {
                w_zeta: zeta_openings.proof(),
                w_zeta_omega: zeta_omega_opening.proof(),
            }
        )
    }

    fn compute_big_quotient(
        &self,
        solution: &PlonkSolution<P::ScalarField>,
        z: &DensePolynomial<P::ScalarField>,
        z_shifted: &DensePolynomial<P::ScalarField>,
        perm_argument: &PermutationArgument<P::ScalarField>,
        transcript_protocol: &TranscriptProtocol<P>,
    ) -> BigQuotientPoly<P::ScalarField> {
        let gate_check_poly = self.compute_gate_check_poly(solution);

        let perm_numerator_poly = perm_argument.numerator().combined() * z;
        let perm_denominator_poly = perm_argument.denominator().combined() * z_shifted;

        let z_poly_m1 = (z - DensePolynomial::from_coefficients_slice(&[P::ScalarField::one()])) * &self.lagrange_1;

        let t = self.divide_by_vanishing(&gate_check_poly)
            + self.divide_by_vanishing(&(perm_numerator_poly - perm_denominator_poly)) * transcript_protocol.get_alpha()
            + self.divide_by_vanishing(&z_poly_m1) * transcript_protocol.get_alpha().square();

        BigQuotientPoly::create_for_domain(t, self.domain.len())
    }

    fn divide_by_vanishing(&self, poly: &DensePolynomial<P::ScalarField>) -> DensePolynomial<P::ScalarField> {
        let Zh = &DenseOrSparsePolynomial::from(&self.Zh);
        let res = DenseOrSparsePolynomial::from(poly).divide_with_q_and_r(Zh);

        let (q, r) = res.unwrap();

        if !r.is_zero() {
            panic!("Remainder is not zero");
        }

        q
    }

    fn compute_gate_check_poly(
        &self,
        solution: &PlonkSolution<P::ScalarField>,
    ) -> DensePolynomial<P::ScalarField> {
        &solution.a * &solution.b * &self.circuit.qm
            + &solution.a * &self.circuit.ql
            + &solution.b * &self.circuit.qr
            + &solution.c * &self.circuit.qo
            + &solution.public_witness.pi_combined + &self.circuit.qc
    }

    fn compute_openings(
        &self,
        solution: &PlonkSolution<P::ScalarField>,
        z_shifted: &DensePolynomial<P::ScalarField>,
        transcript_protocol: &TranscriptProtocol<P>,
    ) -> Openings<P::ScalarField> {
        let zeta = transcript_protocol.get_zeta();

        Openings {
            a: solution.a.evaluate(&zeta),
            b: solution.b.evaluate(&zeta),
            c: solution.c.evaluate(&zeta),
            s_sigma_1: self.circuit.s_sigma_1.evaluate(&zeta),
            s_sigma_2: self.circuit.s_sigma_2.evaluate(&zeta),
            z_shifted: z_shifted.evaluate(&zeta),
        }
    }

    fn linearization_poly(
        &self,
        openings: &Openings<P::ScalarField>,
        solution: &PlonkSolution<P::ScalarField>,
        permutation_argument: &PermutationArgument<P::ScalarField>,
        z_poly: &DensePolynomial<P::ScalarField>,
        big_quotient_poly: &BigQuotientPoly<P::ScalarField>,
        transcript_protocol: &TranscriptProtocol<P>,
    ) -> DensePolynomial<P::ScalarField> {
        let zeta = transcript_protocol.get_zeta();

        let gate_check_poly_linearized = &self.circuit.qm * openings.a * openings.b
            + &self.circuit.ql * openings.a
            + &self.circuit.qr * openings.b
            + &self.circuit.qo * openings.c
            + const_poly(solution.public_witness.pi_combined.evaluate(&zeta))
            + &self.circuit.qc;

        let perm_numerator_poly_linearized = z_poly * permutation_argument.numerator().evaluate(&zeta);

        let denominator = permutation_argument.denominator();
        let perm_denominator_poly_linearized = denominator.linearize_out_on_witness(&zeta)
            * denominator.evaluate_left(&zeta)
            * denominator.evaluate_right(&zeta)
            * openings.z_shifted;

        let perm_start_linearized = (z_poly - const_poly(P::ScalarField::one())) * self.lagrange_1.evaluate(&zeta);

        let linearized_vanishing_t = big_quotient_poly.linearize_on_shifts(&zeta) * self.Zh.evaluate(&zeta);

        let alpha = transcript_protocol.get_alpha();

        gate_check_poly_linearized
            + perm_numerator_poly_linearized * alpha
            - perm_denominator_poly_linearized * alpha
            + perm_start_linearized * alpha.square()
            - linearized_vanishing_t
    }
}

fn shift_poly<F: Field>(poly: &DensePolynomial<F>, scalar: &F) -> DensePolynomial<F> {
    DensePolynomial::from_coefficients_vec(
        poly.coeffs
            .iter()
            .enumerate()
            .map(|(i, c)| *c * scalar.pow([i as u64]))
            .collect(),
    )
}

#[cfg(test)]
mod tests {
    use ark_ec::{AffineRepr, CurveGroup, PrimeGroup};
    use ark_ec::pairing::Pairing;
    use ark_ff::Field;
    use ark_poly::{DenseUVPolynomial, Polynomial};
    use ark_poly::univariate::{SparsePolynomial, DensePolynomial};
    use ark_std::{One, Zero};
    use ark_test_curves::bls12_381::{Bls12_381, Fr, G1Projective};
    use crate::kzg::{BatchOpening, MultipointOpening, KZG};
    use crate::plonk::blinder::{blind_solution};
    use crate::plonk::circuit::{preprocess_circuit, CompiledCircuit};
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::plonk::circuit::PlonkSolution;
    use crate::plonk::domain::PlonkDomain;
    use crate::plonk::permutation::PermutationArgument;
    use crate::plonk::proof::{Commitments, OpeningProofs, Openings, Proof};
    use crate::plonk::protocol::{Party};
    use crate::plonk::prover::{const_poly, shift_poly, PlonkProver};
    use crate::plonk::test_utils::test_utils::{get_test_circuit, get_test_kzg, get_test_solution, hash_permutation_poly};
    use crate::plonk::transcript_protocol::TranscriptProtocol;

    struct TestEnv<'a> {
        kzg: KZG<Bls12_381>,
        domain: &'a PlonkDomain<Fr>,
        test_circuit: CompiledCircuit<Fr>,
        solution: PlonkSolution<Fr>,
        transcript: TranscriptProtocol<Bls12_381>,
        Zh: SparsePolynomial<Fr>,
        omega: Fr,
    }

    fn get_test_domain() -> PlonkDomain<Fr> {
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 4 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);

        domain
    }

    fn prepare_test_environment(domain: &PlonkDomain<Fr>) -> TestEnv {
        let test_circuit = get_test_circuit(&domain);
        let solution = get_test_solution(&domain);
        let Zh = domain.get_vanishing_polynomial();
        let kzg = get_test_kzg::<Bls12_381>(domain.len());
        let transcript = get_test_transcript::<Bls12_381>();

        TestEnv {
            omega: domain.generator(),
            domain,
            kzg,
            test_circuit,
            solution,
            transcript,
            Zh,
        }
    }

    fn get_test_transcript<P: Pairing>() -> TranscriptProtocol<P> {
        let mut transcript = TranscriptProtocol::<P>::new(P::ScalarField::one(), &[]);

        transcript.append_abc_commitments(P::G1Affine::generator(), P::G1Affine::generator(), P::G1Affine::generator());
        transcript.append_z_commitment(P::G1Affine::generator());
        transcript.append_t(P::G1Affine::generator(), P::G1Affine::generator(), P::G1Affine::generator());
        transcript.append_openings(&Openings {
            a: P::ScalarField::one(),
            b: P::ScalarField::one(),
            c: P::ScalarField::one(),
            z_shifted: P::ScalarField::one(),
            s_sigma_2: P::ScalarField::one(),
            s_sigma_1: P::ScalarField::one(),
        });
        transcript.append_opening_proofs(P::G1Affine::generator(), P::G1Affine::generator());

        transcript
    }

    #[test]
    fn compute_big_quotient_test() {
        let domain = get_test_domain();
        let TestEnv { kzg, domain, test_circuit, solution, transcript, Zh, .. } = prepare_test_environment(&domain);
        let preprocessed_circuit= preprocess_circuit(&test_circuit, &kzg);
        
        let (beta, gamma) = transcript.get_beta_gamma();
        let perm_argument = PermutationArgument::new(&domain, beta, gamma, &test_circuit, &solution);

        let z = perm_argument.z_poly();
        let z_shifted = shift_poly(&z, &domain.generator());

        let prover = PlonkProver::new(Party::new(
            &kzg,
            &domain,
            &preprocessed_circuit,
        ));

        let big_q = prover.compute_big_quotient(
            &solution,
            &z,
            &z_shifted,
            &perm_argument,
            &transcript,
        );

        let gate_check_poly = prover.compute_gate_check_poly(&solution);

        let perm_numerator_poly = perm_argument.numerator().combined() * z.clone();
        let perm_denominator_poly = perm_argument.denominator().combined() * z_shifted;
        let lagrange_base_1 = domain.lagrange_polys().first().unwrap();
        let z_poly_m1 = (z - DensePolynomial::from_coefficients_slice(&[Fr::one()])) * lagrange_base_1;
        let alpha = transcript.get_alpha();

        let restored = big_q.into_combined();

        assert_eq!(
            restored * DensePolynomial::from(Zh),
            gate_check_poly
                + (perm_numerator_poly - perm_denominator_poly) * alpha
                + z_poly_m1 * alpha.square()
        );
    }

    #[test]
    fn test_linearization_poly() {
        let domain = get_test_domain();
        let TestEnv { kzg, domain, test_circuit, solution, transcript, .. } = prepare_test_environment(&domain);
        let preprocessed_circuit= preprocess_circuit(&test_circuit, &kzg);

        let (beta, gamma) = transcript.get_beta_gamma();
        let perm_argument = PermutationArgument::new(&domain, beta, gamma, &test_circuit, &solution);

        let z = perm_argument.z_poly();
        let z_shifted = shift_poly(&z, &domain.generator());

        let prover = PlonkProver::new(Party::new(
            &kzg,
            &domain,
            &preprocessed_circuit,
        ));

        let big_q = prover.compute_big_quotient(
            &solution,
            &z,
            &z_shifted,
            &perm_argument,
            &transcript
        );

        let openings = prover.compute_openings(&solution, &z_shifted, &transcript);
        let r_poly = prover.linearization_poly(
            &openings,
            &solution,
            &perm_argument,
            &z,
            &big_q,
            &transcript,
        );

        assert_eq!(r_poly.evaluate(&transcript.get_zeta()), Fr::zero());
    }

    #[test]
    fn test_prove_verify() {
        let domain = get_test_domain();
        let TestEnv { kzg, domain, test_circuit, solution, Zh, omega, .. } = prepare_test_environment(&domain);
        let preprocessed_circuit= preprocess_circuit(&test_circuit, &kzg);
        let mut transcript = TranscriptProtocol::<Bls12_381>::new(domain.generator(), &solution.public_witness.inputs_vector);
        transcript.append_abc_commitments(kzg.commit(&solution.a).into_affine(), kzg.commit(&solution.b).into_affine(), kzg.commit(&solution.c).into_affine());
        let (beta, gamma) = transcript.get_beta_gamma();

        let solution = blind_solution(solution, &Zh);
        let perm_argument = PermutationArgument::new(&domain, beta, gamma, &test_circuit, &solution);
        let z = perm_argument.z_poly();
        transcript.append_z_commitment(kzg.commit(&z).into_affine());
        let alpha = transcript.get_alpha();
        let z_shifted = shift_poly(&z, &domain.generator());

        let lagrange_base_1 = domain.lagrange_polys().first().unwrap();

        let prover = PlonkProver::new(Party::new(
            &kzg,
            &domain,
            &preprocessed_circuit,
        ));

        let big_q = prover.compute_big_quotient(
            &solution,
            &z,
            &z_shifted,
            &perm_argument,
            &transcript,
        );
        let (lo, mid, hi) = big_q.get_splitted_polys();
        transcript.append_t(
            kzg.commit(&lo).into_affine(),
            kzg.commit(&mid).into_affine(),
            kzg.commit(&hi).into_affine(),
        );
        let zeta = transcript.get_zeta();

        let openings = prover.compute_openings(&solution, &z_shifted, &transcript);
        let r_poly = prover.linearization_poly(
            &openings,
            &solution,
            &perm_argument,
            &z,
            &big_q,
            &transcript,
        );

        transcript.append_openings(&openings);
        let vi = transcript.get_vi();

        let zeta_openings = kzg.batch_open(
            &[&r_poly, &solution.a, &solution.b, &solution.c, &test_circuit.s_sigma_1, &test_circuit.s_sigma_2],
            &zeta,
            &vi,
        );
        let zeta_omega_opening = kzg.open(
            &z_shifted,
            &zeta,
        );

        let proof = Proof::new(
            openings.clone(),
            Commitments::<Bls12_381> {
                a: kzg.commit(&solution.a),
                b: kzg.commit(&solution.b),
                c: kzg.commit(&solution.c),
                z: kzg.commit(&z),
                t_lo: kzg.commit(&lo),
                t_mid: kzg.commit(&mid),
                t_hi: kzg.commit(&hi),
            },
            OpeningProofs {
                w_zeta: zeta_openings.proof(),
                w_zeta_omega: zeta_omega_opening.proof(),
            }
        );

        transcript.append_opening_proofs(zeta_openings.proof().into_affine(), zeta_omega_opening.proof().into_affine());
        let u = transcript.get_u();

        // poly

        let mut testing_poly = &test_circuit.qm * openings.a * openings.b
            + &test_circuit.ql * openings.a
            + &test_circuit.qr * openings.b
            + &test_circuit.qo * openings.c
            + const_poly(solution.public_witness.pi_combined.evaluate(&zeta))
            + &test_circuit.qc;

        let perm_numerator_poly_linearized = z.clone() * (
            (openings.a + beta * test_circuit.sid_1.evaluate(&zeta) + gamma)
                * (openings.b + beta * test_circuit.sid_2.evaluate(&zeta) + gamma)
                * (openings.c + beta * test_circuit.sid_3.evaluate(&zeta) + gamma)
        );

        let sigma_a_linearized = openings.a + beta * openings.s_sigma_1 + gamma;
        let sigma_b_linearized = openings.b + beta * openings.s_sigma_2 + gamma;
        let sigma_c_lin_poly = const_poly(openings.c)
            + &test_circuit.s_sigma_3 * beta
            + const_poly(gamma);
        let perm_denominator_poly_linearized = sigma_c_lin_poly * sigma_a_linearized * sigma_b_linearized * openings.z_shifted;

        let perm_start_linearized = (z.clone() - const_poly(Fr::one())) * lagrange_base_1.evaluate(&zeta);

        let linearized_vanishing_t = (
            lo.clone()
                + mid.clone() * zeta.pow([domain.len() as u64])
                + hi.clone() * zeta.pow([domain.len() as u64 * 2])
        ) * Zh.evaluate(&zeta);

        testing_poly = testing_poly + perm_numerator_poly_linearized * alpha;
        testing_poly = testing_poly - perm_denominator_poly_linearized * alpha;
        testing_poly = testing_poly + perm_start_linearized * alpha.square();
        testing_poly = testing_poly - linearized_vanishing_t;

        assert_eq!(testing_poly.evaluate(&zeta), Fr::zero());
        assert_eq!(testing_poly, r_poly);

        let testing_open = kzg.open(&testing_poly, &zeta);

        let mut D = kzg.commit(&test_circuit.qm) * openings.a * openings.b
            + kzg.commit(&test_circuit.ql) * openings.a
            + kzg.commit(&test_circuit.qr) * openings.b
            + kzg.commit(&test_circuit.qo) * openings.c
            + kzg.commit(&const_poly(solution.public_witness.pi_combined.evaluate(&zeta)))
            + kzg.commit(&test_circuit.qc);

        D += proof.commitments.z * alpha * (openings.a + beta * test_circuit.sid_1.evaluate(&zeta) + gamma)
            * (openings.b + beta * test_circuit.sid_2.evaluate(&zeta) + gamma)
            * (openings.c + beta * test_circuit.sid_3.evaluate(&zeta) + gamma);

        D -= (G1Projective::generator() * openings.c + kzg.commit(&test_circuit.s_sigma_3) * beta + G1Projective::generator() * gamma) * alpha
            * (openings.a + beta * openings.s_sigma_1 + gamma)
            * (openings.b + beta * openings.s_sigma_2 + gamma) * openings.z_shifted;

        D += (proof.commitments.z - G1Projective::generator() * Fr::one()) * lagrange_base_1.evaluate(&zeta) * alpha.square();

        D -= (
            proof.commitments.t_lo
                + proof.commitments.t_mid * zeta.pow([domain.len() as u64])
                + proof.commitments.t_hi * zeta.pow([domain.len() as u64 * 2])
        ) * Zh.evaluate(&zeta);

        let mut eval = test_circuit.qm.evaluate(&zeta) * openings.a * openings.b
            + test_circuit.ql.evaluate(&zeta) * openings.a
            + test_circuit.qr.evaluate(&zeta) * openings.b
            + test_circuit.qo.evaluate(&zeta) * openings.c
            + solution.public_witness.pi_combined.evaluate(&zeta)
            + test_circuit.qc.evaluate(&zeta);

        eval += z.evaluate(&zeta) * alpha * (openings.a + beta * test_circuit.sid_1.evaluate(&zeta) + gamma)
            * (openings.b + beta * test_circuit.sid_2.evaluate(&zeta) + gamma)
            * (openings.c + beta * test_circuit.sid_3.evaluate(&zeta) + gamma);

        eval -= (openings.a + beta * openings.s_sigma_1 + gamma)
            * (openings.b + beta * openings.s_sigma_2 + gamma)
            * (openings.c + beta * test_circuit.s_sigma_3.evaluate(&zeta) + gamma) * alpha * openings.z_shifted;

        eval += alpha.square() * (z.evaluate(&zeta) - Fr::one()) * lagrange_base_1.evaluate(&zeta);

        eval -= Zh.evaluate(&zeta)
            * (
            lo.evaluate(&zeta)
                + mid.evaluate(&zeta) * zeta.pow([domain.len() as u64])
                + hi.evaluate(&zeta) * zeta.pow([domain.len() as u64 * 2])
        );

        assert_eq!(testing_open.evaluation, eval);

        let is_valid = kzg.check_multipoint(
            &[
                MultipointOpening {
                    batch_opening: &BatchOpening::new(
                        vec![Fr::zero()],
                        testing_open.proof(),
                    ),
                    linearization_scalar: &vi,
                    opening_point: &zeta,
                    commitments: &[D],
                }
            ],
            &u,
        );

        assert!(is_valid);

        assert_eq!(z_shifted.evaluate(&zeta), z.evaluate(&(zeta * omega)));
    }

    #[test]
    fn test_linearized_parts() {
        let domain = get_test_domain();
        let TestEnv { kzg, domain, test_circuit, solution, Zh, transcript, .. } = prepare_test_environment(&domain);
        let preprocessed_circuit= preprocess_circuit(&test_circuit, &kzg);
        let solution = blind_solution(solution, &Zh);
        let (beta, gamma) = transcript.get_beta_gamma();
        let perm_argument = PermutationArgument::new(&domain, beta, gamma, &test_circuit, &solution);
        let z = perm_argument.z_poly();
        let z_shifted = shift_poly(&z, &domain.generator());

        let prover = PlonkProver::new(Party::new(
            &kzg,
            &domain,
            &preprocessed_circuit,
        ));

        let openings = prover.compute_openings(&solution, &z_shifted, &transcript);

        let zeta = transcript.get_zeta();

        let perm_denominator_poly = perm_argument.denominator().combined() * z_shifted;

        let sigma_a_linearized = openings.a + beta * openings.s_sigma_1 + gamma;
        let sigma_b_linearized = openings.b + beta * openings.s_sigma_2 + gamma;
        let sigma_c_lin_poly = const_poly(openings.c)
            + &test_circuit.s_sigma_3 * beta
            + const_poly(gamma);
        let perm_denominator_poly_linearized = sigma_c_lin_poly * sigma_a_linearized * sigma_b_linearized * openings.z_shifted;

        assert_eq!(openings.a, solution.a.evaluate(&zeta));
        assert_eq!(openings.s_sigma_1, test_circuit.s_sigma_1.evaluate(&zeta));

        assert_eq!(sigma_a_linearized, hash_permutation_poly(&solution.a, &test_circuit.s_sigma_1, beta, gamma).evaluate(&zeta));

        assert_eq!(perm_denominator_poly_linearized.evaluate(&zeta), perm_denominator_poly.evaluate(&zeta));
    }
}