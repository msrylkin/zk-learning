use ark_ec::pairing::Pairing;
use crate::plonk::proof::Openings;
use crate::transcript::Transcript;

pub struct TranscriptProtocol<P: Pairing> {
    beta: Option<P::ScalarField>,
    gamma: Option<P::ScalarField>,
    alpha: Option<P::ScalarField>,
    zeta: Option<P::ScalarField>,
    vi: Option<P::ScalarField>,
    u: Option<P::ScalarField>,
    round: usize,
    transcript: Transcript<P::ScalarField>,
}

impl<P: Pairing> TranscriptProtocol<P> {
    pub fn new(omega: P::ScalarField, public_inputs: &[P::ScalarField]) -> Self {
        let mut transcript = Transcript::new(b"plonk");
        transcript.append_f(omega);

        for pi in public_inputs {
            transcript.append_f(*pi);
        }

        Self {
            beta: None,
            gamma: None,
            alpha: None,
            zeta: None,
            vi: None,
            u: None,
            transcript,
            round: 0,
        }
    }

    pub fn get_beta_gamma(&self) -> (P::ScalarField, P::ScalarField) {
        (self.beta.unwrap(), self.gamma.unwrap())
    }

    pub fn get_alpha(&self) -> P::ScalarField {
        self.alpha.unwrap()
    }

    pub fn get_zeta(&self) -> P::ScalarField {
        self.zeta.unwrap()
    }

    pub fn get_vi(&self) -> P::ScalarField {
        self.vi.unwrap()
    }

    pub fn get_u(&self) -> P::ScalarField {
        self.u.unwrap()
    }

    pub fn append_abc_commitments(&mut self, a: P::G1Affine, b: P::G1Affine, c: P::G1Affine) {
        assert_eq!(self.round, 0);

        self.transcript.append_g1(a);
        self.transcript.append_g1(b);
        self.transcript.append_g1(c);

        let beta = self.transcript.derive();
        self.transcript.append_f(beta);
        let gamma = self.transcript.derive();
        self.transcript.append_f(gamma);

        self.beta = Some(beta);
        self.gamma = Some(gamma);
        self.round = 1;
    }

    pub fn append_z_commitment(&mut self, z: P::G1Affine) {
        assert_eq!(self.round, 1);

        self.transcript.append_g1(z);

        let alpha = self.transcript.derive();
        self.transcript.append_f(alpha);

        self.alpha = Some(alpha);
        self.round = 2;
    }

    pub fn append_t(&mut self, t_lo: P::G1Affine, t_mid: P::G1Affine, t_hi: P::G1Affine) {
        assert_eq!(self.round, 2);

        self.transcript.append_g1(t_lo);
        self.transcript.append_g1(t_mid);
        self.transcript.append_g1(t_hi);

        let zeta = self.transcript.derive();
        self.transcript.append_f(zeta);

        self.zeta = Some(zeta);
        self.round = 3;
    }

    pub fn append_openings(
        &mut self,
        openings: &Openings<P::ScalarField>
    ) {
        assert_eq!(self.round, 3);

        self.transcript.append_f(openings.a);
        self.transcript.append_f(openings.b);
        self.transcript.append_f(openings.c);
        self.transcript.append_f(openings.s_sigma_1);
        self.transcript.append_f(openings.s_sigma_2);
        self.transcript.append_f(openings.z_shifted);

        let vi = self.transcript.derive();
        self.transcript.append_f(vi);

        self.vi = Some(vi);
        self.round = 4;
    }

    pub fn append_opening_proofs(&mut self, zeta_opening: P::G1Affine, zeta_omega_opening: P::G1Affine) {
        assert_eq!(self.round, 4);

        self.transcript.append_g1(zeta_opening);
        self.transcript.append_g1(zeta_omega_opening);

        let u = self.transcript.derive();
        self.transcript.append_f(u);

        self.u = Some(u);
        self.round = 5;
    }
}