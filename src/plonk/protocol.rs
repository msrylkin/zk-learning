use ark_ec::pairing::Pairing;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use crate::kzg::KZG;
use crate::plonk::domain::PlonkDomain;

struct PlonkProtocol<P: Pairing> {
    kzg: KZG<P>,
    domain: PlonkDomain<P::ScalarField>,
    Zh: SparsePolynomial<P::ScalarField>,
    lagrange_1: DensePolynomial<P::ScalarField>,
}

impl<P: Pairing> PlonkProtocol<P> {
    pub fn new(
        kzg: KZG<P>,
        domain: PlonkDomain<P::ScalarField>,
    ) -> Self {
        Self {
            kzg,
            Zh: domain.get_vanishing_polynomial(),
            lagrange_1: domain.lagrange_polys().first().unwrap().clone(),
            domain,
        }
    }

    pub fn prove(

    ) {

    }
}