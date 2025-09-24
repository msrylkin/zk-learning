use ark_ec::pairing::Pairing;
use ark_poly::univariate::{DensePolynomial, SparsePolynomial};
use crate::kzg::KZG;
use crate::plonk::circuit::CompiledCircuit;
use crate::plonk::domain::PlonkDomain;

pub struct Party<'a, P: Pairing> {
    pub kzg: &'a KZG<P>,
    pub domain: &'a PlonkDomain<P::ScalarField>,
    pub Zh: SparsePolynomial<P::ScalarField>,
    pub circuit: &'a CompiledCircuit<'a, P::ScalarField>,
    pub lagrange_1: DensePolynomial<P::ScalarField>,
}

impl<'a, P: Pairing> Party<'a, P> {
    pub fn new(
        kzg: &'a KZG<P>,
        domain: &'a PlonkDomain<P::ScalarField>,
        circuit: &'a CompiledCircuit<'a, P::ScalarField>,
    ) -> Self {
        Self {
            kzg,
            Zh: domain.get_vanishing_polynomial(),
            lagrange_1: domain.lagrange_polys().first().unwrap().clone(),
            domain,
            circuit,
        }
    }
}

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