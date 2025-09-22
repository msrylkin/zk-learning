use std::ops::Deref;
use ark_ff::{Field, PrimeField};
use crate::evaluation_domain::MultiplicativeSubgroup;

pub struct PlonkDomain<'a, F: PrimeField> {
    domain: &'a MultiplicativeSubgroup<F>,
    k1: F,
    k2: F,
}

impl<F: PrimeField> Deref for PlonkDomain<'_, F> {
    type Target = MultiplicativeSubgroup<F>;

    fn deref(&self) -> &Self::Target {
        &self.domain
    }
}

impl<'a, F: PrimeField> IntoIterator for &'a PlonkDomain<'_, F> {
    type Item = &'a F;
    
    type IntoIter = std::slice::Iter<'a, F>;
    
    fn into_iter(self) -> Self::IntoIter {
        self.domain.into_iter()
    }
}

impl<'a, F: PrimeField> PlonkDomain<'a, F> {
    pub fn new(domain: &'a MultiplicativeSubgroup<F>) -> Self {
        let (k1, k2) = Self::pick_coset_shifters(domain);
        
        Self { domain, k1, k2 }
    }
    
    pub fn k1(&self) -> F {
        self.k1
    }

    pub fn k2(&self) -> F {
        self.k2
    }

    fn pick_coset_shifters(domain: &[F]) -> (F, F) {
        let mut i = 2;
        let n = domain.len();
        let (k1, k2) = loop {
            let k1 = F::from(i);
            let k2 = k1.square();

            let k1_n = k1.pow(&[n as u64]);
            let k2_n  = k2.pow(&[n as u64]);
            let k1_over_k2 = (k1 / k2).pow(&[n as u64]);

            if !k1_n.is_one() && !k2_n.is_one() && !k1_over_k2.is_one() {
                break (k1, k2);
            }

            i += 1;
        };

        (k1, k2)
    }
}