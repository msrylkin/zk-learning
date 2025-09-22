use std::ops::Deref;
use ark_ff::{Field, PrimeField};
use crate::evaluation_domain::MultiplicativeSubgroup;

pub struct PlonkDomain<F: PrimeField> {
    domain: MultiplicativeSubgroup<F>,
    k1: F,
    k2: F,
}

impl<F: PrimeField> Deref for PlonkDomain<F> {
    type Target = MultiplicativeSubgroup<F>;

    fn deref(&self) -> &Self::Target {
        &self.domain
    }
}

impl<'a, F: PrimeField> IntoIterator for &'a PlonkDomain<F> {
    type Item = &'a F;

    type IntoIter = std::slice::Iter<'a, F>;

    fn into_iter(self) -> Self::IntoIter {
        self.domain.into_iter()
    }
}

impl<F: PrimeField> PlonkDomain<F> {
    fn new(domain: MultiplicativeSubgroup<F>, k1: F, k2: F) -> Self {
        Self { domain, k1, k2 }
    }

    pub fn k1(&self) -> F {
        self.k1
    }

    pub fn k2(&self) -> F {
        self.k2
    }

    pub fn create_from_subgroup(domain: MultiplicativeSubgroup<F>) -> Self {
        let (k1, k2) = Self::pick_coset_shifters(&domain);

        Self::new(domain, k1, k2)
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

#[cfg(test)]
mod tests {
    use ark_test_curves::bls12_381::Fr;
    use crate::evaluation_domain::generate_multiplicative_subgroup;
    use crate::plonk::domain::PlonkDomain;

    #[test]
    fn pick_coset_shifters_test() {
        let domain = generate_multiplicative_subgroup::<{ 1u64 << 6 }, Fr>();
        let domain = PlonkDomain::create_from_subgroup(domain);

        let coset_k1 = domain.iter().map(|e| domain.k1() * e).collect::<Vec<_>>();
        let coset_k2 = domain.iter().map(|e| domain.k2() * e).collect::<Vec<_>>();

        for h in &domain {
            for k1h in &coset_k1 {
                for k2h in &coset_k2 {
                    assert_ne!(h, k1h);
                    assert_ne!(k1h, k2h);
                }
            }
        }
    }
}