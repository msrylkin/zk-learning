use std::hash::{DefaultHasher, Hash, Hasher};
use std::marker::PhantomData;
use ark_ec::{AffineRepr};
use ark_ff::Field;
use ark_ff::field_hashers::{DefaultFieldHasher, HashToField};
use sha2::Sha256;

pub struct Transcript<F: Field> {
    current: Vec<u8>,
    hasher: DefaultFieldHasher<Sha256>,
    phantom: PhantomData<F>,
}

impl<F: Field> Transcript<F> {
    pub fn new(init: &[u8]) -> Self {
        let hasher = <DefaultFieldHasher<Sha256> as HashToField<F>>::new(init);

        Self {
            current: init.to_vec(),
            hasher,
            phantom: PhantomData,
        }
    }

    pub fn append_g1<P: AffineRepr>(&mut self, element: P) {
        let mut hasher = DefaultHasher::new();

        element.x().hash(&mut hasher);
        element.y().hash(&mut hasher);

        let value = hasher.finish();
        self.current.extend(value.to_be_bytes());
    }

    pub fn append_f(&mut self, element: F) {
        let mut hasher = DefaultHasher::new();

        element.hash(&mut hasher);

        let value = hasher.finish();
        self.current.extend(value.to_be_bytes());
    }

    pub fn derive(&self) -> F {
        let [e]: [F; 1] = self.hasher.hash_to_field(&self.current);

        e
    }
}

#[cfg(test)]
mod tests {
    use ark_ec::AffineRepr;
    use ark_test_curves::bls12_381::{Fr, G1Affine};
    use crate::transcript::Transcript;

    #[test]
    pub fn test_transcript() {
        let mut transcript = Transcript::<Fr>::new(&[7,7,7]);

        assert_eq!(transcript.current, vec![7,7,7]);

        let h_1 = transcript.derive();

        transcript.append_f(Fr::from(999));
        let h_2 = transcript.derive();

        assert_ne!(h_1, h_2);

        transcript.append_g1(G1Affine::generator());
        let h_3 = transcript.derive();

        assert_ne!(h_2, h_3);
        assert_ne!(h_1, h_3);
    }
}