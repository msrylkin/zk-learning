use std::cell::RefCell;
use ark_ff::Field;

pub trait RandomOracle {
    type Item: Field;

    fn get_randomness(&self, n: usize) -> Vec<Self::Item>;
}

pub struct FixedRandomOracle<F: Field> {
    current_index: RefCell<usize>,
    values: Vec<F>
}

impl<F: Field> FixedRandomOracle<F> {
    pub fn new(values: Vec<F>) -> Self {
        assert_ne!(values.len(), 0);

        Self {
            current_index: RefCell::new(0),
            values
        }
    }
}

impl<F: Field> RandomOracle for FixedRandomOracle<F> {
    type Item = F;

    fn get_randomness(&self, n: usize) -> Vec<F> {
        let mut remaining = n;
        let mut randomness = vec![];
        let mut current_index = self.current_index.borrow_mut();

        while remaining > 0 {
            randomness.push(self.values[*current_index]);
            remaining -= 1;
            *current_index += 1;
            
            if *current_index == self.values.len() {
                *current_index = 0;
            }
        }
        
        randomness
    }
}

struct ProxyRandomOracle<'a, P> {
    parent: &'a P,
}

impl<'a, F: Field, P: RandomOracle<Item = F>> ProxyRandomOracle<'a, P> {
    pub fn new(parent_oracle: &'a P) -> Self {
        Self {
            parent: parent_oracle,
        }
    }
}

impl<F: Field, P: RandomOracle<Item = F>> RandomOracle for ProxyRandomOracle<'_, P> {
    type Item = F;

    fn get_randomness(&self, n: usize) -> Vec<F> {
        self.parent.get_randomness(n)
    }
}

#[cfg(test)]
mod test {
    use ark_test_curves::bls12_381::Fr;
    use crate::poly_utils::to_f;
    use crate::random_oracle::{FixedRandomOracle, RandomOracle};

    #[test]
    fn test_fixed_random_oracle() {
        let fixed_random_oracle = FixedRandomOracle::new(to_f::<Fr>(vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 100]));
        
        let n1 = fixed_random_oracle.get_randomness(1);
        assert_eq!(n1, vec![Fr::from(1)]);
        
        let n5 = fixed_random_oracle.get_randomness(5);
        assert_eq!(n5, to_f(vec![2,3,4,5,6]));
        
        let n12 = fixed_random_oracle.get_randomness(12);
        assert_eq!(n12, to_f(vec![7, 8, 9, 100, 1, 2, 3, 4, 5, 6, 7 ,8]));
    }
}