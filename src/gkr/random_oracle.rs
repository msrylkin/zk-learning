use std::cell::RefCell;
use ark_ff::Field;

pub trait RandomOracle<F> {
    fn get_randomness(&self, n: usize) -> Vec<F>;
}

pub struct FixedRandomOracle<F: Field> {
    current_index: RefCell<usize>,
    values: Vec<F>
}

impl <F: Field> FixedRandomOracle<F> {
    pub fn new(values: Vec<F>) -> Self {
        assert_ne!(values.len(), 0);

        Self {
            current_index: RefCell::new(0),
            values
        }
    }
}

impl <F: Field> RandomOracle<F> for FixedRandomOracle<F> {
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