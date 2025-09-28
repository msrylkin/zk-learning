#![allow(non_snake_case)]
#![allow(dead_code)]
#![warn(missing_docs)]

//! # Educational ZKP Crate
//!
//! This crate provides implementations of zero-knowledge proof protocols for educational purposes.
//!
//! ## GKR
//! Proves satisfiability of layered arithmetic circuits using the SumCheck protocol for multivariate polynomials.
//!
//! ## PLONK
//! zk-SNARK protocol for arbitrary circuits via a single polynomial `T(X)` checked against a vanishing polynomial.

mod kzg;
mod schnorr;
mod sumcheck;
mod poly_utils;
mod gkr;
mod random_oracle;
mod plonk;
mod evaluation_domain;
mod transcript;

pub use gkr::*;
pub use sumcheck::*;
pub use kzg::*;
pub use plonk::*;