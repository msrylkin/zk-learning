# **Educational stuff on ZKP**
___
## GKR 

A protocol for proving circuit satisfiability for layered arithmetic circuits.  
It is based on the SumCheck protocol for multivariate polynomials.

Example usage:

```RUST
use ark_test_curves::bls12_381::Fr;
use crate::gkr::{CircuitBuilder, GKRProtocol, Layer};
use crate::random_oracle::FixedRandomOracle;

fn main() {
    let circuit = CircuitBuilder::new(vec![10, 200, 20, 300].into_iter().map(Fr::from).collect())
        .add_layer(
            Layer::new()
                .add_addition_gate((0, 1))
                .add_addition_gate((2, 3))
        )
        .add_layer(
            Layer::new()
                .add_multiplication_gate((0, 1))
        )
        .build();

    let solution = circuit.solve();

    let random_oracle = FixedRandomOracle::new(vec![1,2,3,4,5,6,7,8,9,10].into_iter().map(Fr::from).collect());
    let gkr_protocol = GKRProtocol::new(&random_oracle);

    let gkr_proof = gkr_protocol.prove(&circuit, &solution);

    gkr_protocol.verify(&circuit, &gkr_proof);
}
```

## PLONK

A zk-SNARK protocol for proving circuit satisfiability for arbitrary circuits.  
The protocol constructs a single polynomial `T(X)` that encodes all circuit checks and proves that

`Zh(X) ∣ T(X)`

where `Zh(X)` is the vanishing polynomial over a chosen multiplicative subgroup of the field.

```rust
use ark_test_curves::bls12_381::{Bls12_381, Fr};
use crate::evaluation_domain::generate_multiplicative_subgroup;
use crate::kzg::{setup, KZG};
use crate::plonk::{CircuitDescription, PlonkDomain, PlonkProtocol};

fn get_circuit() -> CircuitDescription<Fr> {
    let mut circuit = CircuitDescription::new();

    let a = circuit.add_variable();
    let b = circuit.add_variable();
    let c = circuit.constant_var(Fr::from(82));
    circuit.make_public(a);

    let mul_result_1 = circuit.multiplication_gate(a, b);
    let mul_result_2 = circuit.multiplication_gate(mul_result_1, c);

    circuit.make_output(mul_result_2);

    circuit
}

fn setup_kzg(n: usize) -> KZG<Bls12_381> {
    let tau = Fr::from(777);
    let config = setup::<Bls12_381>(n * 2, tau);
    let kzg = KZG::new(config);

    kzg
}


fn main() {
    let domain = PlonkDomain::create_from_subgroup(generate_multiplicative_subgroup::<{ 1 << 4 }, Fr>());
    let circuit = get_circuit();
    let public = vec![Fr::from(90)];
    let private = vec![Fr::from(893)];
    let solution = circuit.solve(&public, &private, &domain);

    let kzg = setup_kzg(domain.len());

    let public_input = solution.public_witness.clone();

    let compiled_circuit = circuit.compile(&domain);

    let protocol = PlonkProtocol::new(kzg, domain);
    let instance = protocol.create_instance(&compiled_circuit);

    let proof = instance.prove(solution);
    instance.verify(&proof, &public_input);
}
```