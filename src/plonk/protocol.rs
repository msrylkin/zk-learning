use ark_ec::pairing::Pairing;
use crate::kzg::KZG;

struct PlonkProtocol<P: Pairing> {
    kzg: KZG<P>,
}