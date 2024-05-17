use p3_air::{Air, AirBuilder, BaseAir};
use p3_baby_bear::{BabyBear, DiffusionMatrixBabyBear};
use p3_challenger::DuplexChallenger;
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::{Field, PrimeField32, PrimeField64};
use p3_field::extension::BinomialExtensionField;
use p3_fri::TwoAdicFriPcs;
use p3_matrix::dense::{DenseMatrix, RowMajorMatrix};
use p3_matrix::Matrix;
use p3_merkle_tree::FieldMerkleTreeMmcs;
use p3_poseidon2::{Poseidon2, Poseidon2ExternalMatrixGeneral};
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_uni_stark::StarkConfig;
use std::borrow::Borrow;

const PLONK_GATE_WIDTH: usize = 8;

struct PlonkRow<F> {
    pub q_l: F,
    pub q_r: F,
    pub q_o: F,
    pub q_m: F,
    pub q_c: F,
    pub a: F,
    pub b: F,
    pub c: F,
}

impl<F> Borrow<PlonkRow<F>> for [F] {
    fn borrow(&self) -> &PlonkRow<F> {
        debug_assert_eq!(self.len(), PLONK_GATE_WIDTH);
        let (prefix, shorts, suffix) = unsafe { self.align_to::<PlonkRow<F>>() };
        debug_assert!(prefix.is_empty(), "Alignment should match");
        debug_assert!(suffix.is_empty(), "Alignment should match");
        debug_assert_eq!(shorts.len(), 1);
        &shorts[0]
    }
}

impl<F: PrimeField64> Default for PlonkRow<F> {
    fn default() -> Self {
        PlonkRow {
            q_l: F::zero(),
            q_r: F::zero(),
            q_o: F::zero(),
            q_m: F::zero(),
            q_c: F::zero(),
            a: F::zero(),
            b: F::zero(),
            c: F::zero(),
        }
    }
}

impl<F> BaseAir<F> for PlonkAir {
    fn width(&self) -> usize {
        PLONK_GATE_WIDTH
    }
}

struct PlonkAir {}

impl<AB: AirBuilder> Air<AB> for PlonkAir {
    fn eval(&self, builder: &mut AB) {
        let main = builder.main();
        let local = main.row_slice(0);
        let local: &PlonkRow<AB::Var> = (*local).borrow();
        let mut when_transition = builder.when_transition();

        when_transition.assert_zero((local.q_l * local.a) + (local.q_r * local.b)
            + (local.q_o * local.c) + (local.q_m * (local.a * local.b)) + (local.q_c * local.c));
    }
}

fn generate_add_gate<F: PrimeField64>(a: F, b: F, c: F) -> PlonkRow<F> {
    PlonkRow {
        a,
        b,
        c,
        q_l: F::one(),
        q_r: F::one(),
        q_o: F::neg_one(),
        q_m: F::zero(),
        q_c: F::zero(),
    }
}

fn generate_mul_gate<F: PrimeField64>(a: F, b: F, c: F) -> PlonkRow<F> {
    PlonkRow {
        a,
        b,
        c,
        q_l: F::zero(),
        q_r: F::zero(),
        q_o: F::neg_one(),
        q_m: F::one(),
        q_c: F::zero(),
    }
}

fn generate_trace<F: PrimeField64>(n: usize) -> RowMajorMatrix<F> {
    let mut trace = RowMajorMatrix::new(vec![F::zero(); n * PLONK_GATE_WIDTH], PLONK_GATE_WIDTH);
    let (prefix, rows, suffix) = unsafe {
        trace.values.align_to_mut::<PlonkRow<F>>()
    };

    assert!(prefix.is_empty(), "Alignment check! Ethereum aligned??");
    assert!(suffix.is_empty(), "Alignment check! Ethereum aligned??");

    rows[0] = generate_add_gate(F::one(), F::one(), F::one() + F::one());
    rows[1] = generate_mul_gate(F::from_canonical_u32(2), F::one(), F::from_canonical_u32(2));

    for i in 2..n {
        rows[i] = PlonkRow::default();
    }

    trace
}

type Val = BabyBear;
type Perm = Poseidon2<Val, Poseidon2ExternalMatrixGeneral, DiffusionMatrixBabyBear, 16, 7>;
type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
type MyCompress = TruncatedPermutation<Perm, 2, 8, 16>;
type ValMmcs =
FieldMerkleTreeMmcs<<Val as Field>::Packing, <Val as Field>::Packing, MyHash, MyCompress, 8>;
type Challenge = BinomialExtensionField<Val, 4>;
type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
type Challenger = DuplexChallenger<Val, Perm, 16, 8>;
type Dft = Radix2DitParallel;
type Pcs = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;
type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;

#[cfg(test)]
mod test {
    use p3_baby_bear::{BabyBear, DiffusionMatrixBabyBear};
    use p3_field::AbstractField;
    use p3_fri::FriConfig;
    use p3_keccak_air::generate_trace_rows;
    use p3_poseidon2::Poseidon2ExternalMatrixGeneral;
    use p3_uni_stark::{prove, verify};
    use p3_util::array_serialization::serialize;
    use rand::thread_rng;
    use serde::Serialize;
    use crate::{ChallengeMmcs, Challenger, Dft, generate_trace, MyCompress, MyConfig, MyHash, Pcs, Perm, PlonkAir, Val, ValMmcs};

    #[test]
    fn test_sdf() {
        let perm = Perm::new_from_rng_128(
            Poseidon2ExternalMatrixGeneral,
            DiffusionMatrixBabyBear,
            &mut thread_rng(),
        );
        let hash = MyHash::new(perm.clone());
        let compress = MyCompress::new(perm.clone());
        let val_mmcs = ValMmcs::new(hash, compress);
        let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());
        let dft = Dft {};
        let trace = generate_trace(4);
        let fri_config = FriConfig {
            log_blowup: 2,
            num_queries: 28,
            proof_of_work_bits: 8,
            mmcs: challenge_mmcs,
        };
        let pcs = Pcs::new(dft, val_mmcs, fri_config);
        let config = MyConfig::new(pcs);
        let mut challenger = Challenger::new(perm.clone());
        let proof = prove(&config, &PlonkAir {}, &mut challenger, trace, &vec![]);

        let mut challenger = Challenger::new(perm);
        verify(&config, &PlonkAir {}, &mut challenger, &proof, &vec![]).expect("verification failed");
    }
}