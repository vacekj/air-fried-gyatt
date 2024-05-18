
use p3_air::{Air, AirBuilder, BaseAir};





use p3_field::AbstractField;
use p3_field::{Field, PrimeField64};

use p3_matrix::dense::{RowMajorMatrix};
use p3_matrix::Matrix;




use std::borrow::Borrow;

const PLONK_GATE_WIDTH: usize = 15;

mod circuit_checker;

struct PlonkRow<F> {
    pub q_l: F,
    pub q_r: F,
    pub q_o: F,
    pub q_m: F,
    pub q_c: F,
    pub a: F,
    pub b: F,
    pub c: F,
    pub copy: CopyConstraints<F>,
}

struct CopyConstraints<F> {
    pub id_0: F,
    pub id_1: F,
    pub id_2: F,
    pub sigma_0: F,
    pub sigma_1: F,
    pub sigma_2: F,
    pub acc: F,
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

impl<F: PrimeField64> Default for CopyConstraints<F> {
    fn default() -> Self {
        CopyConstraints {
            id_0: F::zero(),
            id_1: F::zero(),
            id_2: F::zero(),
            sigma_0: F::zero(),
            sigma_1: F::zero(),
            sigma_2: F::zero(),
            acc: F::zero(),
        }
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
            copy: CopyConstraints::default(),
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
        let (row, shift) = (main.row_slice(0), main.row_slice(1));
        let row: &PlonkRow<AB::Var> = (*row).borrow();
        let shift: &PlonkRow<AB::Var> = (*shift).borrow();
        let mut when_transition = builder.when_transition();

        when_transition.assert_zero(
            (row.q_l * row.a)
                + (row.q_r * row.b)
                + (row.q_o * row.c)
                + (row.q_m * (row.a * row.b))
                + (row.q_c * row.c),
        );

        build_copy_constraints(builder, row, shift);
    }
}

fn build_copy_constraints<AB: AirBuilder>(
    builder: &mut AB,
    row: &PlonkRow<AB::Var>,
    shift: &PlonkRow<AB::Var>,
) {
    // TOOD: require randomness
    let one = AB::Expr::one();
    let numerator = (row.a + row.copy.id_0 + one.clone())
        * (row.b + row.copy.id_1 + one.clone())
        * (row.c + row.copy.id_2 + one.clone());
    let denominator = (row.a + row.copy.sigma_0 + one.clone())
        * (row.b + row.copy.sigma_1 + one.clone())
        * (row.c + row.copy.sigma_2 + one.clone());

    builder.assert_zero((row.copy.acc * numerator) - (shift.copy.acc * denominator));
}

// We write in empty copy constraints, as we will fill them in later
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
        copy: CopyConstraints::default(),
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
        copy: CopyConstraints::default(),
    }
}

fn generate_trace<F: PrimeField64>(n: usize) -> RowMajorMatrix<F> {
    let mut trace = RowMajorMatrix::new(vec![F::zero(); n * PLONK_GATE_WIDTH], PLONK_GATE_WIDTH);
    let (prefix, rows, suffix) = unsafe { trace.values.align_to_mut::<PlonkRow<F>>() };

    assert!(prefix.is_empty(), "Alignment check! Ethereum aligned??");
    assert!(suffix.is_empty(), "Alignment check! Ethereum aligned??");

    rows[0] = generate_add_gate(F::one(), F::one(), F::one() + F::one());
    rows[1] = generate_mul_gate(F::from_canonical_u32(2), F::one(), F::from_canonical_u32(2));

    for i in 2..n {
        rows[i] = PlonkRow::default();
    }

    // Generate sub groups for copy constraints
    let n_f = F::from_canonical_usize(n);
    let two = F::from_canonical_u16(2);
    let three = F::from_canonical_usize(3);

    for (i, row)  in rows.iter_mut().enumerate() {
        let i_f = F::from_canonical_usize(i);
        row.copy.id_0 = i_f;
        row.copy.id_1 = n_f.mul(two) + i_f;
        row.copy.id_2 = (n_f.mul(three)) + i_f;

        // Same sigmas for now
        row.copy.sigma_0 = i_f;
        row.copy.sigma_1 = n_f.mul(two) + i_f;
        row.copy.sigma_2 = (n_f.mul(three)) + i_f;
    }

    // Calculate grand product polynomial
    // Probably just resize a vector unsafe?
    let mut numerator = [F::zero(); 4];
    let mut denominator = [F::zero(); 4];

    for i in 0..n {
        numerator[i] = calculate_numerator(&rows[i]);
        // TODO: batch inverse
        denominator[i] = calculate_denominator(&rows[i]);
    }

    for i in 0..n - 1 {
        // Calculate running numerator and denominator products
        numerator[i + 1] *= numerator[i];
        denominator[i + 1] *= denominator[i];
    }

    // Calculate grand product accumulator
    rows[0].copy.acc = F::one();
    for i in 1..n {
        rows[i].copy.acc = numerator[i] * denominator[i];
    }

    trace
}

fn calculate_numerator<F: PrimeField64>(row: &PlonkRow<F>) -> F {
    (row.a + row.copy.id_0 + F::one())
        * (row.b + row.copy.id_1 + F::one())
        * (row.c + row.copy.id_2 + F::one())
}

fn calculate_denominator<F: PrimeField64>(row: &PlonkRow<F>) -> F {
    ((row.a + row.copy.sigma_0 + F::one())
        * (row.b + row.copy.sigma_1 + F::one())
        * (row.c + row.copy.sigma_2 + F::one()))
    .inverse()
}

#[cfg(test)]
mod test {
    use crate::{circuit_checker::check_constraints, generate_trace, PlonkAir};
    use p3_baby_bear::{BabyBear, DiffusionMatrixBabyBear};
    use p3_challenger::DuplexChallenger;
    use p3_commit::ExtensionMmcs;
    use p3_dft::Radix2DitParallel;
    use p3_field::{extension::BinomialExtensionField, Field};
    use p3_fri::{FriConfig, TwoAdicFriPcs};
    use p3_matrix::dense::DenseMatrix;
    use p3_merkle_tree::FieldMerkleTreeMmcs;
    use p3_poseidon2::{Poseidon2, Poseidon2ExternalMatrixGeneral};
    use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
    use p3_uni_stark::{prove, verify, StarkConfig};
    use rand::thread_rng;

    type Val = BabyBear;
    type Perm = Poseidon2<Val, Poseidon2ExternalMatrixGeneral, DiffusionMatrixBabyBear, 16, 7>;
    type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
    type MyCompress = TruncatedPermutation<Perm, 2, 8, 16>;
    type ValMmcs = FieldMerkleTreeMmcs<
        <Val as Field>::Packing,
        <Val as Field>::Packing,
        MyHash,
        MyCompress,
        8,
    >;
    type Challenge = BinomialExtensionField<Val, 4>;
    type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
    type Challenger = DuplexChallenger<Val, Perm, 16, 8>;
    type Dft = Radix2DitParallel;
    type Pcs = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;
    type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;

    fn prove_and_verify(trace: &DenseMatrix<BabyBear>) {
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
        let fri_config = FriConfig {
            log_blowup: 2,
            num_queries: 28,
            proof_of_work_bits: 8,
            mmcs: challenge_mmcs,
        };
        let pcs = Pcs::new(dft, val_mmcs, fri_config);
        let config = MyConfig::new(pcs);
        let mut challenger = Challenger::new(perm.clone());

        check_constraints(&PlonkAir {}, trace, &vec![]);

        let proof = prove(
            &config,
            &PlonkAir {},
            &mut challenger,
            trace.clone(),
            &vec![],
        );

        let mut challenger = Challenger::new(perm);
        verify(&config, &PlonkAir {}, &mut challenger, &proof, &vec![])
            .expect("verification failed");
    }

    #[test]
    fn test_simple_add_gates() {
        let trace = generate_trace(4);
        prove_and_verify(&trace);
    }
}
