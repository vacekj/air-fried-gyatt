use std::borrow::Borrow;
use std::fmt;
use std::fmt::{Display, Formatter};

use acir::circuit::Opcode;
use acir::native_types::WitnessMap;
use acir::FieldElement;
use p3_air::{Air, AirBuilder, BaseAir};
pub use p3_baby_bear::BabyBear;
use p3_baby_bear::DiffusionMatrixBabyBear;
use p3_challenger::DuplexChallenger;
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::{extension::BinomialExtensionField, AbstractField, Field, PrimeField64};
use p3_fri::{FriConfig, TwoAdicFriPcs};
use p3_matrix::dense::DenseMatrix;
use p3_matrix::dense::RowMajorMatrix;
use p3_matrix::Matrix;
use p3_merkle_tree::FieldMerkleTreeMmcs;
use p3_poseidon2::{Poseidon2, Poseidon2ExternalMatrixGeneral};
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_uni_stark::{prove, verify, StarkConfig};
use rand::thread_rng;

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

use crate::circuit_checker::check_constraints;

const PLONK_GATE_WIDTH: usize = 15;

mod circuit_checker;

#[derive(Clone)]
pub struct PlonkRow<F> {
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

impl<F: Display> Display for PlonkRow<F> {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} * ({} * {}) + {} * {} + {} * {} + {} * {} + {} = 0",
            self.q_m,
            self.a,
            self.b,
            self.q_l,
            self.a,
            self.q_r,
            self.b,
            self.q_o,
            self.c,
            self.q_c,
        )
    }
}

impl<F: PrimeField64> PlonkRow<F> {
    fn set_linear_term(&mut self, x: F, witness: F) {
        println!(
            "x: {}, witness: {}",
            x.as_canonical_u64(),
            witness.as_canonical_u64()
        );
        if self.a == F::zero() {
            self.a = witness;
            self.q_l = x;
        } else if self.b == F::zero() {
            self.b = witness;
            self.q_r = x;
        } else if self.c == F::zero() {
            self.c = witness;
            self.q_o = x;
        } else {
            unreachable!(
                "Maximum linear term amount reached. Make sure you are using --expression-width 3"
            );
        }
    }
}

#[derive(Clone)]
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

struct PlonkAir;

pub struct PlonkBuilder<F: PrimeField64> {
    rows: Vec<PlonkRow<F>>,
}

/// The modulus of the field.
pub const P: u32 = 15 * (1 << 27) + 1;

/// The modulus of the field as a u64.
const P_U64: u64 = P as u64;

pub fn ftbb(field: FieldElement) -> BabyBear {
    BabyBear::from_canonical_u32(field.to_u128() as u32)
}

impl<F: PrimeField64> PlonkBuilder<F> {
    fn new() -> PlonkBuilder<F> {
        PlonkBuilder { rows: vec![] }
    }

    pub fn from_acir_program(
        witness_map: &WitnessMap,
        opcodes: &[Opcode],
    ) -> PlonkBuilder<BabyBear> {
        let mut rows: Vec<PlonkRow<BabyBear>> = vec![];
        for opcode in opcodes {
            match opcode {
                Opcode::AssertZero(exp) => {
                    let mut row: PlonkRow<BabyBear> = PlonkRow::default();

                    // Handle mul part of the gate
                    if !exp.mul_terms.is_empty() {
                        let mul_term = &exp.mul_terms[0];
                        row.q_m = ftbb(mul_term.0);

                        let w_l = &mul_term.1;
                        row.a = ftbb(*witness_map.get(w_l).unwrap());

                        let w_r = &mul_term.2;
                        row.b = ftbb(*witness_map.get(w_r).unwrap());
                    }

                    // Handle linear combination part of the gate
                    for term in &exp.linear_combinations {
                        row.set_linear_term(ftbb(term.0), ftbb(*witness_map.get(&term.1).unwrap()));
                    }

                    // Set constant term
                    row.q_c = ftbb(exp.q_c);

                    rows.push(row);
                }

                Opcode::BrilligCall { .. } => {
                    println!(
                        "Ignored brillig (unconstrained) call in opcode interpreter (TODO: Remove)"
                    );
                }
                op @ (Opcode::BlackBoxFuncCall(_)
                | Opcode::Directive(_)
                | Opcode::MemoryOp { .. }
                | Opcode::MemoryInit { .. }
                | Opcode::Call { .. }) => {
                    panic!("Unhandled op {:?}", op);
                }
            };
        }

        PlonkBuilder { rows }
    }

    pub fn compile(self) -> Vec<PlonkRow<F>> {
        self.rows
    }
}

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
                + row.q_c,
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

pub fn generate_trace_from_plonk_rows<F: PrimeField64>(
    plonk_rows: &mut Vec<PlonkRow<F>>,
) -> RowMajorMatrix<F> {
    // NB: N needs to be a power of two
    let mut n = plonk_rows.len().next_power_of_two();
    if n == 1 {
        n = 2;
    }
    let mut trace = RowMajorMatrix::new(vec![F::zero(); n * PLONK_GATE_WIDTH], PLONK_GATE_WIDTH);
    let (prefix, rows, suffix) = unsafe { trace.values.align_to_mut::<PlonkRow<F>>() };

    assert!(prefix.is_empty(), "Alignment check! Ethereum aligned??");
    assert!(suffix.is_empty(), "Alignment check! Ethereum aligned??");

    /* Insert the rows */
    for i in 0..plonk_rows.len() {
        rows[i] = plonk_rows[i].clone();
    }

    // Generate sub groups for copy constraints
    let n_f = F::from_canonical_usize(n);
    let two = F::from_canonical_u16(2);
    let three = F::from_canonical_usize(3);

    for (i, row) in rows.iter_mut().enumerate() {
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
    let mut numerator = Vec::with_capacity(n);
    numerator.fill(F::zero());
    unsafe { numerator.set_len(n) }
    let mut denominator = Vec::with_capacity(n);
    denominator.fill(F::zero());
    unsafe { denominator.set_len(n) }

    for i in 0..n {
        numerator[i] = calculate_numerator(&rows[i]);
        // TODO: batch inverse
        denominator[i] = calculate_denominator(&rows[i]);
    }

    for i in 0..n - 1 {
        // Calculate running numerator and denominator products
        let num = numerator[i];
        numerator[i + 1] *= num;

        let den = denominator[i];
        denominator[i + 1] *= den;
    }

    // Calculate grand product accumulator
    rows[0].copy.acc = F::one();
    for i in 1..n {
        rows[i].copy.acc = numerator[i] * denominator[i];
    }

    check_constraints(&PlonkAir, &trace, &vec![]);

    trace
}

pub fn generate_trace<F: PrimeField64>(n: usize) -> RowMajorMatrix<F> {
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

    for (i, row) in rows.iter_mut().enumerate() {
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

pub fn prove_and_verify(trace: &DenseMatrix<BabyBear>) {
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
    verify(&config, &PlonkAir, &mut challenger, &proof, &vec![]).expect("verification failed");
}

#[cfg(test)]
mod test {
    use super::{generate_trace, generate_trace_from_plonk_rows, prove_and_verify, PlonkBuilder};
    use acir::circuit::Opcode;
    use acir::native_types::WitnessMap;
    use p3_baby_bear::BabyBear;

    #[test]
    fn test_simple_add_gates() {
        let trace = generate_trace(4);
        prove_and_verify(&trace);
    }

    #[test]
    fn test_acir() {
        let opcodes: Vec<Opcode> = vec![];
        let witness_map = WitnessMap::new();
        let dc_builder: PlonkBuilder<BabyBear> =
            PlonkBuilder::<BabyBear>::from_acir_program(&witness_map, &opcodes);
        let mut rows = dc_builder.compile();
        let trace = generate_trace_from_plonk_rows(&mut rows);
        prove_and_verify(&trace);
    }
}
