use acir::circuit::Opcode;
use acir::native_types::WitnessStack;
use circuit_checker::{generate_trace_from_plonk_rows, prove_and_verify, BabyBear, PlonkBuilder};
use clap::Parser;
use nargo::artifacts::program::ProgramArtifact;
use serde_json;
use std::fs;

#[derive(Parser, Debug)]
#[clap(about = "Air Fried Booty Noir Backend")]
struct Args {
    #[clap(index = 1)]
    circuit_path: String,
    #[clap(index = 2)]
    acir_path: String,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("Breathe Air");

    let args = Args::parse();
    println!("args: {:?}", args);

    let contents = fs::read(args.circuit_path)?;
    let witness_stack: WitnessStack = contents.as_slice().try_into()?;
    let witnesses = &witness_stack.peek().expect("Empty witness stack").witness;

    let contents = fs::read_to_string(args.acir_path)?;
    let program: ProgramArtifact = serde_json::from_str(contents.as_str())?;

    let func = &program.bytecode.functions[0];
    dbg!(witnesses);
    let opcodes = &func.opcodes;
    dbg!(opcodes);
    // for op in opcodes.clone() {
    //     match op {
    //         Opcode::AssertZero(exp) => {
    //             dbg!(exp);
    //         }
    //         _ => {}
    //     }
    // }

    let builder = PlonkBuilder::<BabyBear>::from_acir_program(witnesses, opcodes);
    let mut gates = builder.compile();
    println!("q_m * (a * b) + q_l * a + q_r * b + q_o * c + q_c = 0");
    for gate in &gates {
        println!("{}", gate);
    }

    let trace = generate_trace_from_plonk_rows(&mut gates);

    prove_and_verify(&trace);
    println!("BRO WE VERIFIED");
    Ok(())
}
