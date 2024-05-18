use acir::circuit::Opcode;
use acir::native_types::WitnessStack;
use clap::Parser;
use nargo::artifacts::program::ProgramArtifact;
use serde_json;
use std::fs;

#[derive(Parser, Debug)]
#[clap(about = "Vape Disposal 6000 Noir Backend")]
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
    let witnesses = witness_stack.peek().expect("Empty witness stack");

    dbg!(&witnesses);

    let contents = fs::read_to_string(args.acir_path)?;
    let program: ProgramArtifact = serde_json::from_str(contents.as_str())?;

    let func = &program.bytecode.functions[0];

    for op in &func.opcodes {
        match op {
            Opcode::AssertZero(expr) => {
                dbg!(&expr);
            }
            Opcode::BrilligCall { .. } => {}
            _ => todo!(),
        }
    }

    Ok(())
}
