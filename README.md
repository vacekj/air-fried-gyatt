# Air Fried Gyatt

## How To Run

1. Install `noirup`: https://github.com/noir-lang/noirup
2. Clone pinned noir fork: `https://github.com/Philogy/noir-snapshot`
3. CD into the cloned Noir repo and run `noirup -p .`
4. Compile & generate witnesses for your circuit using `nargo` (cd into `mini_circuit` for our example, CD back out after)
    1. `nargo compile --expression-width 3`
    2. `nargo execute witness.gz`
5. Run backend with `cargo run --bin gyatt-fryer-cli -- mini_circuit/target/wittness.gz mini_circuit/target/mini_circuit.json`
