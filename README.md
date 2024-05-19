# Air Fried Gyatt

## How To Run

1. Install `noirup`: https://github.com/noir-lang/noirup
2. Clone pinned noir fork: `https://github.com/Philogy/noir-snapshot`
3. CD into the cloned Noir repo and run `noirup -p .`
4. Compile & generate witnesses for your circuit using `nargo`
    1. `nargo compile --expression-width 3`
    2. `nargo execute witness.gz`
5. Run backend with `cargo run --release --bin gyatt-fryer-cli -- target/wittness.gz target/mini_circuit.json`

# Notes

You need to pass `--expression-width 3` to `nargo compile`.
