# Air Fried Gyatt

## How To Run

1. Install `noirup`: https://github.com/noir-lang/noirup
2. Clone pinned noir fork: `https://github.com/Philogy/noir-snapshot`
3. CD into the cloned Noir repo and run `noirup -p .`
4. Compile & generate witnesses for our circuit under `mini_circuit` using `./scripts/update-noir.sh`
5. Run backend with `cargo run --bin gyatt-fryer-cli -- mini_circuit/target/witty.gz mini_circuit/target/mini_circuit.json`

# Notes

You need to pass `--expression-width 3` to `nargo compile`.