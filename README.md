# Bloody Mary 600 Disposable

## How To Run

1. Install `noirup`: https://github.com/noir-lang/noirup
2. Clone pinned noir fork: `https://github.com/Philogy/noir-snapshot`
3. CD into the repo and run `noirup -p .`
4. Compile & generate witnesses for our circuit under `mini_circuit` using `./scripts/update-noir.sh`
5. Run backend with `cargo run --bin smoke-stack -- mini_circuit/target/witty.gz mini_circuit/target/mini_circuit.json`

