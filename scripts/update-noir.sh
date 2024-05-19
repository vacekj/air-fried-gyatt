#!/bin/bash
cd mini_circuit || exit
nargo compile --expression-width 3
nargo execute witty.gz
cd ..
