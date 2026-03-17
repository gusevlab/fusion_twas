#!/bin/bash

echo "=== Step 1: Simulate data ==="
Rscript simulate_data.R --seed 1 --n_genes 30 --n_causal_genes 20 --N_gwas 200e3

echo "=== Step 2: Run FUSION pipeline ==="
Rscript run_twas.R --coloc_susie_P 0.05 --PATH_gcta /Users/alexandergusev/Downloads/gcta-1.95.1-macOS-arm64/bin/gcta64

echo "=== Pipeline complete ==="
