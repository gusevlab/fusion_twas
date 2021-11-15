# Utilities for FUSION analyses

## `compare_models.R` : Compute model correlations
Reads in two lists of FUSION models and computes pairwise correlations between all pairs of nearby models. Execution is similar to main FUSION anlayses: `Rscript compare_models.R --pos1 TCGA-BRCA.GE.TUMOR.pos --pos2 TCGA-BRCA.GE.NORMAL.pos --chr 1 --ref_ld_chr LDREF/1000G.EUR. --window 100000`. Where `--pos1` and `--pos2` point to model position files that must contain the labeled columns (`CHR`, `P1`, `P2,` `WGT`); `ref_ld_chr` points to the choromosome level LD reference data (as with main FUSION analyses); `--window` specifies the window to place around model boundaries to consider "overlapping". The output is printed to screen and includes: `Model 1`, `Model 2`, `Distance`, `Correlation`. 

## `make_score.R` : Make individual-level predictors
Reads in a FUSION model and generates a score file for individual level prediction. `Rscript make_score.R [wgt.RDat file] > [SCORE_FILE]`. To predict into individual level data, the following command can then be used: `plink --bfile [GENOTYPES] --score [SCORE_FILE] 1 2 4`.

## `FUSION.profile_wgt.R` : Profile models
Reads a list of FUSION models and outputs various metrics. Run as `Rscript FUSION.profile_wgt.R [LIST] > [OUTPUT]` where `[LIST]` is a file listing `*.wgt.RDat` files for which parameters are to be extracted. The run will output a table of results to `[OUTPUT]` as well as a summary of means to the standard error.
