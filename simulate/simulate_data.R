# simulate_data.R
# -----------------------------------------------------------------------
# Simulates heritable gene expression and GWAS summary statistics for
# TWAS validation.
#
# Model:
#   expr_i  = G_causal_i %*% w  +  eps_expr       (h² = hsq)
#   y_i     = sum_g( beta_g * expr_g_i ) + eps_y   (for causal genes)
#   Z_j     = Z_j_ref + (sqrt(N_gwas) - sqrt(N_ref)) * (LD %*% b)_j
#             where b_j = sum_g( w_jg * beta_g )    (true per-allele effects)
#
# Output (all written to --out_dir):
#   gene_list.txt          - gene manifest for FUSION
#   weights/GENE*.wgt.RDat - TRUE weight files (sparse, known causal SNPs)
#   sumstats.txt           - GWAS summary statistics
#   expr_matrix.rds        - N_ref × G simulated expression matrix
#   simulation_truth.rds   - ground-truth object for evaluation
# -----------------------------------------------------------------------

suppressMessages(library(optparse))
local({
  f = grep("--file=", commandArgs(FALSE), value = TRUE)
  d = if (length(f)) { me = normalizePath(sub("--file=","",f[1])); dd = dirname(me); if (basename(dd) %in% c("utils","simulate")) dirname(dd) else dd } else getwd()
  source(file.path(d, "utils", "plink_utils.R"))
})

option_list = list(
  make_option("--bfile", default="../LDREF/1000G.EUR.22", type='character',
              help="PLINK binary prefix used to derive chr22/chr21 reference panels [default: %default]"),
  make_option("--out_dir", default="output", type='character',
              help="Output directory [default: %default]"),
  make_option("--n_genes", default=30L, type='integer',
              help="Number of gene windows to simulate [default: %default]"),
  make_option("--n_causal_genes", default=20L, type='integer',
              help="Number of genes with a causal effect on the trait [default: %default]"),
  make_option("--n_causal_eqtl", default=3L, type='integer',
              help="Number of causal cis-eQTLs per gene [default: %default]"),
  make_option("--hsq", default=0.10, type='double',
              help="Expression heritability (h²) [default: %default]"),
  make_option("--beta_gene", default=0.10, type='double',
              help="SD of N(0, beta_gene²) from which causal gene effects are sampled [default: %default]"),
  make_option("--N_gwas", default=500000L, type='integer',
              help="GWAS sample size used to scale summary-statistic noise [default: %default]"),
  make_option("--window_size", default=100L, type='integer',
              help="Number of SNPs per gene window [default: %default]"),
  make_option("--seed", default=42L, type='integer',
              help="Random seed for reproducibility [default: %default]")
)

opt = parse_args(OptionParser(option_list = option_list))
set.seed(opt$seed)

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$out_dir, "weights"), recursive = TRUE, showWarnings = FALSE)

chr_prefix = function(bfile, chr) {
  if (grepl("\\.[0-9XYM]+$", bfile)) {
    out = sub("\\.[0-9XYM]+$", paste0(".", chr), bfile)
  } else {
    out = paste0(bfile, ".", chr)
  }
  out
}

load_chr_panel = function(chr) {
  bfile = chr_prefix(opt$bfile, chr)
  panel = read_plink(bfile, impute = "avg")
  G = scale(panel$bed)
  na_cols = which(apply(is.na(G), 2, any))
  if (length(na_cols) > 0) G[, na_cols] = 0
  list(chr = chr, bfile = bfile, fam = panel$fam, bim = panel$bim, G = G)
}

causal_chr = 22L
null_chr = 21L

# -----------------------------------------------------------------------
# 1. Load reference genotypes
# -----------------------------------------------------------------------
cat("Loading reference panel genotypes ...\n")
panels = list(
  "22" = load_chr_panel(causal_chr),
  "21" = load_chr_panel(null_chr)
)

N_ref = nrow(panels[["22"]]$G)
if (nrow(panels[["21"]]$G) != N_ref) {
  stop("chr21 and chr22 reference panels have different sample counts")
}
if (!all(panels[["21"]]$fam[, 1] == panels[["22"]]$fam[, 1]) ||
    !all(panels[["21"]]$fam[, 2] == panels[["22"]]$fam[, 2])) {
  stop("chr21 and chr22 reference panels do not share the same sample ordering")
}

cat(sprintf("  %d samples, %d SNPs on chr22\n", N_ref, ncol(panels[["22"]]$G)))
cat(sprintf("  %d samples, %d SNPs on chr21\n", N_ref, ncol(panels[["21"]]$G)))

# -----------------------------------------------------------------------
# 2. Define non-overlapping gene windows
# -----------------------------------------------------------------------
opt$n_causal_genes = min(opt$n_causal_genes, opt$n_genes)
n_null_genes = opt$n_genes - opt$n_causal_genes

make_windows = function(P, n) {
  if (n == 0L) {
    return(list(starts = integer(0), ends = integer(0)))
  }
  n_windows = floor(P / opt$window_size)
  if (n > n_windows) {
    stop(sprintf("Requested %d windows but only %d fit on this chromosome", n, n_windows))
  }
  starts = floor(seq(1, P - opt$window_size, length.out = n))
  list(starts = starts, ends = starts + opt$window_size - 1L)
}

windows_22 = make_windows(ncol(panels[["22"]]$G), opt$n_causal_genes)
windows_21 = make_windows(ncol(panels[["21"]]$G), n_null_genes)

causal_gene_idx   = seq_len(opt$n_causal_genes)
causal_gene_betas = rnorm(opt$n_causal_genes, mean = 0, sd = opt$beta_gene)
cat(sprintf("  Causal genes: %d  |  beta ~ N(0, %.2f²)  |  range [%.3f, %.3f]\n",
            opt$n_causal_genes, opt$beta_gene,
            min(causal_gene_betas), max(causal_gene_betas)))

# Gene manifest (FUSION weights-list format: WGT, ID, CHR, P0, P1, N)
gene_chr = c(rep(causal_chr, opt$n_causal_genes), rep(null_chr, n_null_genes))
window_starts = c(windows_22$starts, windows_21$starts)
window_ends = c(windows_22$ends, windows_21$ends)
p0 = c(panels[["22"]]$bim[windows_22$starts, 4], panels[["21"]]$bim[windows_21$starts, 4])
p1 = c(panels[["22"]]$bim[windows_22$ends, 4], panels[["21"]]$bim[windows_21$ends, 4])
gene_list = data.frame(
  WGT = paste0("weights/GENE", seq_len(opt$n_genes), ".wgt.RDat"),
  ID  = paste0("GENE", seq_len(opt$n_genes)),
  CHR = gene_chr,
  P0  = p0,
  P1  = p1,
  N   = N_ref,
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------
# 3. Simulate gene expression
# -----------------------------------------------------------------------
cat(sprintf("Simulating expression for %d genes (h²=%.2f, %d causal eQTLs each) ...\n",
            opt$n_genes, opt$hsq, opt$n_causal_eqtl))

# Full per-chromosome true per-SNP trait effects (accumulated across causal genes)
b_true = list(
  "22" = numeric(ncol(panels[["22"]]$G)),
  "21" = numeric(ncol(panels[["21"]]$G))
)

expr_matrix  = matrix(NA_real_, nrow = N_ref, ncol = opt$n_genes,
                      dimnames = list(NULL, gene_list$ID))
true_weights = vector("list", opt$n_genes)

for (g in seq_len(opt$n_genes)) {
  chr_key = as.character(gene_list$CHR[g])
  ws    = window_starts[g]
  we    = window_ends[g]
  G_chr = panels[[chr_key]]$G
  snp_info = panels[[chr_key]]$bim
  G_win = G_chr[, ws:we, drop = FALSE]    # N x window_size
  P_win = ncol(G_win)

  # ---- Pick causal eQTLs ----
  n_caus = min(opt$n_causal_eqtl, P_win)
  causal_local = sort(sample(P_win, n_caus))
  G_caus = G_win[, causal_local, drop = FALSE]

  # ---- Sample weights and scale to exact h² ----
  w_raw = rnorm(n_caus)
  g_pred = as.vector(G_caus %*% w_raw)
  w_scaled = w_raw * sqrt(opt$hsq / var(g_pred))

  # ---- Simulate expression (standardise to mean=0, sd=1) ----
  signal  = as.vector(G_caus %*% w_scaled)          # Var ≈ hsq
  epsilon = rnorm(N_ref, 0, sqrt(1 - opt$hsq))
  expr    = scale(signal + epsilon)
  expr_matrix[, g] = as.vector(expr)

  # ---- Store sparse full-window weight vector ----
  w_full = numeric(P_win)
  w_full[causal_local] = w_scaled

  snps = snp_info[ws:we, ]
  true_weights[[g]] = list(
    gene_id       = gene_list$ID[g],
    chr           = gene_list$CHR[g],
    window_snps   = ws:we,
    snp_ids       = snps[, 2],
    causal_local  = causal_local,         # local indices within window
    causal_global = (ws:we)[causal_local],# global indices in chromosome-specific G
    w_full        = w_full                # length = window_size, sparse
  )

  # ---- Accumulate per-SNP trait effects for causal genes ----
  which_causal = which(causal_gene_idx == g)
  if (length(which_causal) > 0) {
    b_true[[chr_key]][ws:we] = b_true[[chr_key]][ws:we] + w_full * causal_gene_betas[which_causal]
  }

  # ---- Save TRUE weight file in FUSION .wgt.RDat format ----
  wgt.matrix = matrix(w_full, nrow = P_win, ncol = 1,
                      dimnames = list(snps[, 2], "top1"))
  cv.performance = matrix(c(opt$hsq, 1e-4), nrow = 2, ncol = 1,
                          dimnames = list(c("rsq", "pval"), "top1"))
  hsq   = c(opt$hsq, 0.01)
  hsq.pv = 1e-4
  N.tot  = N_ref
  save(wgt.matrix, snps, cv.performance, hsq, hsq.pv, N.tot,
       file = file.path(opt$out_dir, gene_list$WGT[g]))
}

# -----------------------------------------------------------------------
# 4. Simulate GWAS summary statistics
#
# Method (preserves LD-correlated noise):
#   (a) Build individual-level trait from reference panel phenotypes
#   (b) Compute marginal Z-scores in the reference sample (N_ref = 489)
#       Z_ref_j = X_j' y / sqrt(N_ref)  =>  E[Z_ref_j] = sqrt(N_ref) * b_j
#   (c) Amplify to N_gwas while keeping LD-correlated residuals:
#       Z_gwas = Z_ref + (sqrt(N_gwas) - sqrt(N_ref)) * (LD %*% b_true)_j
#             where (LD %*% b_true) is estimated as X' X b_true / N_ref
# -----------------------------------------------------------------------
cat(sprintf("Simulating GWAS summary statistics (N_gwas = %d) ...\n", opt$N_gwas))

# (a) Individual-level trait
genetic_component = as.vector(panels[["22"]]$G %*% b_true[["22"]])
var_genetic = var(genetic_component)
var_env     = max(1 - var_genetic, 0.01)
y  = genetic_component + rnorm(N_ref, 0, sqrt(var_env))
y  = as.vector(scale(y))

# (b/c) Reference-panel marginal Z-scores and signal amplification by chromosome
sumstats_list = list()
signal_snps = 0L
top_gwas_z = 0
for (chr_key in names(panels)) {
  G_chr = panels[[chr_key]]$G
  snp_info = panels[[chr_key]]$bim
  Z_ref = as.vector(t(G_chr) %*% y) / sqrt(N_ref)
  b_ld = as.vector(t(G_chr) %*% (G_chr %*% b_true[[chr_key]])) / N_ref
  Z_gwas = Z_ref + (sqrt(opt$N_gwas) - sqrt(N_ref)) * b_ld
  signal_snps = signal_snps + sum(b_true[[chr_key]] != 0)
  top_gwas_z = max(top_gwas_z, max(abs(Z_gwas)))
  sumstats_list[[chr_key]] = data.frame(
    SNP = snp_info[, 2],
    CHR = snp_info[, 1],
    A1  = snp_info[, 5],
    A2  = snp_info[, 6],
    Z   = round(Z_gwas, 4),
    stringsAsFactors = FALSE
  )
}
sumstats = do.call(rbind, sumstats_list)

cat(sprintf("  Signal SNPs (|b_true|>0): %d\n", signal_snps))
cat(sprintf("  Top GWAS |Z|: %.2f\n", top_gwas_z))
write.table(sumstats, file = file.path(opt$out_dir, "sumstats.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

# -----------------------------------------------------------------------
# 5. Save remaining outputs
# -----------------------------------------------------------------------
write.table(gene_list, file = file.path(opt$out_dir, "gene_list.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")

gene_effects_signed              = setNames(numeric(opt$n_genes), gene_list$ID)
gene_effects_signed[causal_gene_idx] = causal_gene_betas

truth = list(
  gene_list        = gene_list,
  causal_gene_idx  = causal_gene_idx,
  causal_gene_ids  = gene_list$ID[causal_gene_idx],
  gene_effects     = gene_effects_signed,   # signed: ±beta_gene for causal, 0 for null
  true_weights     = true_weights,
  window_starts    = window_starts,
  window_ends      = window_ends,
  panels           = lapply(panels, function(x) list(chr = x$chr, bfile = x$bfile, bim = x$bim)),
  b_true           = b_true,
  opt              = opt
)
saveRDS(truth,       file = file.path(opt$out_dir, "simulation_truth.rds"))
saveRDS(expr_matrix, file = file.path(opt$out_dir, "expr_matrix.rds"))

cat("\n=== Simulation complete ===\n")
cat(sprintf("  Genes simulated  : %d\n",   opt$n_genes))
cat(sprintf("  Causal genes     : %s\n",   paste(gene_list$ID[causal_gene_idx], collapse=", ")))
cat(sprintf("  h² (expression)  : %.2f\n", opt$hsq))
cat(sprintf("  beta_gene        : %.2f\n", opt$beta_gene))
cat(sprintf("  N_gwas           : %d\n",   opt$N_gwas))
cat(sprintf("  Output directory : %s\n",   opt$out_dir))
