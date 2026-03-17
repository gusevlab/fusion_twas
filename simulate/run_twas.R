# run_twas.R
# -----------------------------------------------------------------------
# End-to-end TWAS validation pipeline using the ACTUAL FUSION scripts:
#
#   Stage 1 — FUSION.compute_weights.R
#     For each simulated gene, extracts the window genotypes via PLINK,
#     writes the simulated expression as a phenotype file, then calls
#     FUSION.compute_weights.R to fit top1 / enet / lasso models with
#     5-fold CV and estimate heritability via GCTA REML.
#
#   Stage 2 — FUSION.assoc_test.R
#     Calls FUSION.assoc_test.R with the fitted weights and GWAS sumstats
#     to compute TWAS Z-scores, exactly as in production.
#
#   Stage 3 — Comparison to ground truth
#     Reads FUSION.assoc_test.R's output table, computes expected TWAS Z
#     from the known true weights, and reports TP/FP/FN at Bonferroni α.
#     Saves a scatter plot: Expected vs Observed TWAS Z.
#
# Prerequisites:
#   plink 1.9  : --PATH_plink (default /tmp/plink_mac/plink)
#   gcta 1.95  : --PATH_gcta  (default /tmp/gcta_wrapper.sh)
#   R packages : plink2R, optparse, glmnet
#
# Usage:
#   cd simulate/
#   Rscript run_twas.R [options]
# -----------------------------------------------------------------------

suppressMessages(library(optparse))
local({
  f = grep("--file=", commandArgs(FALSE), value = TRUE)
  d = if (length(f)) { me = normalizePath(sub("--file=","",f[1])); dd = dirname(me); if (basename(dd) %in% c("utils","simulate")) dirname(dd) else dd } else getwd()
  source(file.path(d, "utils", "plink_utils.R"))
})

option_list = list(
  make_option("--bfile", default = "../LDREF/1000G.EUR.22", type = 'character',
              help = "PLINK binary prefix for the chr22 reference panel [default: %default]"),
  make_option("--out_dir", default = "output", type = 'character',
              help = "Directory produced by simulate_data.R [default: %default]"),
  make_option("--fusion_dir", default = "..", type = 'character',
              help = "Root of FUSION repository (contains FUSION.compute_weights.R) [default: %default]"),
  make_option("--PATH_plink", default = "plink", type = 'character',
              help = "Path to PLINK 1.9 executable [default: %default]"),
  make_option("--PATH_gcta", default = "gcta", type = 'character',
              help = "Path to GCTA executable (or wrapper) [default: %default]"),
  make_option("--models", default = "top1,enet,lasso", type = 'character',
              help = "Models for FUSION.compute_weights.R [default: %default]"),
  make_option("--hsq_p", default = 1.0, type = 'double',
              help = "Heritability p-value threshold (1.0 = keep all genes) [default: %default]"),
  make_option("--hsq_set", default = 0.10, type = 'double',
              help = "Fix h² to this value in FUSION.compute_weights.R (skips slow GCTA REML). Set to NA to run GCTA. [default: %default]"),
  make_option("--crossval", default = 5L, type = 'integer',
              help = "Cross-validation folds [default: %default]"),
  make_option("--min_r2pred", default = 0.7, type = 'double',
              help = "Min LD-imputation r² in FUSION.assoc_test.R [default: %default]"),
  make_option("--coloc_susie_P", default = NA, type = 'double', action = "store",
              help = "P-value threshold below which to run coloc.susie in FUSION.assoc_test.R [default: off]"),
  make_option("--coloc_P", default = NA, type = 'double', action = "store",
              help = "P-value threshold below which to run coloc in FUSION.assoc_test.R [default: off]"),
  make_option("--seed", default = 42L, type = 'integer',
              help = "Random seed passed to FUSION.compute_weights.R [default: %default]")
)

opt = parse_args(OptionParser(option_list = option_list))

chr_prefix = function(bfile, chr) {
  if (grepl("\\.[0-9XYM]+$", bfile)) {
    out = sub("\\.[0-9XYM]+$", paste0(".", chr), bfile)
  } else {
    out = paste0(bfile, ".", chr)
  }
  out
}

rbind_fill = function(tbls) {
  cols = unique(unlist(lapply(tbls, names)))
  tbls = lapply(tbls, function(df) {
    miss = setdiff(cols, names(df))
    for (nm in miss) df[[nm]] = NA
    df[, cols, drop = FALSE]
  })
  do.call(rbind, tbls)
}

# -----------------------------------------------------------------------
# Paths to the two FUSION scripts
# -----------------------------------------------------------------------
script_weights = file.path(opt$fusion_dir, "FUSION.compute_weights.R")
script_assoc   = file.path(opt$fusion_dir, "FUSION.assoc_test.R")

for (sc in c(script_weights, script_assoc)) {
  if (!file.exists(sc)) stop("Cannot find FUSION script: ", sc)
}

# Validate external tools
for (tool in c(opt$PATH_plink, opt$PATH_gcta)) {
  rc = system(paste(tool, "--help"), ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (rc != 0) {
    if (!file.exists(tool)) stop("Tool not found: ", tool)
  }
}
cat(sprintf("Using plink : %s\n", opt$PATH_plink))
cat(sprintf("Using gcta  : %s\n", opt$PATH_gcta))

# -----------------------------------------------------------------------
# Load simulation ground truth + reference panel sample IDs
# -----------------------------------------------------------------------
cat("\nLoading simulation data ...\n")
truth       = readRDS(file.path(opt$out_dir, "simulation_truth.rds"))
expr_matrix = readRDS(file.path(opt$out_dir, "expr_matrix.rds"))
gene_list   = truth$gene_list
n_genes     = nrow(gene_list)
gene_chrs   = unique(as.character(gene_list$CHR))

# We need the FAM file sample IDs to build phenotype files
ref_bfile = chr_prefix(opt$bfile, gene_chrs[1])
genos_obj = read_plink(ref_bfile, impute = "avg")
fam       = genos_obj$fam          # N x 6 (FID, IID, ...)
N_ref     = nrow(fam)

# -----------------------------------------------------------------------
# Stage 1 — Run FUSION.compute_weights.R for every gene
# -----------------------------------------------------------------------
cat(sprintf("\n=== Stage 1: FUSION.compute_weights.R (%d genes) ===\n", n_genes))

plink_dir  = file.path(opt$out_dir, "plink_genes")
wgt_dir    = file.path(opt$out_dir, "weights_fusion")
tmp_base   = file.path(opt$out_dir, "tmp")
dir.create(plink_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(wgt_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_base,  recursive = TRUE, showWarnings = FALSE)

for (g in seq_len(n_genes)) {
  gene_id = gene_list$ID[g]
  gene_chr = as.character(gene_list$CHR[g])
  gene_bfile = chr_prefix(opt$bfile, gene_chr)
  window_snps = truth$true_weights[[g]]$snp_ids

  gene_plink_prefix = file.path(plink_dir, gene_id)
  gene_pheno_file   = paste0(gene_plink_prefix, ".pheno")
  gene_wgt_prefix   = file.path(wgt_dir, gene_id)
  gene_tmp_prefix   = file.path(tmp_base, gene_id)

  # --- 1a. Write SNP-list file for this window ---
  snp_list_file = paste0(gene_plink_prefix, ".snplist")
  write.table(window_snps, file = snp_list_file,
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  # --- 1b. Extract window genotypes with PLINK ---
  plink_cmd = paste(
    opt$PATH_plink,
    "--allow-no-sex",
    "--bfile", gene_bfile,
    "--extract", snp_list_file,
    "--make-bed",
    "--out", gene_plink_prefix
  )
  rc = system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  if (rc != 0 || !file.exists(paste0(gene_plink_prefix, ".bed"))) {
    warning(sprintf("PLINK failed for %s, skipping.", gene_id))
    next
  }

  # --- 1c. Write phenotype file (FID IID expr) ---
  pheno_tbl = data.frame(FID = fam[, 1], IID = fam[, 2],
                         PHENO = expr_matrix[, g])
  write.table(pheno_tbl, file = gene_pheno_file,
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  # --- 1d. Call FUSION.compute_weights.R ---
  hsq_set_flag = if (!is.na(opt$hsq_set)) paste("--hsq_set", opt$hsq_set) else ""

  fusion_wgt_cmd = paste(
    "Rscript", script_weights,
    "--bfile",      gene_plink_prefix,
    "--pheno",      gene_pheno_file,
    "--out",        gene_wgt_prefix,
    "--tmp",        gene_tmp_prefix,
    "--PATH_plink", opt$PATH_plink,
    "--PATH_gcta",  opt$PATH_gcta,
    "--models",     opt$models,
    "--hsq_p",      opt$hsq_p,
    hsq_set_flag,
    "--crossval",   opt$crossval,
    "--verbose 0",
    "--save_hsq"
  )
  rc = system(fusion_wgt_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  if (file.exists(paste0(gene_wgt_prefix, ".wgt.RDat"))) {
    if (g %% 10 == 0 || g == n_genes)
      cat(sprintf("  [%d/%d] %s — weights computed\n", g, n_genes, gene_id))
  } else {
    cat(sprintf("  [%d/%d] %s — skipped by FUSION (low h² or convergence)\n",
                g, n_genes, gene_id))
  }
}

# -----------------------------------------------------------------------
# Build a gene list containing only genes whose weight files exist
# -----------------------------------------------------------------------
gene_list_fusion = gene_list
gene_list_fusion$WGT = paste0(gene_list$ID, ".wgt.RDat")   # relative path
gene_list_fusion$exists = file.exists(
  file.path(wgt_dir, gene_list_fusion$WGT))

cat(sprintf("\nWeight files produced: %d / %d\n",
            sum(gene_list_fusion$exists), n_genes))

# Write FUSION-format gene list (WGT paths are relative to weights_dir)
gene_list_file = file.path(opt$out_dir, "fusion_gene_list.txt")
write.table(gene_list_fusion[gene_list_fusion$exists,
                              c("WGT","ID","CHR","P0","P1","N")],
            file = gene_list_file,
            quote = FALSE, row.names = FALSE, sep = "\t")

# -----------------------------------------------------------------------
# Stage 2 — Run FUSION.assoc_test.R
# -----------------------------------------------------------------------
cat("\n=== Stage 2: FUSION.assoc_test.R ===\n")

assoc_out = file.path(opt$out_dir, "twas_assoc")
assoc_tables = list()
susie_summary_tables = list()
susie_result_tables = list()

for (assoc_chr in unique(as.character(gene_list_fusion$CHR[gene_list_fusion$exists]))) {
  gene_list_chr = gene_list_fusion[gene_list_fusion$exists &
                                     as.character(gene_list_fusion$CHR) == assoc_chr,
                                   c("WGT","ID","CHR","P0","P1","N")]
  gene_list_chr_file = file.path(opt$out_dir, paste0("fusion_gene_list.chr", assoc_chr, ".txt"))
  write.table(gene_list_chr, file = gene_list_chr_file,
              quote = FALSE, row.names = FALSE, sep = "\t")

  assoc_out_chr = paste0(assoc_out, ".chr", assoc_chr)
  ref_ld_chr = sub(sprintf("\\.%s$", assoc_chr), ".", chr_prefix(opt$bfile, assoc_chr))
  if (identical(ref_ld_chr, chr_prefix(opt$bfile, assoc_chr))) {
    stop("Could not derive --ref_ld_chr prefix from --bfile: ", opt$bfile)
  }

  fusion_assoc_cmd = paste(
    "Rscript", script_assoc,
    "--sumstats",   file.path(opt$out_dir, "sumstats.txt"),
    "--weights",    gene_list_chr_file,
    "--weights_dir", wgt_dir,
    "--ref_ld_chr", ref_ld_chr,
    "--chr",        assoc_chr,
    "--min_r2pred", opt$min_r2pred,
    if (!is.na(opt$coloc_susie_P)) paste("--coloc_susie_P", opt$coloc_susie_P) else "",
    if (!is.na(opt$coloc_P)) paste("--coloc_P", opt$coloc_P) else "",
    "--GWASN",      truth$opt$N_gwas,
    "--out",        assoc_out_chr
  )
  cat("Running:\n ", fusion_assoc_cmd, "\n")
  rc = system(fusion_assoc_cmd)
  if (rc != 0) warning("FUSION.assoc_test.R returned non-zero exit code for chr", assoc_chr)

  if (!file.exists(assoc_out_chr)) {
    stop("FUSION.assoc_test.R produced no output at: ", assoc_out_chr)
  }

  assoc_tables[[assoc_chr]] = read.table(assoc_out_chr, head = TRUE, as.is = TRUE, sep = "\t")
  susie_summary_file = paste0(assoc_out_chr, ".coloc.susie.summary")
  susie_results_file = paste0(assoc_out_chr, ".coloc.susie.results")
  if (file.exists(susie_summary_file)) {
    susie_summary_tables[[assoc_chr]] = read.table(susie_summary_file, head = TRUE, as.is = TRUE, sep = "\t")
  }
  if (file.exists(susie_results_file)) {
    susie_result_tables[[assoc_chr]] = read.table(susie_results_file, head = TRUE, as.is = TRUE, sep = "\t")
  }
}

twas_out = rbind_fill(assoc_tables)
write.table(twas_out, file = assoc_out, quote = FALSE, row.names = FALSE, sep = "\t")
if (length(susie_summary_tables) > 0) {
  write.table(rbind_fill(susie_summary_tables),
              file = paste0(assoc_out, ".coloc.susie.summary"),
              quote = FALSE, row.names = FALSE, sep = "\t")
}
if (length(susie_result_tables) > 0) {
  write.table(rbind_fill(susie_result_tables),
              file = paste0(assoc_out, ".coloc.susie.results"),
              quote = FALSE, row.names = FALSE, sep = "\t")
}

# -----------------------------------------------------------------------
# Stage 3 — Parse output and compare to ground truth
# -----------------------------------------------------------------------
cat("\n=== Stage 3: Results ===\n")

# Annotate causal status
twas_out$CAUSAL = twas_out$ID %in% truth$causal_gene_ids

# Compute theoretical expected TWAS Z for each gene using TRUE weights
# E[TWAS.Z] = beta_gene * sqrt(N_gwas) * sqrt(w_true' * LD * w_true)
genos_by_chr = lapply(gene_chrs, function(chr) {
  g = read_plink(chr_prefix(opt$bfile, chr), impute = "avg")
  G = scale(g$bed)
  na_cols = which(apply(is.na(G), 2, any))
  if (length(na_cols) > 0) G[, na_cols] = 0
  list(G = G, bim = g$bim)
})
names(genos_by_chr) = gene_chrs

twas_out$TWAS.Z.EXPECTED = NA_real_

for (g in seq_len(nrow(twas_out))) {
  g_idx = match(twas_out$ID[g], gene_list$ID)
  if (is.na(g_idx)) next
  tw       = truth$true_weights[[g_idx]]
  beta_g   = truth$gene_effects[g_idx]   # signed: ±beta_gene or 0
  if (beta_g == 0) {
    twas_out$TWAS.Z.EXPECTED[g] = 0
    next
  }
  # Project true weight SNPs into the reference LD
  w_global  = tw$w_full                        # length = window_size
  chr_key   = as.character(tw$chr)
  G_all     = genos_by_chr[[chr_key]]$G
  bim_chr   = genos_by_chr[[chr_key]]$bim
  snp_ids   = tw$snp_ids                       # SNP names in this window
  m_in_ref  = match(snp_ids, bim_chr[, 2])     # column indices in chr-specific G_all
  m_ok      = !is.na(m_in_ref)
  if (sum(m_ok) < 2) next
  G_win     = G_all[, m_in_ref[m_ok], drop = FALSE]
  LD_win    = t(G_win) %*% G_win / (N_ref - 1)
  w_ok      = w_global[m_ok]
  quad      = as.numeric(w_ok %*% LD_win %*% w_ok)
  twas_out$TWAS.Z.EXPECTED[g] =
    round(beta_g * sqrt(truth$opt$N_gwas) * sqrt(max(quad, 0)), 2)
}

# Sort by p-value
twas_out = twas_out[order(as.numeric(twas_out$TWAS.P), na.last = TRUE), ]

cat("\nTop 20 TWAS hits (sorted by p-value):\n")
print(
  head(twas_out[, c("ID","CAUSAL","MODEL","MODELCV.R2",
                    "TWAS.Z","TWAS.Z.EXPECTED","TWAS.P")], 20),
  digits = 3, row.names = FALSE
)

# Detection at Bonferroni threshold
n_tested    = nrow(twas_out)
bonf_thresh = 0.05 / n_tested
sig         = twas_out$ID[!is.na(as.numeric(twas_out$TWAS.P)) &
                           as.numeric(twas_out$TWAS.P) < bonf_thresh]
tp = sum(sig %in% truth$causal_gene_ids)
fp = sum(!sig %in% truth$causal_gene_ids)
fn = sum(!truth$causal_gene_ids %in% sig)

cat(sprintf("\n--- Detection (Bonferroni p < %.2e across %d genes) ---\n",
            bonf_thresh, n_tested))
cat(sprintf("  True positives  : %d / %d causal\n", tp, length(truth$causal_gene_ids)))
cat(sprintf("  False positives : %d\n", fp))
cat(sprintf("  False negatives : %d\n", fn))

cat("\n--- Observed vs Expected TWAS Z (causal genes only) ---\n")
causal_rows = twas_out[twas_out$CAUSAL, ]
print(causal_rows[, c("ID","TWAS.Z","TWAS.Z.EXPECTED","MODEL","MODELCV.R2","HSQ")],
      digits = 3, row.names = FALSE)

# Save annotated results
write.table(twas_out,
            file  = file.path(opt$out_dir, "twas_results_fusion.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")
cat(sprintf("\nAnnotated results saved to %s/twas_results_fusion.txt\n", opt$out_dir))

if ("COLOC.SUSIE.PP4" %in% names(twas_out)) {
  parse_max_pp4 = function(x) {
    if (is.na(x) || !nzchar(x)) return(NA_real_)
    vals = suppressWarnings(as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]))
    vals = vals[!is.na(vals)]
    if (length(vals) == 0L) return(NA_real_)
    max(vals)
  }
  twas_out$COLOC.SUSIE.PP4.MAX = vapply(twas_out$COLOC.SUSIE.PP4, parse_max_pp4, numeric(1))
}

# -----------------------------------------------------------------------
# Scatter plot: Expected vs Observed TWAS Z-scores
# -----------------------------------------------------------------------
cat("\nGenerating scatter plot ...\n")

plot_df = twas_out[!is.na(twas_out$TWAS.Z.EXPECTED) &
                   !is.na(suppressWarnings(as.numeric(twas_out$TWAS.Z))), ]
plot_df$TWAS.Z.num = as.numeric(plot_df$TWAS.Z)

# Pearson r for causal and null genes
r_causal = if (sum(plot_df$CAUSAL, na.rm = TRUE) >= 2) {
  cor(plot_df$TWAS.Z.num[plot_df$CAUSAL],
      plot_df$TWAS.Z.EXPECTED[plot_df$CAUSAL],
      use = "complete.obs")
} else {
  NA_real_
}
r_null   = if (sum(!plot_df$CAUSAL, na.rm = TRUE) >= 2) {
  cor(plot_df$TWAS.Z.num[!plot_df$CAUSAL],
      plot_df$TWAS.Z.EXPECTED[!plot_df$CAUSAL],
      use = "complete.obs")
} else {
  NA_real_
}
r_all    = cor(plot_df$TWAS.Z.num, plot_df$TWAS.Z.EXPECTED,
               use = "complete.obs")

# Bonferroni Z threshold (two-sided)
bonf_z = qnorm(1 - 0.05 / (2 * n_tested))

plot_file = file.path(opt$out_dir, "twas_z_scatter.pdf")
pdf(plot_file, width = 6, height = 6)

  col_causal = "#D62728"
  col_null   = "#AAAAAA"
  col_vec    = ifelse(plot_df$CAUSAL, col_causal, col_null)
  pch_vec    = ifelse(plot_df$CAUSAL, 19, 1)
  cex_vec    = ifelse(plot_df$CAUSAL, 0.85, 0.55)

  axis_lim = ceiling(max(abs(c(plot_df$TWAS.Z.EXPECTED,
                                plot_df$TWAS.Z.num)), na.rm = TRUE) * 1.05)

  plot(plot_df$TWAS.Z.EXPECTED, plot_df$TWAS.Z.num,
       col  = col_vec, pch = pch_vec, cex = cex_vec,
       xlim = c(-axis_lim, axis_lim),
       ylim = c(-axis_lim, axis_lim),
       xlab = "Expected TWAS Z  (true weights)",
       ylab = "Observed TWAS Z  (FUSION estimated weights)",
       main = sprintf("TWAS Validation — Expected vs Observed Z\n%d genes, %d causal  |  h\u00b2=%.2f  |  N=%s",
                      nrow(plot_df), sum(plot_df$CAUSAL),
                      truth$opt$hsq,
                      formatC(truth$opt$N_gwas, format = "d", big.mark = ",")))

  abline(0, 1,  lty = 2, col = "navy",        lwd = 1.5)
  abline(h = 0, lty = 1, col = "gray70",      lwd = 0.8)
  abline(v = 0, lty = 1, col = "gray70",      lwd = 0.8)
  abline(h =  bonf_z, lty = 3, col = "forestgreen", lwd = 1.2)
  abline(h = -bonf_z, lty = 3, col = "forestgreen", lwd = 1.2)

  legend_entries = c(sprintf("Causal  (r = %s)",
                             ifelse(is.na(r_causal), "NA", sprintf("%.3f", r_causal))),
                     "y = x  (ideal)",
                     "Bonferroni \u00b1Z")
  legend_cols = c(col_causal, "navy", "forestgreen")
  legend_pch = c(19, NA, NA)
  legend_lty = c(NA, 2, 3)
  legend_lwd = c(NA, 1.5, 1.2)
  if (!is.na(r_null)) {
    legend_entries = append(legend_entries,
                            sprintf("Null    (r = %.3f)", r_null),
                            after = 1)
    legend_cols = append(legend_cols, col_null, after = 0)
    legend_pch = append(legend_pch, 1, after = 0)
    legend_lty = append(legend_lty, NA, after = 0)
    legend_lwd = append(legend_lwd, NA, after = 0)
  }
  legend("topleft",
         legend = legend_entries,
         col  = legend_cols,
         pch  = legend_pch,
         lty  = legend_lty,
         lwd  = legend_lwd,
         bty  = "n", cex = 0.8)

  mtext(sprintf("Overall r = %.3f", r_all), side = 3, adj = 1, line = 0.3, cex = 0.78)

dev.off()
cat(sprintf("Scatter plot saved to %s\n", plot_file))

# -----------------------------------------------------------------------
# Plot: max coloc.susie PP4 by gene causal status
# -----------------------------------------------------------------------
if ("COLOC.SUSIE.PP4.MAX" %in% names(twas_out)) {
  cat("\nGenerating coloc.susie PP4 plot ...\n")

  pp4_df = twas_out[!is.na(twas_out$COLOC.SUSIE.PP4.MAX), ]
  pp4_file = file.path(opt$out_dir, "coloc_susie_max_pp4.pdf")

  pdf(pp4_file, width = 7, height = 5)

    causal_pp4 = pp4_df$COLOC.SUSIE.PP4.MAX[pp4_df$CAUSAL]
    null_pp4 = pp4_df$COLOC.SUSIE.PP4.MAX[!pp4_df$CAUSAL]
    grp = factor(ifelse(pp4_df$CAUSAL, "Causal", "Non-causal"),
                 levels = c("Causal", "Non-causal"))
    cols = c("Causal" = "#D62728", "Non-causal" = "#7F7F7F")

    boxplot(COLOC.SUSIE.PP4.MAX ~ grp, data = pp4_df,
            col = cols[levels(grp)],
            border = cols[levels(grp)],
            ylim = c(0, 1),
            ylab = "Max coloc.susie PP4 per gene",
            xlab = "",
            main = sprintf("coloc.susie max PP4 by causal status\n%d genes, %d causal, %d non-causal",
                           nrow(pp4_df), sum(pp4_df$CAUSAL), sum(!pp4_df$CAUSAL)))
    stripchart(COLOC.SUSIE.PP4.MAX ~ grp, data = pp4_df,
               vertical = TRUE, method = "jitter", jitter = 0.15,
               pch = 19, cex = 0.7, add = TRUE,
               col = adjustcolor(cols[as.character(grp)], alpha.f = 0.55))
    abline(h = c(0.1, 0.5, 0.9), col = "gray85", lty = 3)

    legend_entries = character(0)
    if (length(causal_pp4) > 0) {
      legend_entries = c(legend_entries,
                         sprintf("Causal median = %.3f", median(causal_pp4)))
    }
    if (length(null_pp4) > 0) {
      legend_entries = c(legend_entries,
                         sprintf("Non-causal median = %.3f", median(null_pp4)))
    }
    if (length(legend_entries) > 0) {
      legend("topright", legend = legend_entries, bty = "n", cex = 0.85)
    }

  dev.off()
  cat(sprintf("coloc.susie PP4 plot saved to %s\n", pp4_file))
}
