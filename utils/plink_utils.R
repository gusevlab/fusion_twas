# utils/plink_utils.R
# -----------------------------------------------------------------------
# Drop-in replacement for plink2R::read_plink().
# Zero external dependencies — uses only base R.
#
# Returns list(bed, bim, fam) with the same structure as plink2R:
#   $bed  N×P numeric matrix; dosage = copies of A1 (first allele in .bim)
#   $bim  P×6 data.frame  (CHR, SNP, CM, BP, A1, A2)
#   $fam  N×6 data.frame  (FID, IID, PID, MID, SEX, PHENO)
#
# impute = "avg"  : replace NA with per-SNP column mean (same as plink2R)
# impute = "none" : leave NA as-is
# -----------------------------------------------------------------------

read_plink = function(prefix, impute = "none") {
  bim = read.table(paste0(prefix, ".bim"), header = FALSE,
                   colClasses = c("integer","character","numeric",
                                  "integer","character","character"))
  fam = read.table(paste0(prefix, ".fam"), header = FALSE,
                   colClasses = c("character","character","character",
                                  "character","integer","numeric"))
  n   = nrow(fam)
  p   = nrow(bim)
  bps = ceiling(n / 4L)   # bytes per SNP

  con   = file(paste0(prefix, ".bed"), "rb")
  magic = readBin(con, raw(), 3L)
  if (!identical(magic, as.raw(c(0x6c, 0x1b, 0x01))))
    stop("Not a SNP-major PLINK .bed file: ", prefix, ".bed")
  raw = readBin(con, integer(), p * bps, size = 1L, signed = FALSE)
  close(con)

  # Reshape: bps rows (byte index within SNP) × p columns (SNP index)
  raw = matrix(raw, nrow = bps, ncol = p)

  # PLINK 2-bit genotype encoding (2 LSBs = individual 1, next 2 = individual 2 ...):
  #   00 → 2  (hom A1/A1 — minor homozygous; 2 copies of A1)
  #   01 → NA (missing)
  #   10 → 1  (heterozygous; 1 copy of A1)
  #   11 → 0  (hom A2/A2 — major homozygous; 0 copies of A1)
  # This makes mean(column) = 2 * MAF, consistent with plink2R and FUSION.
  lut = c(2, NA_real_, 1, 0)   # index by (2-bit value + 1)

  # Unpack all 4 bit-pairs per byte in 4 vectorised passes.
  # Pass s handles individual indices seq(s+1, bps*4, 4) across all SNPs.
  bed = matrix(NA_real_, nrow = bps * 4L, ncol = p)
  for (s in 0:3) {
    rows       = seq(s + 1L, bps * 4L, by = 4L)
    bits       = bitwAnd(bitwShiftR(raw, s * 2L), 3L)   # integer matrix bps × p
    bed[rows, ] = matrix(lut[bits + 1L], nrow = bps, ncol = p)
  }

  bed = bed[seq_len(n), , drop = FALSE]   # trim padding individuals

  if (impute == "avg") {
    mu      = colMeans(bed, na.rm = TRUE)
    na_snps = which(colSums(is.na(bed)) > 0L)
    for (j in na_snps)
      bed[is.na(bed[, j]), j] = mu[j]
  }

  colnames(bed) = bim[, 2]
  list(bed = bed, bim = bim, fam = fam)
}
