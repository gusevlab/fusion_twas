# --- make_score.R
# Inputs a FUSION wgt.RDat file and prints to screen a plink format score file for individual-level prediction
# ---

arg = commandArgs(trailingOnly=T)
if ( length(arg) == 0 || !file.exists(arg[1]) ) {
cat("\nUsage: Rscript make_score.R [wgt.RDat file] > [SCORE_FILE]\n")
cat("The SCORE_FILE can be loaded directly into the PLINK software (https://www.cog-genomics.org/plink2)\n")
cat("The best performing model (most significant cross-validation p-value) is always selected, even if it's top-snp\n")
cat("To predict into individual-level data use the following command:\n")
cat("plink --bfile [GENOTYPES] --score [SCORE_FILE] 1 2 4\n")

if ( length(arg) > 0 ) {
cat("ERROR:",arg[1],"does not exist or cannot be read\n")
}

q()
}

load( arg[1] )
best = which.min(cv.performance[2,])

if ( names(best) == "lasso" || names(best) == "enet" ) {
keep = wgt.matrix[,best] != 0
} else if ( names(best) == "top1" ) {
keep = which.max(wgt.matrix[,best]^2)
} else { 
keep = 1:nrow(wgt.matrix)
}
write.table( format(cbind( (snps[,c(2,5,6)]) , wgt.matrix[,best])[keep,],digits=3) , quote=F , row.names=F , col.names=F , sep='\t' )
