# ---
# Code to compute pairwise model correlation across different sets of models (e.g. TWAS/RWAS/CWAS)
# Sample command:
# 	Rscript pairs.R --pos1 TCGA-BRCA.GE.TUMOR.pos --pos2 TCGA-BRCA.GE.NORMAL.pos --chr 1 --ref_ld_chr ../../LDREF/1000G.EUR. --window 100000
# ---

suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))
suppressMessages(library("RColorBrewer"))

allele.qc = function(a1,a2,ref1,ref2) {
	a1 = toupper(a1)
	a2 = toupper(a2)
	ref1 = toupper(ref1)
	ref2 = toupper(ref2)

	ref = ref1
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip1 = flip

	ref = ref2
	flip = ref
	flip[ref == "A"] = "T"
	flip[ref == "T"] = "A"
	flip[ref == "G"] = "C"
	flip[ref == "C"] = "G"
	flip2 = flip;

	snp = list()
	snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
        snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
        snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F

	snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

	return(snp)
}

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.dir <- dirname(script.name)

# Perform permutation test
# opt$perm = FALSE

option_list = list(
  make_option("--pos1", action="store", default=NA, type='character',
              help="Path to list of models [required]"),
  make_option("--pos2", action="store", default=NA, type='character',
              help="Path to list of models to compare to [required]"),
  make_option("--ref_ld_chr", action="store", default=NA, type='character',
              help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
  make_option("--chr", action="store", default=NA, type='character',
              help="Chromosome to analyze, currently only single chromosome analyses are performed [required]"),
  make_option("--window", action="store", default=500e3, type='integer',
              help="Window around pos1 features to take [default: %default]")                   
)

opt = parse_args(OptionParser(option_list=option_list))
options( digits = 3 )

chr = opt$chr
# read all weights in list and compute union of SNPs to analyze
pos1 = read.table(opt$pos1,as.is=T,head=T)
pos2 = read.table(opt$pos2,as.is=T,head=T)

pos1 = pos1[ pos1$CHR == chr , ]
pos2 = pos2[ pos2$CHR == chr , ]

# load in genotype files by chromosome, restrict to matching SNPs and combine
genos = read_plink(paste(opt$ref_ld_chr,chr,sep=''),impute="avg")
MAFS = apply(genos$bed,2,mean)
genos$bed = scale(genos$bed)
N = nrow(genos$fam)

# load in all models in pos1
ge1_g.matrix = matrix(nrow=nrow(genos$bed),ncol=nrow(pos1))

for ( i in 1:nrow(pos1) ) {
	load( pos1$WGT[i] )
	wgt.matrix[is.na(wgt.matrix)] = 0
	# Match up the SNPs and weights
	m = match( snps[,2] , genos$bim[,2] )
	m.keep = !is.na(m)
	snps = snps[m.keep,]
	wgt.matrix = wgt.matrix[m.keep,]
	
	cur.genos = genos$bed[,m[m.keep]]
	cur.bim = genos$bim[m[m.keep],]
	# Flip WEIGHTS for mismatching alleles
	qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
	wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]
	
	# Predict into reference
	mod = which.min(cv.performance[2,])
	ge1_g.matrix[,i] = scale(cur.genos %*% wgt.matrix[ , mod ])
	if ( i %% 100 == 0 ) cat( "Loading pos1",i,"/",nrow(pos1),'\n' , file=stderr() )
}

# load in all models in pos2
ge2_g.matrix = matrix(nrow=nrow(genos$bed),ncol=nrow(pos2))

for ( i in 1:nrow(pos2) ) {
	load( pos2$WGT[i] )
	wgt.matrix[is.na(wgt.matrix)] = 0
	# Match up the SNPs and weights
	m = match( snps[,2] , genos$bim[,2] )
	m.keep = !is.na(m)
	snps = snps[m.keep,]
	wgt.matrix = wgt.matrix[m.keep,]
	
	cur.genos = genos$bed[,m[m.keep]]
	cur.bim = genos$bim[m[m.keep],]
	# Flip WEIGHTS for mismatching alleles
	qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
	wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]
	
	# Predict into reference
	mod = which.min(cv.performance[2,])
	ge2_g.matrix[,i] = scale( cur.genos %*% wgt.matrix[ , mod ] )
	if ( i %% 100 == 0 ) cat( "Loading pos2",i,"/",nrow(pos2),'\n' , file=stderr() )
}

# make comparisons
for ( i in 1:nrow(pos1) ) {
	# Get nearby models
	cur.keep = pos1$CHR[i] == pos2$CHR & pos2$P0 < pos1$P1[i] + opt$window & pos2$P1 > pos1$P0[i] - opt$window
	for ( ii in which(cur.keep) ) {
		dist = mean( c(pos1$P0[i],pos1$P1[i] ) ) - mean( c(pos2$P0[ii],pos2$P1[ii]) )
		cat( pos1$WGT[i] , pos2$WGT[ii] , dist , cor(ge1_g.matrix[,i],ge2_g.matrix[,ii]) , '\n' , sep='\t' )
	}
}	
