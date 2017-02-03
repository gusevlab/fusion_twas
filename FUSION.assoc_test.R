suppressMessages(library('plink2R'))
suppressMessages(library("optparse"))

option_list = list(
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to summary statistics (must have SNP and Z column headers) [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--weights", action="store", default=NA, type='character',
              help="File listing molecular weight RDat files (must have columns WGT,ID,CHR,P0,P1) [required]"),
  make_option("--weights_dir", action="store", default=NA, type='character',
              help="Path to directory where weight files (WGT column) are stored [required]"),
  make_option("--ref_ld_chr", action="store", default=NA, type='character',
              help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
  make_option("--perm", action="store", default=0, type='integer',
              help="Maximum number of permutations to perform for each feature test [default: 0/off]"),
  make_option("--perm_minp", action="store", default=0.05, type='double',
              help="Minimum p-value for which to initiate permutation test, if --perm flag present [default: %default]"),
  make_option("--indep_weights", action="store_true", default=FALSE,
              help="Attempt to identify features with independently significant weights [EXPERIMENTAL]."),
  make_option("--cond_weights", action="store_true", default=FALSE,
              help="Attempt to identify features with independently significant weights [EXPERIMENTAL]."),       		  
  make_option("--chr", action="store", default=NA, type='character',
              help="Chromosome to analyze, currently only single chromosome analyses are performed [required]")					  
)

opt = parse_args(OptionParser(option_list=option_list))

allele.qc = function(a1,a2,ref1,ref2) {
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
	snp[["flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)

	return(snp)
}

# Load in summary stats
sumstat = read.table(opt$sumstats,head=T,as.is=T)

# Load in list of weights
# TODO : TEST FOR NO HEADER HERE
wgtlist = read.table(opt$weights,head=T,as.is=T)
wgtlist = wgtlist[ as.character(wgtlist$CHR) == as.character(opt$chr) , ]
chr = unique(wgtlist$CHR)

N = nrow(wgtlist)
out.tbl = data.frame( "FILE" = character(N) , "ID" = character(N) , "CHR" = numeric(N) , "P0" = numeric(N) , "P1" = numeric(N) ,"HSQ" = numeric(N) , "BEST.GWAS.ID" = character(N) , "BEST.GWAS.Z" = numeric(N) , "EQTL.ID" = character(N) , "EQTL.R2" = numeric(N) , "EQTL.Z" = numeric(N) , "EQTL.GWAS.Z" = numeric(N) , "NSNP" = numeric(N) , "NWGT" = numeric(N) , "MODEL" = character(N) , "MODELCV.R2" = numeric(N) , "MODELCV.PV" = numeric(N) , "TWAS.Z" = numeric(N) , "TWAS.P" = numeric(N) , stringsAsFactors=FALSE )

if ( !is.na(opt$perm) && opt$perm > 0 ) {
	out.tbl$PERM.PV = numeric(N)
	out.tbl$PERM.N = numeric(N)
	out.tbl$PERM.ANL_PV = numeric(N)
	permz = qnorm(opt$perm_minp/2,lower.tail=F)
}

# Load in reference data
genos = read_plink(paste(opt$ref_ld_chr,chr,sep=''),impute="avg")

# Match summary data to input, record NA where summary data is missing
m = match( genos$bim[,2] , sumstat$SNP )
sum.missing = is.na(m)
sumstat = sumstat[m,]
sumstat$SNP = genos$bim[,2]
sumstat$A1[ sum.missing ] = genos$bim[sum.missing,5]
sumstat$A2[ sum.missing ] = genos$bim[sum.missing,6]

# QC / allele-flip the input and output
qc = allele.qc( sumstat$A1 , sumstat$A2 , genos$bim[,5] , genos$bim[,6] )

# Flip Z-scores for mismatching alleles
sumstat$Z[ qc$flip ] = -1 * sumstat$Z[ qc$flip ]
sumstat$A1[ qc$flip ] = genos$bim[qc$flip,5]
sumstat$A2[ qc$flip ] = genos$bim[qc$flip,6]

# Remove strand ambiguous SNPs (if any)
if ( sum(!qc$keep) > 0 ) {
	genos$bim = genos$bim[qc$keep,]
	genos$bed = genos$bed[,qc$keep]
	sumstat = sumstat[qc$keep,]
}

# WARNING if too many NAs in summary stats

## For each wgt file:
for ( w in 1:nrow(wgtlist) ) {
	#cat( unlist(wgtlist[w,]) , '\n' )
	# Load weights
	wgt.file = paste(opt$weights_dir,"/",wgtlist$WGT[w],sep='')
	load(wgt.file)
	# Remove NAs (these should not be here)
	wgt.matrix[is.na(wgt.matrix)] = 0
	
	# Match up the SNPs and weights
	m = match( snps[,2] , genos$bim[,2] )
	m.keep = !is.na(m)
	snps = snps[m.keep,]
	wgt.matrix = wgt.matrix[m.keep,]
	cur.genos = scale(genos$bed[,m[m.keep]])
	cur.bim = genos$bim[m[m.keep],]
	# Flip WEIGHTS for mismatching alleles
	qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
	wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]
	rm(snps)

	cur.FAIL = FALSE

	# Match up the SNPs and the summary stats
	m = match(cur.bim[,2] , sumstat$SNP)
	cur.Z = sumstat$Z[m]

	# Identify the best model
	mod.best = (which.max(cv.performance[1,]))
	if ( names(mod.best) == "top1" ) {
		# cat( "WARNING: top eQTL is the best predictor for this gene, continuing with 2nd-best model\n" )
		mod.best = names( which.max(cv.performance[1,colnames(cv.performance)!="top1"]) )
		mod.best = which( colnames(cv.performance) == mod.best )
	}

	# Compute LD matrix
	if ( length(cur.Z) == 0 ) {
		cat( "WARNING : " , unlist(wgtlist[w,]) , " had no overlapping SNPs\n")
		cur.FAIL = TRUE
		out.tbl$NSNP[w] = NA
	} else {
		cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)	
		out.tbl$NSNP[w] = nrow(cur.LD)
		cur.miss = is.na(cur.Z)
		# Impute missing Z-scores
		if ( sum(cur.miss) != 0 ) {
		if ( sum(!cur.miss) == 0 ) {
		cat( "WARNING : " , unlist(wgtlist[w,]) , " had no overlapping GWAS Z-scores\n")
		cur.FAIL = TRUE
		} else {
		cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
		cur.impz = cur.wgt %*% cur.Z[!cur.miss]
		cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
		cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)
		}
		}

	if ( !cur.FAIL ) {
		# Compute TWAS Z-score
		cur.twasz = wgt.matrix[,mod.best] %*% cur.Z
		cur.twasr2pred = wgt.matrix[,mod.best] %*% cur.LD %*% wgt.matrix[,mod.best]
		
		if ( opt$indep_weights ) {
			# the lasso (sparsest) model is always used
			mod.lasso = which(colnames(wgt.matrix) == "lasso")
			mod.assoc = which(colnames(wgt.matrix) == "top1")
			# do informed LD-pruning
			keep = wgt.matrix[,mod.lasso] != 0
			cur.ld = t(cur.genos[,keep]) %*% cur.genos[,keep] / (nrow(cur.genos)-1)
			pruned = rep(F,ncol(cur.ld))
			# walk down the p-value list
			for ( i in order(wgt.matrix[keep,mod.assoc]^2,decreasing=T) ) {
				if ( !pruned[i] ) {
					# remove anything in LD
					pruned[ cur.ld[ i , ]^2 > 0.1 ] = T
					pruned[i] = F
				}
			}
			if ( sum(!pruned) > 2 ) {
				pruned = which( keep )[ !pruned ]			
				# compute indep TWAS Z-stat
				indep.twasz = wgt.matrix[pruned,mod.lasso] %*% cur.Z[pruned,drop=F]
				indep.twasr2pred = wgt.matrix[pruned,mod.lasso] %*% cur.LD[pruned,pruned] %*% wgt.matrix[pruned,mod.lasso]
				# compute SMR test here
				smrstat = (wgt.matrix[pruned,mod.assoc]^2 * cur.Z[pruned,drop=F]^2) / (wgt.matrix[pruned,mod.assoc]^2 + cur.Z[pruned,drop=F]^2)
				cat( wgt.file, paste(format(smrstat,digits=3),sep=',',collapse=',') , paste(format(wgt.matrix[pruned,mod.assoc],digits=3),sep=',',collapse=',') , paste(format(cur.Z[pruned,drop=F],digits=3),sep=',',collapse=',') , cor( wgt.matrix[pruned,mod.assoc]  , cur.Z[pruned,drop=F] ) , cur.twasz/sqrt(cur.twasr2pred) , indep.twasz/sqrt(indep.twasr2pred) , '\n' )
			}
		} else if ( opt$cond_weights ) {
			mod.assoc = which(colnames(wgt.matrix) == "top1")
			# do conditional analysis to get independent effects
			marg.z = wgt.matrix[,mod.assoc]
			cond.z = wgt.matrix[,mod.assoc]
			cond.z.save = rep(0,length(cond.z))
			wgt.keep = rep(F,length(cond.z))
			zthresh = qnorm( 0.05 / length(cond.z) / 2,lower.tail=F)
			
			while ( sum(cond.z^2 > zthresh^2) != 0 ) {
				# add most conditionally significant feature 
				wgt.keep[ which.max(cond.z^2) ] = T
				cond.z.save[ which.max(cond.z^2) ] = cond.z[ which.max(cond.z^2) ]
				cur.dinv = solve(cur.LD[wgt.keep,wgt.keep])
				prev.max = max(cond.z^2)
				for ( i in 1:length(cond.z) ) {
					if ( wgt.keep[i] ) {
						cond.z[i] = 0
					} else if ( max(cur.LD[i,wgt.keep]^2) > 0.5 ) {
						cond.z[i] = 0
					} else {
						# estimate conditional effect size
						cur.b = marg.z[i] - cur.LD[i,wgt.keep,drop=F] %*% (cur.dinv %*% marg.z[wgt.keep,drop=F])
						cur.b.var = ( 1 - cur.LD[i,wgt.keep,drop=F] %*% ( cur.dinv ) %*% t(cur.LD[i,wgt.keep,drop=F]) )
						if ( cur.b.var < 0 ) cond.z[i] = 0
						else cond.z[i] = cur.b / sqrt( cur.b.var )
					}
				}
				
				# TODO: Better quantifying + reporting of overfit here.
				if ( max(cond.z^2) > prev.max ) break()
			}
			
			if ( sum(wgt.keep) > 2 ) {
				tst = cor.test( marg.z[ wgt.keep ] , cur.Z[ wgt.keep ] )
				cat( wgt.file , sum(wgt.keep) , paste( round(marg.z[wgt.keep],1) , collapse=',' ) , paste( round(cur.Z[wgt.keep],1) , collapse=',' ) , tst$p.value , '\n' )
			}
		}		
				
		if ( cur.twasr2pred > 0 ) {
			cur.twas = cur.twasz / sqrt(cur.twasr2pred)
			# Perform the permutation test
			if ( !is.na(opt$perm) && opt$perm > 0 && cur.twas^2 > permz^2 ) {
				perm.twas = rep(NA,opt$perm)
				perm.pval = NA
				for ( i in 1:opt$perm ) {
					perm.wgt = wgt.matrix[ sample( nrow(wgt.matrix) ) , mod.best ]
					perm.twas[i] = perm.wgt %*% cur.Z / sqrt( perm.wgt %*% cur.LD %*% perm.wgt )
				
					# adaptive permutation, stop if 10 instances were observed
					# see: Che et al. PMC4070098 for derivations
					if ( sum(perm.twas[1:i]^2 > cur.twas[1]^2) > 10 ) {
						perm.pval = (sum(perm.twas[1:i]^2 > cur.twas[1]^2) + 1) / (i+1)
						perm.N = i+1
						break()
					}
				}
				
				if ( is.na(perm.pval) ) {
					perm.pval = sum(perm.twas^2 > cur.twas[1]^2) / opt$perm
					perm.N = opt$perm
					# also estimate analytical p-value based on mu+sd of the null			
					anal.zscore = ( cur.twas[1] - 0 ) / sd( perm.twas , na.rm=T )
				} else {
					anal.zscore = NA
				}
			
				out.tbl$PERM.PV[w] = perm.pval
				out.tbl$PERM.N[w] = perm.N
				out.tbl$PERM.ANL_PV[w] = 2*pnorm(abs(anal.zscore),lower.tail=F)
			}
		} else {
			cur.FAIL=T
		}
	}
	}

	# populate the output
	out.tbl$FILE[w] = wgt.file
	out.tbl$CHR[w] = wgtlist$CHR[w]
	out.tbl$P0[w] = wgtlist$P0[w]
	out.tbl$P1[w] = wgtlist$P1[w]
	out.tbl$ID[w] = wgtlist$ID[w]
	out.tbl$HSQ[w] = hsq[1]
	out.tbl$MODEL[w] = colnames( cv.performance )[ mod.best ]
	out.tbl$MODELCV.R2[w] = cv.performance[1,mod.best]
	out.tbl$MODELCV.PV[w] = cv.performance[2,mod.best]

	eqtlmod = colnames(wgt.matrix) == "top1"
	if ( cur.FAIL || length(eqtlmod) == 0 ) {
		out.tbl$EQTL.ID[w] = NA
		out.tbl$EQTL.R2[w] = NA
		out.tbl$EQTL.Z[w] =  NA
		out.tbl$EQTL.GWAS.Z[w] = NA
	} else {
		topeqtl = which.max( wgt.matrix[,eqtlmod]^2 )
		out.tbl$EQTL.ID[w] = names( topeqtl )
		out.tbl$EQTL.R2[w] = cv.performance[1,eqtlmod]
		out.tbl$EQTL.Z[w] = wgt.matrix[ topeqtl , eqtlmod ]
		out.tbl$EQTL.GWAS.Z[w] = cur.Z[ topeqtl ]
	}
	
	topgwas = which.max( cur.Z^2 )
	if ( !cur.FAIL && length(topgwas) != 0 && !is.na(topgwas) ) {
		out.tbl$BEST.GWAS.ID[w] = rownames(wgt.matrix)[ topgwas ]
		out.tbl$BEST.GWAS.Z[w] = cur.Z[ topgwas ]
	} else {
		out.tbl$BEST.GWAS.ID[w] = NA
		out.tbl$BEST.GWAS.Z[w] = NA
	}
	
	if ( !cur.FAIL ) {
		out.tbl$NWGT[w] = sum( wgt.matrix[,mod.best] != 0 )
		out.tbl$TWAS.Z[w] = cur.twas
	} else {
		out.tbl$TWAS.Z[w] = NA
	}
}

# compute p-value
out.tbl$TWAS.P = 2*(pnorm( abs(out.tbl$TWAS.Z) , lower.tail=F))

# WRITE MHC TO SEPARATE FILE
mhc = out.tbl$CHR == 6 & out.tbl$P0 > 26e6 & out.tbl$P1 < 34e6
if ( sum( mhc ) > 0 ) {
	write.table( format( out.tbl[mhc,] , digigs=3 ) , quote=F , row.names=F , sep='\t' , file=paste(opt$out,".MHC",sep='') )
}
write.table( format( out.tbl[!mhc,] , digits=3 ) , quote=F , row.names=F , sep='\t' , file=opt$out )
