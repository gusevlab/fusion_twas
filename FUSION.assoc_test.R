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
  make_option("--force_model", action="store", default=NA, type='character',
              help="Force specific predictive model to be used, no flag (default) means select most significant cross-val. Options: blup,lasso,top1,enet"),
  make_option("--caviar", action="store_true", default=FALSE,
              help="Generate eCAVIAR-format (Z,LD) files for fine-mapping [default off]"),
  make_option("--jlim", action="store_true", default=FALSE,
              help="NOT IMPLEMENTED: Compute JLIM statistic [Chun et al Nat Genet 2017].\nRequires jlimR library installed. [default: %default]"),			  
  make_option("--max_impute", action="store", default=0.5 , type='double',
              help="Maximum fraction of SNPs allowed to be missing per gene (will be imputed using LD). [default: %default]"),			  
  make_option("--min_r2pred", action="store", default=0.7 , type='double',
              help="Minimum average LD-based imputation accuracy allowed for expression weight SNP Z-scores. [default: %default]"),			  
  make_option("--perm", action="store", default=0, type='integer',
              help="Maximum number of permutations to perform for each feature test [default: 0/off]"),
  make_option("--perm_minp", action="store", default=0.05, type='double',
              help="Minimum p-value for which to initiate permutation test, if --perm flag present [default: %default]"),    		  
  make_option("--chr", action="store", default=NA, type='character',
              help="Chromosome to analyze, currently only single chromosome analyses are performed [required]"),
  make_option("--coloc_P", action="store", default=NA, type='double',
              help="P-value below which to compute COLOC statistic [Giambartolomei et al PLoS Genet 2013]\nRequires coloc library installed and --GWASN flag. [default NA/off]"),
  make_option("--GWASN", action="store", default=NA, type='integer',
              help="Total GWAS/sumstats sample size for inference of standard GWAS effect size."),
  make_option("--PANELN", action="store", default=NA, type='character',
              help="File listing sample size for each panel for inference of standard QTL effect size, cross-referenced against 'PANEL' column in weights file")      
)

opt = parse_args(OptionParser(option_list=option_list))

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

# Load in summary stats
sumstat = read.table(opt$sumstats,head=T,as.is=T)

# Load in list of weights
# TODO : TEST FOR NO HEADER HERE
wgtlist = read.table(opt$weights,head=T,as.is=T)
wgtlist = wgtlist[ as.character(wgtlist$CHR) == as.character(opt$chr) , ]
chr = unique(wgtlist$CHR)

N = nrow(wgtlist)
out.tbl = data.frame( "PANEL" = rep(NA,N) , "FILE" = character(N) , "ID" = character(N) , "CHR" = numeric(N) , "P0" = character(N) , "P1" = character(N) ,"HSQ" = numeric(N) , "BEST.GWAS.ID" = character(N) , "BEST.GWAS.Z" = numeric(N) , "EQTL.ID" = character(N) , "EQTL.R2" = numeric(N) , "EQTL.Z" = numeric(N) , "EQTL.GWAS.Z" = numeric(N) , "NSNP" = numeric(N) , "NWGT" = numeric(N) , "MODEL" = character(N) , "MODELCV.R2" = character(N) , "MODELCV.PV" = character(N) , "TWAS.Z" = numeric(N) , "TWAS.P" = numeric(N) , stringsAsFactors=FALSE )

if ( opt$jlim ) {
	suppressMessages(library('jlimR'))
	jlim.r2res = 0.8
	jlim.thresholdingP = 0.1
	out.tbl$JLIM.P = numeric(N)
	out.tbl$JLIM.STAT = numeric(N)
}

if ( !is.na(opt$coloc_P) ) {
	if ( is.na(opt$GWASN) || opt$GWASN < 1 ) {
		cat("ERROR : --GWASN flag required to be positive integer for COLOC analysis\n")
		q()
	}
	if ( sum(names(wgtlist) == "N") == 0 ) {
		if ( sum(names(wgtlist) == "PANEL") == 0 || is.na(opt$PANELN) ) {
			cat("ERROR : 'N' field needed in weights file or 'PANEL' column and --PANELN flag required for COLOC analysis\n")
			q()
		} else { 
			paneln = read.table(opt$PANELN,as.is=T,head=T,sep='\t')
			m = match( wgtlist$PANEL , paneln$PANEL )
			wgtlist$N = paneln$N[ m ]
		}
	}
	suppressMessages(library('coloc'))
	out.tbl$COLOC.PP0 = as.numeric(rep(NA,N))
	out.tbl$COLOC.PP1 = as.numeric(rep(NA,N))
	out.tbl$COLOC.PP2 = as.numeric(rep(NA,N))
	out.tbl$COLOC.PP3 = as.numeric(rep(NA,N))
	out.tbl$COLOC.PP4 = as.numeric(rep(NA,N))
}

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

# TODO: WARNING if too many NAs in summary stats

FAIL.ctr = 0

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
	wgt.matrix = wgt.matrix[m.keep,,drop=F]
	cur.genos = scale(genos$bed[,m[m.keep]])
	cur.bim = genos$bim[m[m.keep],]
	# Flip WEIGHTS for mismatching alleles
	qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
	wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]

	cur.FAIL = FALSE

	# Match up the SNPs and the summary stats
	m = match(cur.bim[,2] , sumstat$SNP)
	cur.Z = sumstat$Z[m]

	# which rows have rsq
	row.rsq = grep( "rsq" , rownames(cv.performance) )
	# which rows have p-values
	row.pval = grep( "pval" , rownames(cv.performance) )	
	
	# Identify the best model
	if ( !is.na(opt$force_model) ) {
		mod.best = which( colnames(wgt.matrix) == opt$force_model )
		if ( length(mod.best) == 0 ) {
			cat( "WARNING : --force_model" , mod.best ,"does not exist for", unlist(wgtlist[w,]) , "\n")
			cur.FAIL = TRUE
		}	
	} else {
		# get the most significant model
		mod.best = which.min(apply(cv.performance[row.pval,,drop=F],2,min,na.rm=T))
	}
	if ( length(mod.best) == 0 ) {
		cat( "WARNING : " , unlist(wgtlist[w,]) , " did not have a predictive model ... skipping entirely\n" )
		FAIL.ctr = FAIL.ctr + 1
		next
	}

	if ( sum(wgt.matrix[, mod.best] != 0) == 0 ) {
		cat( "WARNING : " , unlist(wgtlist[w,]) , names(cv.performance)[ mod.best ] , "had", length(cur.Z) , "overlapping SNPs, but none with non-zero expression weights, try more SNPS or a different model\n")
		cur.FAIL = TRUE
	}

	# if this is a top1 model, clear out all the other weights
	if ( substr( (colnames(cv.performance))[ mod.best ],1,4) == "top1" ) wgt.matrix[ -which.max(wgt.matrix[,mod.best]^2)  , mod.best] = 0

	# Compute LD matrix
	if ( length(cur.Z) == 0 ) {
		cat( "WARNING : " , unlist(wgtlist[w,]) , " had no overlapping SNPs\n")
		cur.FAIL = TRUE
		out.tbl$NSNP[w] = NA
	} else if ( !cur.FAIL ) {
		cur.LD = t(cur.genos) %*% cur.genos / (nrow(cur.genos)-1)	
		out.tbl$NSNP[w] = nrow(cur.LD)
		cur.miss = is.na(cur.Z)
		# Impute missing Z-scores
		if ( sum(cur.miss) != 0 ) {
			if ( sum(!cur.miss) == 0 ) {
				cat( "WARNING : " , unlist(wgtlist[w,]) , "had no overlapping GWAS Z-scores, skipping this gene\n")
				cur.FAIL = TRUE
			} else if ( mean(cur.miss) > opt$max_impute ) {
				cat( "WARNING : " , unlist(wgtlist[w,]) , "had" , sum(cur.miss) , "/" , length(cur.miss) , "non-overlapping GWAS Z-scores, skipping this gene.\n")
				cur.FAIL = TRUE
			} else {
				cur.wgt =  cur.LD[cur.miss,!cur.miss] %*% solve( cur.LD[!cur.miss,!cur.miss] + 0.1 * diag(sum(!cur.miss)) )
				cur.impz = cur.wgt %*% cur.Z[!cur.miss]
				cur.r2pred = diag( cur.wgt %*% cur.LD[!cur.miss,!cur.miss] %*% t(cur.wgt) )
				cur.Z[cur.miss] = cur.impz / sqrt(cur.r2pred)

				all.r2pred = rep(1,length(cur.Z))
				all.r2pred[ cur.miss ] = cur.r2pred
				if ( sum(is.na(all.r2pred)) != 0 ) {
					cat( "WARNING : " , unlist(wgtlist[w,]) , "had missing GWAS Z-scores that could not be imputed, skipping this gene.\n" )
					cur.FAIL = TRUE
				} else if ( mean( all.r2pred[ wgt.matrix[,mod.best] != 0 ] ) < opt$min_r2pred ) {
					cat( "WARNING : " , unlist(wgtlist[w,]) , "had mean GWAS Z-score imputation r2 of" , mean( all.r2pred[ wgt.matrix[,mod.best] != 0 ] ) , "at expression weight SNPs, skipping this gene.\n")
					cur.FAIL = TRUE
				}
			}
		}
		
		if ( !cur.FAIL ) {
			# Compute TWAS Z-score
			cur.twasz = wgt.matrix[,mod.best] %*% cur.Z
			cur.twasr2pred = wgt.matrix[,mod.best] %*% cur.LD %*% wgt.matrix[,mod.best]
					
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
				cat( "WARNING : " , unlist(wgtlist[w,]) , " had zero predictive accuracy, try a different model.\n")
			}
		}
	}

	# populate the output
	if ( sum(names(wgtlist) == "PANEL") == 1 ) out.tbl$PANEL[w] = wgtlist$PANEL[w]
	out.tbl$FILE[w] = wgt.file
	out.tbl$CHR[w] = wgtlist$CHR[w]
	out.tbl$P0[w] = wgtlist$P0[w]
	out.tbl$P1[w] = wgtlist$P1[w]
	out.tbl$ID[w] = wgtlist$ID[w]
	if ( exists("hsq") ) {
		out.tbl$HSQ[w] = hsq[1]
	}
	out.tbl$MODEL[w] = colnames( cv.performance )[ mod.best ]
	out.tbl$MODELCV.R2[w] = paste(format(cv.performance[row.rsq,mod.best],digits=2,trim=T),collapse=',')
	out.tbl$MODELCV.PV[w] = paste(format(cv.performance[row.pval,mod.best],digits=2,trim=T),collapse=',')

	eqtlmod = colnames(wgt.matrix) == "top1"
	topeqtl = which.max( wgt.matrix[,eqtlmod]^2 )
	
	if ( cur.FAIL || sum(eqtlmod) == 0 || length(topeqtl) == 0 || is.na(topeqtl) ) {
		out.tbl$EQTL.ID[w] = NA
		out.tbl$EQTL.R2[w] = NA
		out.tbl$EQTL.Z[w] =  NA
		out.tbl$EQTL.GWAS.Z[w] = NA
	} else {
		out.tbl$EQTL.ID[w] = rownames(wgt.matrix)[topeqtl]
		out.tbl$EQTL.R2[w] = cv.performance[1,eqtlmod]
		out.tbl$EQTL.Z[w] = wgt.matrix[ topeqtl , eqtlmod ]
		out.tbl$EQTL.GWAS.Z[w] = cur.Z[ topeqtl ]
		
		# write CAVIAR inputs
		if( opt$caviar ) {
			cur.Z = as.matrix(cur.Z,ncol=1)
			rownames(cur.Z) = snps[,2]
			cav.out = paste( opt$out , wgtlist$ID[w] , "CAVIAR" , sep='.' )
			write.table( format(cur.LD,digits=3) , quote=F , col.names=F , row.names=F , file = paste( cav.out , ".LD" , sep='' ) )
			write.table( format(wgt.matrix[,eqtlmod],digits=3) , quote=F , col.names=F , sep='\t' , file = paste( cav.out , ".EQTL.Z" , sep='') )
			write.table( format(cur.Z,digits=3) , quote=F , col.names=F , sep='\t' , file = paste( cav.out , ".GWAS.Z" , sep='') )
		}
		
		# perform JLIM test
		if ( opt$jlim ) {
			jlim.lambda.t = jlimR:::calc.stat( cur.Z , wgt.matrix[,eqtlmod] , cur.LD , cur.LD , jlim.r2res )
			# TODO : Sample permuted eQTL Z-scores from the MvN using LD matrix
			jlim.nulldist = jlimR:::perm.test( cur.Z , wgt.matrix[,eqtlmod], permmat, cur.LD, cur.LD, jlim.thresholdingP, jlim.r2res, jlim.lambda.t)
			permP = sum(jlim.nulldist >= jlim.lambda.t, na.rm = TRUE)/sum(!is.na(NULLDIST))
			out.tbl$JLIM.STAT[w] = jlim.lambda.t
			out.tbl$JLIM.P[w] = permP		
		}
	}
	
	topgwas = which.max( cur.Z^2 )
	if ( !cur.FAIL && length(topgwas) != 0 && !is.na(topgwas) ) {
		out.tbl$BEST.GWAS.ID[w] = snps[ topgwas , 2 ]
		out.tbl$BEST.GWAS.Z[w] = cur.Z[ topgwas ]
	} else {
		out.tbl$BEST.GWAS.ID[w] = NA
		out.tbl$BEST.GWAS.Z[w] = NA
	}
	
	if ( !cur.FAIL ) {
		out.tbl$NWGT[w] = sum( wgt.matrix[,mod.best] != 0 )
		out.tbl$TWAS.Z[w] = cur.twas
		out.tbl$TWAS.P[w] = 2*(pnorm( abs(out.tbl$TWAS.Z[w]) , lower.tail=F))
	} else {
		out.tbl$TWAS.Z[w] = NA
		out.tbl$TWAS.P[w] = NA
	}

	# perform COLOC test
	if ( !is.na(opt$coloc_P) && !is.na(out.tbl$TWAS.Z[w]) && out.tbl$TWAS.P[w] < opt$coloc_P && !is.na(wgtlist$N[w]) ) {
		b1 = wgt.matrix[,eqtlmod] / sqrt(wgtlist$N[w])
		b2 = cur.Z / sqrt(opt$GWASN)

		vb1 = rep(1/wgtlist$N[w],length(b1))
		vb2 = rep(1/opt$GWASN,length(b2))

		err = suppressMessages(capture.output(clc <- coloc.abf(dataset1=list(beta=b1,varbeta=vb1,type="quant",N=wgtlist$N[w],sdY=1),dataset2=list(beta=b2,varbeta=vb2,type="quant",N=opt$GWASN,sdY=1))))
		out.tbl$COLOC.PP0[w] = round(clc$summary[2],3)
		out.tbl$COLOC.PP1[w] = round(clc$summary[3],3)
		out.tbl$COLOC.PP2[w] = round(clc$summary[4],3)
		out.tbl$COLOC.PP3[w] = round(clc$summary[5],3)
		out.tbl$COLOC.PP4[w] = round(clc$summary[6],3)
	}
	if ( cur.FAIL ) FAIL.ctr = FAIL.ctr + 1
}

cat("Analysis completed.\n")
cat("NOTE:",FAIL.ctr,"/",nrow(wgtlist),"genes were skipped\n")
if ( FAIL.ctr / nrow(wgtlist) > 0.1 ) {
cat("If a large number of genes were skipped, verify that your GWAS Z-scores, expression weights, and LDREF data use the same SNPs (or nearly)\n")
cat("Or consider pre-imputing your summary statistics to the LDREF markers using summary-imputation software such as [https://github.com/bogdanlab/fizi]\n")
}
# compute p-value
#out.tbl$TWAS.P = 2*(pnorm( abs(out.tbl$TWAS.Z) , lower.tail=F))

# WRITE MHC TO SEPARATE FILE
mhc = as.numeric(out.tbl$CHR) == 6 & as.numeric(out.tbl$P0) > 26e6 & as.numeric(out.tbl$P1) < 34e6

out.tbl$P0 = apply( as.matrix(out.tbl$P0) , 1 , toString )
out.tbl$P1 = apply( as.matrix(out.tbl$P1) , 1 , toString )

if ( sum( mhc ) > 0 ) {
	cat("Results in the MHC are written to",paste(opt$out,".MHC",sep=''),", evaluate with caution due to complex LD structure\n")
	write.table( format( out.tbl[mhc,] , digits=3 ) , quote=F , row.names=F , sep='\t' , file=paste(opt$out,".MHC",sep='') )
}
write.table( format( out.tbl[!mhc,] , digits=3 ) , quote=F , row.names=F , sep='\t' , file=opt$out )
