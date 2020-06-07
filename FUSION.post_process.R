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
  make_option("--input", action="store", default=NA, type='character',
              help="Path to TWAS test output [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--sumstats", action="store", default=NA, type='character',
              help="Path to LDSC format summary statistics [required]"),
  make_option("--ref_ld_chr", action="store", default=NA, type='character',
              help="Prefix to reference LD files in binary PLINK format by chromosome [required]"),
  make_option("--minp_input", action="store", default=1.0, type='double',
              help="Minimum p-value to include feature in analysis [default: %default]"),
  make_option("--max_r2", action="store", default=0.90, type='double',
              help="Features with r^2 greater than this will be considered identical [default: %default]"),		
  make_option("--min_r2", action="store", default=0.05, type='double',
              help="Features with r^2 less than this will be considered independent [default: %default]"),
  make_option("--locus_win", action="store", default=100e3, type='integer',
              help="How much to expand each feature (in bp) to define contiguous loci [default: %default]"), 
  make_option("--max_cz_increase", action="store", default=1.96, type='double',
              help="Maximum allowed increase in conditional Z-score (can indicate LD mismatch / complex locus) [default: %default]"),                     
  make_option("--plot", action="store_true", default=FALSE,
              help="Generate pdf plots for each locus [default: OFF]"),
  make_option("--plot_legend", action="store", default=NA, type='character',
              help="Add a legend to the plot to color code reference panels [options: all/joint for which genes to include]"),
  make_option("--plot_corr", action="store_true", default=FALSE,
              help="Plot correlation of genetic values for each locus (requires corrplot library) [default: OFF]"),  
  make_option("--plot_individual", action="store_true", default=FALSE,
              help="Plot conditional analyses of individual genes [default: OFF]"),  
  make_option("--plot_eqtl", action="store_true", default=FALSE,
              help="Plot eQTL signal below GWAS signal (requires --plot; --plot_individual when multiple genes associated) [default: OFF]"),
  make_option("--plot_scatter", action="store_true", default=FALSE,
              help="Plot TWAS scatterplot (requires --plot; --plot_individual when multiple genes associated) [default: OFF]"),
  make_option("--report", action="store_true", default=FALSE,
              help="Generate a report document with information on each locus [default: OFF]"),              
  make_option("--omnibus", action="store_true", default=FALSE,
              help="Perform the omnibus test for genes (ID field) with multiple models (FILE field). NOTE: This disables all other tests."),
  make_option("--omnibus_corr", action="store", default=NA, type='character',
              help="Only print the pairwise correlations between reference panels for the specified model [options: top1,blup,bslmm,enet,lasso or best]"),
  make_option("--eqtl_model", action="store", default="top1", type='character',
              help="Name of the predictive for which weights should be used for marginal eQTL plotting (experimental) [default: %default]"),	
  make_option("--ldsc", action="store_true", default=FALSE,
              help="Compute LD-scores across all features. NOTE: This disables all other tests."),             
  make_option("--save_loci", action="store_true", default=FALSE,
              help="Save conditioned GWAS results for each locus [default: %default]"),
  make_option("--chr", action="store", default=NA, type='character',
              help="Chromosome to analyze, currently only single chromosome analyses are performed [required]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),  
  make_option("--zthresh", action="store", default=FALSE, type='double',
              help="Z-score threshold for performing summary conditional analysis [default: %default]")                            
)

opt = parse_args(OptionParser(option_list=option_list))
options( digits = 3 )

# --- TODO
# Better computation of clumps
# Perform feature selection only within each clump and print separate outputs
# ---

chr = opt$chr
# read all weights in list and compute union of SNPs to analyze
wgtlist = read.table(opt$input,as.is=T,head=T)
if( !("TWAS.P" %in% colnames(wgtlist)) ) {
	cat( "ERROR: --input does not contain TWAS.P column header\n" , file=stderr() )
	q()
}

wgtlist = wgtlist[ wgtlist$CHR == chr & wgtlist$TWAS.P < opt$minp_input & !is.na(wgtlist$TWAS.P) , ]

# load in genotype files by chromosome, restrict to matching SNPs and combine
genos = read_plink(paste(opt$ref_ld_chr,chr,sep=''),impute="avg")
MAFS = apply(genos$bed,2,mean)
genos$bed = scale(genos$bed)
N = nrow(genos$fam)

if ( opt$plot ) {
	if ( opt$plot_eqtl && opt$plot_scatter ) {
		cat( "WARNING: both --plot_eqtl and --plot_scatter cannot be enabled, plotting eQTL only\n" , file=stderr() )
		opt$plot_scatter = FALSE
	}
	
	if ( !file.exists("glist-hg19") ) {
		if ( file.exists( paste(script.dir,"/glist-hg19",sep='') ) ){
		glist = read.table(paste(script.dir,"/glist-hg19",sep=''),as.is=T)
		glist = glist[glist[,1] == chr,]
		} else {
		cat( "WARNING: glist-hg19 file listing gene locations (with header: CHR P0 P1 ID) needed for locus plots\n" , file=stderr() )
		glist = matrix(nrow=0,ncol=4)
		}
	} else {
		glist = read.table("glist-hg19",as.is=T)
		glist = glist[glist[,1] == chr,]
	}
	colnames(glist) = c("CHR","P0","P1","ID")	
} else if ( (opt$plot_individual || !is.na(opt$plot_legend) || opt$plot_eqtl) ) {
	cat( "WARNING: plotting flags set without --plot, figures will not be generated\n" , file=stderr() )
}

# --- OMNIBUS TEST ACROSS MATCHING FEATURES
if ( opt$omnibus ) {
	# identify duplicate features
	dup.genes = unique(wgtlist$ID[ duplicated( wgtlist$ID ) ])

	cat( "GENE\tREF1\tREF2\tCORR\n" , file=paste( opt$out , ".CORR.dat" , sep=''))
	if ( is.na(opt$omnibus_corr) ) cat( "GENE\tNUM.REF\tMIN.TWAS.P\tNUM.REF.PRUNED\tOMNIBUS.P\n" , file=paste(opt$out , ".omnibus.pv" , sep=''))
				
	for ( g in dup.genes ) {
		cur.keep = which( wgtlist$ID == g )
		M = length( cur.keep )
		ge_g.matrix = matrix(nrow=nrow(genos$bed),ncol=M)
		ref.names = unlist(lapply(strsplit(dirname(wgtlist$FILE[ cur.keep ]),"/"),tail,1))
		colnames(ge_g.matrix) = ref.names
		cur.drop = rep(F,length(cur.keep))
		
		for ( i in 1:length( cur.keep ) ) {
			load( wgtlist$FILE[ cur.keep[i] ] )
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
	
			if ( !is.na(opt$omnibus_corr) && opt$omnibus_corr != "best" ) {
				mod = which(colnames(wgt.matrix) == opt$omnibus_corr)
			} else {
				mod = which(colnames(wgt.matrix) == wgtlist$MODEL[ cur.keep[i] ])			
			}
			
			# Predict into reference
			if ( length(mod) == 0 ) cur.drop[i] = T
			else ge_g.matrix[,i] = cur.genos %*% wgt.matrix[ , mod ]
		}
		
		cur.keep = cur.keep[ !cur.drop ]
		ge_g.matrix = ge_g.matrix[,!cur.drop]

		if ( sum(cur.keep) > 1 ) {
			ge_g.cor = cor(ge_g.matrix)
			cur.wgt = wgtlist[ cur.keep , ]
			ge_z  = cur.wgt$TWAS.Z

			if ( !is.na(opt$omnibus_corr) ) {
				# just print correlations
				# format the table for output:
				cur.tbl = ge_g.cor
				cur.tbl[lower.tri(cur.tbl,diag=T)] = NA
				cur.tbl = as.data.frame(as.table(cur.tbl))
				cur.tbl = cur.tbl[!is.na(cur.tbl[,3]),]
				write.table( cbind(g,format(cur.tbl,digits=3)) , quote =F , col.names=F , row.names=F , sep='\t' , file=paste( opt$out , ".CORR.dat" , sep='') , append=T )	
			} else {
				# perform the omnibus test 
				# do "informed" LD-pruning to remove highly correlated genes
				pruned = rep(F,M)
				# walk down the p-value list
				for ( i in order(ge_z^2,decreasing=T) ) {
					if ( !pruned[i] ) {
						# remove anything in LD
						pruned[ ge_g.cor[ i , ]^2 > opt$max_r2 ] = T
						pruned[i] = F
					}
				}
				
				if ( sum(!pruned) > 1 ) {
					chisq = t(ge_z[!pruned]) %*% solve(ge_g.cor[!pruned,!pruned]) %*% ge_z[!pruned]
					pv.chi = pchisq( chisq , df=sum(!pruned) , lower.tail=F)
					cur.tbl = ge_g.cor[!pruned,!pruned]
					cur.tbl[lower.tri(cur.tbl,diag=T)] = NA
					cur.tbl = as.data.frame(as.table(cur.tbl))
					cur.tbl = cur.tbl[!is.na(cur.tbl[,3]),]
					write.table( cbind(g,format(cur.tbl,digits=3)) , quote =F , col.names=F , row.names=F , sep='\t' , file=paste( opt$out , ".CORR.dat" , sep='') , append=T )
					cat( g , M , min( cur.wgt$TWAS.P ) , sum(!pruned) , pv.chi , '\n' , sep='\t' , file=paste(opt$out , ".omnibus.pv" , sep='') , append=T )
				}
			}
		}
	}
	q()
# --- DONE OMNIBUS TEST
}

# list of SNPs that overlap loaded features
genos.keep = rep(F,nrow(genos$bim))
# matrix for predicted expression:
ge_g.matrix = matrix(nrow=nrow(genos$bed),ncol=nrow(wgtlist))

# matrix for eQTL
if ( opt$plot_eqtl ) {
	eqtl.z = list()
	eqtl.pos = list()
}

for ( i in 1:nrow(wgtlist) ) {
	load( wgtlist$FILE[i] )
	wgt.matrix[is.na(wgt.matrix)] = 0
	# Match up the SNPs and weights
	m = match( snps[,2] , genos$bim[,2] )
	m.keep = !is.na(m)
	snps = snps[m.keep,]
	wgt.matrix = wgt.matrix[m.keep,]
	genos.keep[ m[m.keep] ] = T
	
	cur.genos = genos$bed[,m[m.keep]]
	cur.bim = genos$bim[m[m.keep],]
	# Flip WEIGHTS for mismatching alleles
	qc = allele.qc( snps[,5] , snps[,6] , cur.bim[,5] , cur.bim[,6] )
	wgt.matrix[qc$flip,] = -1 * wgt.matrix[qc$flip,]
	if ( opt$plot_eqtl ) {
		eqtl.pos[[i]] = snps[,4]
		eqtl.z[[i]] = wgt.matrix[ , which(colnames(wgt.matrix) == opt$eqtl_model) ]
	}
	
	# Predict into reference
	mod = which(colnames(wgt.matrix) == wgtlist$MODEL[i])
	ge_g.matrix[,i] = cur.genos %*% wgt.matrix[ , mod ]
}
ge_g.matrix = scale( ge_g.matrix )

if ( opt$ldsc ) {
	LDSC = t(genos$bed) %*% ge_g.matrix / (N-1)
	LDSC = LDSC^2 - 1/N
	LDSC = apply(LDSC,1,sum)
	f.LDSC <- gzfile(paste(opt$out,".l2.ldscore.gz",sep=''), "w")
	LDSC.df = data.frame( "CHR" = genos$bim[,1] , SNP = genos$bim[,2] , BP = genos$bim[,4] , CM = genos$bim[,3] , "MAF" = MAFS , "LD" = LDSC )
	write.table( format(LDSC.df,digits=3) , quote=F , row.names=F , col.names=T , file=f.LDSC )
	close(f.LDSC)
	cat( ncol(ge_g.matrix) , '\n' , file=paste(opt$out,".l2.M_5_50",sep='') , sep='')
	q()
}

# compute gene LD matrix
ge_g.ld = cor(ge_g.matrix)
ge_g.z = wgtlist$TWAS.Z

if( opt$zthresh ) {
	zthresh = opt$zthresh
} else {
	zthresh = qnorm( 0.05 / nrow( wgtlist ) / 2,lower.tail=F)
}

if( opt$verbose > 1 ) cat( nrow( wgtlist ) , " weights considered, a weight must have Z^2 > ", zthresh^2 , " to be retained in the model\n" , sep='' , file=stderr() ) 

if ( opt$plot_corr ) {
	rownames(ge_g.ld) = paste(wgtlist$PANEL,wgtlist$ID)
	colnames(ge_g.ld) = rep("_",ncol(ge_g.ld))
	library("corrplot")

	pca = eigen(ge_g.ld)
	cat( "pca varexp : " , pca$values / sum(pca$values) , '\n' )

	sz = nrow(ge_g.ld) / 4
	if ( sz < 3 ) sz = 3
	
	pdf( file=paste(opt$out,".corrplot.pdf",sep=''),height=sz,width=sz)

	par(cex = 0.3)
	corrplot( ge_g.ld^2 , method="color" , type="upper" , order="hclust" , addCoef.col="black" , tl.cex=0.5/0.3 , col=colorRampPalette(c("blue","white","red"))(200) ,  cl.lim=c(0,1))
	par(cex = 1 )
	dev.off()
}
	
# zero out low LD
ge_g.ld[ ge_g.ld^2 < opt$min_r2 ] = 0


# --- PERFORM FEATURE SELECTION ACROSS ALL GENES
cond.z = ge_g.z
ge.keep = rep(F,length(ge_g.z))
ge.drop = rep(F,length(ge_g.z))

if ( sum(cond.z^2 > zthresh^2) == 0 ) {
	if( opt$verbose > 0 ) cat( "WARNING: no models had an absolute marginal association statistic higher than " , zthresh , ". Skipping\n" , sep='' , file=stderr() ) 
} else {
	while ( sum(cond.z^2 > zthresh^2) != 0 ) {
		# add most conditionally significant feature 
		ge.keep[ which.max(cond.z^2) ] = T
		if( opt$verbose > 1 ) cat( wgtlist$FILE[ which.max(cond.z^2) ] , " added to model with conditional Z-score of ", cond.z[which.max(cond.z^2)] , "\n" , sep='' , file=stderr() ) 

		cur.dinv = solve(ge_g.ld[ge.keep,ge.keep])
		for ( i in 1:length(cond.z) ) {
			if ( ge.keep[i] || ge.drop[i] ) {
				cond.z[i] = 0
			} else if ( max(ge_g.ld[i,ge.keep]^2) > opt$max_r2 ) {
				cond.z[i] = 0
				if( opt$verbose > 1 ) cat( wgtlist$FILE[ i ] , " dropped from the model due to correlation higher than --max_r2 to other genes in the model\n" , sep='' , file=stderr() ) 
			} else {
				# estimate conditional effect size
				cur.b = ge_g.z[i] - ge_g.ld[i,ge.keep,drop=F] %*% (cur.dinv %*% ge_g.z[ge.keep,drop=F])
				cur.b.var = ( 1 - ge_g.ld[i,ge.keep,drop=F] %*% ( cur.dinv ) %*% t(ge_g.ld[i,ge.keep,drop=F]) )
				if ( cur.b.var < 0 ) {
					cond.z[i] = 0
					if( opt$verbose > 1 ) cat( wgtlist$FILE[ i ] , " dropped from the model due to misspecified correlation with other genes\n" , sep='' , file=stderr() ) 
				} else cond.z[i] = cur.b / sqrt( cur.b.var )
			}
		}
	
		# drop any genes for which the conditional association increased over the marginal
		cur.unstable = cond.z^2 > ge_g.z^2 + (opt$max_cz_increase)^2
		if ( sum(cur.unstable) != 0 ) {
			if( opt$verbose > 1 ) cat( unlist( wgtlist$FILE[ cur.unstable ] ) , " became more significant after conditional analysis and are dropped from the model (this is usually a sign of LD mismatch or a complex locus)\n" , sep='' , file=stderr() ) 
			ge.drop[ cur.unstable ] = T
			cond.z[ cur.unstable ] = 0
		}
	}
	if( opt$verbose > 1 ) cat( "final best conditional Z^2 = " , max(cond.z^2) , "\n" , sep='' , file=stderr() ) 

	# FINAL FEATURE SELECTED ESTIMATES:
	ge.keep = ge.keep & !ge.drop
	joint.keep = ge.keep
	cur.dinv = solve(ge_g.ld[joint.keep,joint.keep])

	# joint estimate for features kept in the model:
	joint.b = ( cur.dinv %*% ge_g.z[joint.keep,drop=F] )
	joint.se = sqrt( diag( cur.dinv ) )
	joint.z = joint.b / joint.se
	joint.p = 2*(pnorm( abs( joint.z ) , lower.tail=F))
	df.out = data.frame( "FILE" = wgtlist$FILE[ge.keep] , "ID" = wgtlist$ID[ge.keep] , "TWAS.Z" = wgtlist$TWAS.Z[ge.keep] , "TWAS.P" = wgtlist$TWAS.P[ge.keep]  , "JOINT.BETA" = joint.b , "JOINT.BETA.SE" = joint.se , "JOINT.Z" = joint.z , "JOINT.P" = joint.p )
	write.table( format(df.out,digits=2) , quote=F , col.names=T , row.names=F , sep='\t' , file=paste(opt$out,".joint_included.dat",sep='') )

	# conditional estimate for features dropped from the model:
	cond.b = ge_g.z[!ge.keep] - ge_g.ld[!ge.keep,ge.keep,drop=F] %*% (cur.dinv %*% ge_g.z[ge.keep])
	cond.b.var = diag( 1 - ge_g.ld[!ge.keep,ge.keep,drop=F] %*% ( cur.dinv ) %*% t(ge_g.ld[!ge.keep,ge.keep,drop=F]) )
	cond.b.var[ cond.b.var < 0 ] = NA
	cond.se = sqrt( cond.b.var )
	cond.z = cond.b / cond.se
	cond.p = 2*(pnorm( abs( cond.z ) , lower.tail=F))
	df.out = data.frame( "FILE" = wgtlist$FILE[!ge.keep] , "ID" = wgtlist$ID[!ge.keep] , "TWAS.Z" = wgtlist$TWAS.Z[!ge.keep] , "TWAS.P" = wgtlist$TWAS.P[!ge.keep] , "COND.BETA" = cond.b , "COND.BETA.SE" = cond.se , "COND.Z" = cond.z , "COND.P" = cond.p )
	write.table( format(df.out,digits=2) , quote=F , col.names=T , row.names=F , sep='\t' , file=paste(opt$out,".joint_dropped.dat",sep='') )
	# --- DONE FEATURE SELECTION
}

# --- PERFORM CONDITIONAL GWAS
# load in summary statistics and flip appropriately
sumstat = read.table(opt$sumstats,head=T,as.is=T)

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

# Identify contiguous loci
runs = rle(genos.keep)
if ( runs$val[1] ) {
	loc.starts = c(1,cumsum(runs$lengths)[ !runs$val ]+1)
	loc.ends = cumsum(runs$lengths)[ runs$val ]
	loc.starts = loc.starts[ 1:length(loc.ends) ]
} else {
	loc.starts = cumsum(runs$lengths)[ !runs$val ]+1
	loc.ends = cumsum(runs$lengths)[ runs$val ]
	loc.starts = loc.starts[ 1:length(loc.ends) ]
}
loc.starts = genos$bim[loc.starts,4] - opt$locus_win
loc.ends = genos$bim[loc.ends,4] + opt$locus_win
if( opt$verbose > 0 ) cat( length(loc.starts) , " strictly non-overlapping loci\n" , sep='' , file=stderr() )

# Consolidate overlapping loci
cons.loc.starts = loc.starts[1]
cons.loc.ends = loc.ends[1]
loc.ctr = 1

if ( length(loc.starts) > 1 ) {
for ( i in 2:length(loc.starts) ) {
	if ( loc.starts[i] < cons.loc.ends[ loc.ctr ] ) {
		cons.loc.ends[ loc.ctr ] = max( cons.loc.ends[ loc.ctr ] , loc.ends[i] )
	} else {
		cons.loc.starts = c(cons.loc.starts,loc.starts[i])
		cons.loc.ends = c(cons.loc.ends,loc.ends[i])
		loc.ctr = loc.ctr+1
	}
}
}

if( opt$verbose > 0 ) cat( "consolidated to ", length(cons.loc.starts) , " non-overlapping loci with ", opt$locus_win , " bp buffer\n" , sep='' , file=stderr() )

if ( opt$report ) {
	file.report = paste(opt$out,".report",sep='')
	cat( "FILE" , "CHR" , "P0" , "P1" , "HIT.GENES" , "JOINT.GENES" , "BEST.TWAS.P" , "BEST.SNP.P" , "COND.SNP.P" , "VAR.EXP\n" , sep='\t' , file=file.report )
}

# iterate over loci
for ( i in 1:length(cons.loc.starts) ) {

	# get overlapping features
	ge.keep = wgtlist$P0 < cons.loc.ends[i] & wgtlist$P1 > cons.loc.starts[i]
	
	# --- PERFORM FEATURE SELECTION IN THIS CLUMP
	marg.z = ge_g.z
	cond.z = ge_g.z
	cond.z[ !ge.keep ] = 0
	joint.keep = rep(F,length(ge_g.z))

	if ( sum(cond.z^2 > zthresh^2) == 0 ) {
		if( opt$verbose > 0 ) cat( "WARNING: no models in CLUMP ", i , " had an absolute marginal association statistic higher than " , zthresh , ". Skipping\n" , sep='' , file=stderr() ) 
	} else {
		while ( sum(cond.z^2 > zthresh^2) != 0 ) {
			# add most conditionally significant feature 
			joint.keep[ which.max(cond.z^2) ] = T
			if( opt$verbose > 1 ) cat( (wgtlist$FILE[ge.keep])[ which.max(cond.z^2) ] , " added to model with conditional Z-score of ", cond.z[which.max(cond.z^2)] , "\n" , sep='' , file=stderr() ) 

			cur.dinv = solve(ge_g.ld[joint.keep,joint.keep])
			for ( ii in which(cond.z != 0) ) {
				if ( joint.keep[ii] ) {
					cond.z[ii] = 0
				} else if ( max(ge_g.ld[ii,joint.keep]^2) > opt$max_r2 ) {
					cond.z[ii] = 0
				} else {
					# estimate conditional effect size
					cur.b = marg.z[ii] - ge_g.ld[ii,joint.keep,drop=F] %*% (cur.dinv %*% marg.z[joint.keep,drop=F])
					cur.b.var = ( 1 - ge_g.ld[ii,joint.keep,drop=F] %*% ( cur.dinv ) %*% t(ge_g.ld[ii,joint.keep,drop=F]) )
					if ( cur.b.var < 0 ) cond.z[ii] = 0
					else cond.z[ii] = cur.b / sqrt( cur.b.var )
				}
			}
		}
		cur.dinv = solve(ge_g.ld[joint.keep,joint.keep])
		# ----

		cur.keep = genos$bim[,4] > cons.loc.starts[i] & genos$bim[,4] < cons.loc.ends[i]
		snp.z = sumstat$Z[ cur.keep ]
		snp.ge.ld = t(genos$bed[,cur.keep]) %*% ge_g.matrix[,joint.keep] / ( N - 1 )

		# TODO (optional) : impute any missing z-scores
		snp.cond.b = snp.z - snp.ge.ld %*% (cur.dinv %*% ge_g.z[joint.keep])
		snp.cond.se = diag( sqrt( 1 - snp.ge.ld %*% ( cur.dinv ) %*% t(snp.ge.ld) ) )

		# zero out any SNPs with > max_r2 LD to a gene
		autocor = apply(snp.ge.ld^2,1,max,na.rm=T) > opt$max_r2
		snp.cond.b[ autocor ] = 0
		snp.cond.se[ autocor ] = 1
		snp.cond.z = snp.cond.b / snp.cond.se
		pv = 2*(pnorm( abs(snp.z) , lower.tail=F))

		# generate before/after manhattan plot
		if ( opt$plot ) {
			# get overlapping genes
			cur.glist = glist[ apply(glist[,2:3],1,min) < cons.loc.ends[i] & apply(glist[,2:3],1,max) > cons.loc.starts[i] , ]
			# identify and remove genes with same name as features
			m = match( cur.glist[,4] , wgtlist$ID[ ge.keep ] )
			cur.glist = cur.glist[ is.na( m ), ]

			tot = nrow(cur.glist) + sum(ge.keep)
		

			pdf( file=paste(opt$out,".loc_",i,".pdf",sep='') , height = 2 + log(tot)/2 )
			par( oma = c(5,4,0,0)+0.1 , mar=c(0,0,1,1)+0.1 , xpd=NA , las=1 )
		
			lay.eqtl = matrix( c( rep(1,ceiling(log(tot)/3)),2,2,3,3) , ncol=1 )
			lay.gwas = matrix( c( rep(1,ceiling(log(tot)/2)),2,2,2) , ncol=1 )
			lay.scatter = cbind(lay.gwas,lay.gwas,rep(3,nrow(lay.gwas)))
		
			if ( opt$plot_eqtl && sum(ge.keep) == 1 ) {
				layout(lay.eqtl)
			} else if ( opt$plot_scatter && sum(ge.keep) == 1 ) {
				layout(lay.scatter)
			} else {
				layout(lay.gwas)
			}
		
			plot( 0 , 0 , type="n" , ylim = c(-0.1,1.1) , xlim = range( genos$bim[ cur.keep , 4 ] / 1e6 ) , bty="n" , xlab="" , ylab="" , xaxt="n" , yaxt="n" )

			# --- Gene names and positions
	
			# color code for weights reference
			if ( !is.na(opt$plot_legend) ) {
			if ( opt$plot_legend == "all" ) {
				ref.names = unlist(lapply(strsplit(dirname(wgtlist$FILE[ ge.keep ]),"/"),tail,1))
				uni.ref.names = sort( unique(ref.names) )
				# clr.ref = rainbow( length(uni.ref.names) )
				clr.ref = brewer.pal( max(3,length(uni.ref.names)) , "Set3" )
				m = match( ref.names , uni.ref.names )
				clr.leg = c( rep(NA,nrow(cur.glist)) , clr.ref[ m ] )
				clr.num = c( rep("",nrow(cur.glist)) , m )
			} else if ( opt$plot_legend == "joint" ) {
				ref.names = unlist(lapply(strsplit(dirname(wgtlist$FILE[ joint.keep & ge.keep ]),"/"),tail,1))
				uni.ref.names = sort( unique(ref.names) )
				#clr.ref = rainbow( length(uni.ref.names) )
				clr.ref = brewer.pal( max(3,length(uni.ref.names)) , "Set3" )

				m = match( ref.names , uni.ref.names )
				clr.leg = rep(NA , sum(ge.keep))
				clr.leg[ joint.keep[ge.keep] ] = clr.ref[ m ]
		
				clr.num = rep("",sum(ge.keep))			
				clr.num[ joint.keep[ge.keep] ] = m
		
				clr.leg = c( rep(NA,nrow(cur.glist)) , clr.leg )
				clr.num = c( rep("",nrow(cur.glist)) , clr.num )		
			}
			}

			clr.fg = rep("gray30",nrow(cur.glist))
			clr.bg = rep("gray",nrow(cur.glist))
			clr.fg2 = rep("#3182bd", sum(ge.keep) )
			clr.bg2 = rep("#deebf7", sum(ge.keep) )
			clr.bg2[ joint.keep[ ge.keep ] ] = "#66bd63"
			clr.fg2[ joint.keep[ ge.keep ] ] = "#1a9850"
			clr.fg = c( clr.fg , clr.fg2 )
			clr.bg = c( clr.bg , clr.bg2 )	
				
			cur.glist = rbind( cur.glist , wgtlist[ ge.keep , c("CHR","P0","P1","ID") ] )	
		
			ord = order(apply(cur.glist[,2:3],1,min))
			cur.glist = cur.glist[ ord , ]	
			clr.fg = clr.fg[ ord ]
			clr.bg = clr.bg[ ord ]
			if ( !is.na( opt$plot_legend) ) {
				clr.leg = clr.leg[ ord ]
				clr.num = clr.num[ ord ]
			}
	
			cur.gstart = apply(cur.glist[,2:3],1,min)/1e6
			cur.gend = apply(cur.glist[,2:3],1,max)/1e6
			# compute size with text
			g.size = strwidth(cur.glist[,4]) * 1.5
			g.ysize = strheight(cur.glist[1,4])
			g.ypos = rep(0,nrow(cur.glist))
			# iterate over each gene and drop row if it overlaps
			for ( g in 1:nrow(cur.glist) ) {
				cur.row = 1
				while ( sum(cur.gstart[g] - g.size < cur.gend & cur.gend[g] > cur.gstart - g.size & g.ypos == cur.row) != 0 ) cur.row = cur.row + 1
				g.ypos[g] = cur.row
			}
	
			# compute text scaling
			txt.scale = min( 1 , ( 1 / (1+max(g.ypos)) ) / g.ysize )
			# rescale ypos to between 0 and 1
			if ( 1 / max(g.ypos) < g.ysize ) {
				g.ypos = (g.ypos - min(g.ypos)) / (max(g.ypos) - min(g.ypos))
			} else {
				g.ypos = (g.ypos - 1) * g.ysize * 1.5
			}
			# --- done with gene names and positions
					
			g.ctr = 1
			cur.gene = which(ge.keep)[ g.ctr ]		
			while ( TRUE ) {
				# --- GENE PLOTTING
				if ( g.ctr > 1 ) {
					if ( opt$plot_eqtl ) {
						layout(lay.eqtl)
					} else if ( opt$plot_scatter ) {
						layout(lay.scatter)
					} else {		
						layout(lay.gwas)
					}	
					plot( 0 , 0 , type="n" , ylim = c(-0.1,1.1) , xlim = range( genos$bim[ cur.keep , 4 ] / 1e6 ) , bty="n" , xlab="" , ylab="" , xaxt="n" , yaxt="n" )
					clr.fg = rep("gray30",nrow(cur.glist))
					clr.bg = rep("gray",nrow(cur.glist))
					clr.bg[ cur.glist$ID == wgtlist$ID[ cur.gene ] ] = "#66bd63"
					clr.fg[ cur.glist$ID == wgtlist$ID[ cur.gene ] ] = "#1a9850"		
				}
			
				rect( cur.gstart , g.ypos - txt.scale*g.ysize*0.45 , cur.gend , g.ypos + txt.scale*g.ysize*0.45 , border=clr.fg , col=clr.bg )
				text( cur.gstart , g.ypos , cur.glist[,4] , pos=2 , cex=txt.scale , col=clr.fg )
				if ( !is.na( opt$plot_legend) ) {
					# legend dots
					points( cur.gend + strwidth('*') , g.ypos , pch=19 , col=clr.leg , cex=txt.scale )
					text( cur.gend + strwidth('*') , g.ypos , clr.num , pch=19 , cex=txt.scale * 0.5 )
					legend( "topleft" , legend=paste(1:(length(uni.ref.names)),uni.ref.names), pch=19 , cex=txt.scale , col=clr.ref , bty="n" )
				}
				# --- DONE GENE PLOTTING

				# Manhattan plots:
				if ( opt$plot_eqtl && ( sum(ge.keep) == 1 || g.ctr > 1 ) ) {
					plot( genos$bim[ cur.keep , 4 ] / 1e6 , -log10(pv) , pch=19 , cex=0.5 , col="gray" , xlab="" , ylab="-log10(P-value)" , bty="n" )				
				} else {			
					plot( genos$bim[ cur.keep , 4 ] / 1e6 , -log10(pv) , pch=19 , cex=0.5 , col="gray" , xlab=paste("chr ", chr , " physical position (MB)",sep='') , ylab="-log10(P-value)" , bty="n" )
				}
				if ( g.ctr == 1 ) {
					pv.cond = 2*(pnorm( abs(snp.cond.z) , lower.tail=F))
				} else {
					pv.cond = 2*(pnorm( abs(cur.snp.cond.z) , lower.tail=F))
				}
				points( genos$bim[ cur.keep , 4 ] / 1e6 , -log10(pv.cond) , pch=19 , cex=0.5 , col="#3182bd" )
			
				# EQTL plots
				if ( sum(ge.keep) == 1 || g.ctr > 1 ) {
					if ( opt$plot_scatter ) {
						plot( snp.ge.ld , snp.z  , pch=19 , cex=0.5 , xlab="Corr. to TWAS" , ylab="GWAS Z-score" , xlim=c(-1,1) , ylim=c(-1*max(abs(snp.z)),max(abs(snp.z))) )	
						par(xpd=F)
						abline( 0 , ge_g.z[cur.gene] , lty=2 )
						par(xpd=NA)					
					} else if ( opt$plot_eqtl ) {
						pv.eqtl = 2*(pnorm( abs(eqtl.z[[ cur.gene ]] ) , lower.tail=F ))
						plot( eqtl.pos[[cur.gene]] / 1e6 , -log10(pv.eqtl) , xlim=range( genos$bim[ cur.keep , 4 ] / 1e6 ) , pch=19 , cex=0.5 , col="#3182bd" , xlab=paste("chr ", chr , " physical position (MB)",sep='') , ylab="-log10(P-value)" , bty="n" )
					} 
				}
			
				# Now plot each individual gene or exit
				if ( sum(ge.keep) > 1 && opt$plot_individual && g.ctr <= sum(ge.keep) ) {
					cur.gene = which(ge.keep)[ g.ctr ]
					snp.ge.ld = t(genos$bed[,cur.keep]) %*% ge_g.matrix[,cur.gene] / ( N - 1 )

					snp.cond.b = snp.z - snp.ge.ld %*% ge_g.z[cur.gene]
					snp.cond.se = diag( sqrt( 1 - snp.ge.ld %*% t(snp.ge.ld) ) )

					# zero out any SNPs with > max_r2 LD to a gene
					autocor = apply(snp.ge.ld^2,1,max,na.rm=T) > opt$max_r2
					snp.cond.b[ autocor ] = NA
					snp.cond.se[ autocor ] = NA
				
					# save conditional z-score for this gene separately from main conditional result
					cur.snp.cond.z = snp.cond.b / snp.cond.se
					g.ctr = g.ctr + 1
				} else {
					break()
				}
			}
			dev.off()
		}

		pv.cond = 2*(pnorm( abs(snp.cond.z) , lower.tail=F))
	
		# print conditional Z-scores	
		if ( opt$save_loci ) {
			df.out = data.frame( "SNP" = genos$bim[cur.keep,2] , "POS" = genos$bim[cur.keep,4] , "GWAS.Z" = snp.z , "GWAS.P" = pv , "GWAS_cond.Z" = snp.cond.z , "GWAS_cond.P" = pv.cond )
			write.table( format(df.out,digits=2) , quote=F , col.names=T , row.names=F , sep='\t' , file=paste(opt$out,".loc_",i,".cond",sep='') )
		}
		if ( opt$report ) {
			df.out = data.frame( "SNP" = genos$bim[cur.keep,2] , "POS" = genos$bim[cur.keep,4] , "GWAS.LOGP" = -log10(pv) , "GWAS_cond.LOGP" = -log10(pv.cond) )
			write.table( format(df.out,digits=2,trim=T) , quote=F , col.names=T , row.names=F , sep=',' , file=paste(opt$out,".loc_",i,".cond.csv",sep='') )
		}
	
		cat( "locus " , i , " best GWAS Chisq\t" , max(snp.z^2,na.rm=T) , '\n' , sep='')
		cat( "locus " , i , " best GWAS Chisq conditioned\t" , snp.cond.z[which.max(snp.z^2)]^2 , '\n' , sep='' )
		cat( "locus " , i , " best conditioned Chisq\t" , max(snp.cond.z^2,na.rm=T) , '\n' , sep='' )
	
		if ( opt$report ) {
			best.snp.chisq = max(snp.z^2,na.rm=T)
			cond.snp.chisq = snp.cond.z[which.max(snp.z^2)]^2
			cat( paste(opt$out,".loc_",i,sep='') , opt$chr , range( genos$bim[ cur.keep , 4 ] ) , length(unique(wgtlist$ID[ge.keep])) , sum( joint.keep[ ge.keep ] ) , min(wgtlist$TWAS.P[ ge.keep ],na.rm=T) , 2*pnorm(sqrt(best.snp.chisq), lower.tail=F) , 2*pnorm(sqrt(cond.snp.chisq), lower.tail=F) , 1 - cond.snp.chisq / best.snp.chisq , '\n' , sep='\t' , file=file.report , append=T )
		
			cur.wgt = wgtlist[ ge.keep , ]
			cur.wgt$JOINT = joint.keep[ ge.keep ]
		
			# compute correlation of each gene to the top SNP
			cur.wgt$TOP.SNP.COR = round(t(t( (genos$bed[,cur.keep])[,which.max(snp.z^2),drop=F] ) %*% ge_g.matrix[,ge.keep] / ( N - 1 )),2)
			write.table( cur.wgt , quote=F , row.names=F , col.names=T , sep='\t' , file=paste(opt$out,".loc_",i,".genes",sep='')  )
		}
	}
}
# --- DONE CONDITIONAL ANALYSIS
