# ==== TODO
# * Make sure BLUP/BSLMM weights are being scaled properly based on MAF

suppressMessages(library("optparse"))
suppressMessages(library('plink2R'))
suppressMessages(library('glmnet'))
suppressMessages(library('methods'))

option_list = list(
  make_option("--bfile", action="store", default=NA, type='character',
              help="Path to PLINK binary input file prefix (minus bed/bim/fam) [required]"),
  make_option("--out", action="store", default=NA, type='character',
              help="Path to output files [required]"),
  make_option("--tmp", action="store", default=NA, type='character',
              help="Path to temporary files [required]"),
  make_option("--pheno", action="store", default=NA, type='character',
              help="Path to molecular phenotype file (PLINK format) [optional, taken from bfile otherwise]"),
  make_option("--PATH_plink", action="store", default="plink", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gcta", action="store", default="gcta_nr_robust", type='character',
              help="Path to plink executable [%default]"),
  make_option("--PATH_gemma", action="store", default="gemma", type='character',
              help="Path to plink executable [%default]"),
  make_option("--covar", action="store", default=NA, type='character',
              help="Path to quantitative covariates (PLINK format) [optional]"),
  make_option("--resid", action="store_true", default=FALSE,
              help="Also regress the covariates out of the genotypes [default: %default]"),              
  make_option("--hsq_p", action="store", default=0.01, type='double',
              help="Minimum heritability p-value for which to compute weights [default: %default]"),
  make_option("--hsq_set", action="store", default=NA, type='double',
              help="Skip heritability estimation and set hsq estimate to this value [optional]"),
  make_option("--crossval", action="store", default=5, type='double',
              help="How many folds of cross-validation, 0 to skip [default: %default]"),
  make_option("--verbose", action="store", default=1, type="integer",
              help="How much chatter to print: 0=nothing; 1=minimal; 2=all [default: %default]"),
  make_option("--noclean", action="store_true", default=FALSE,
              help="Do not delete any temporary files (for debugging) [default: %default]"),
  make_option("--rn", action="store_true", default=FALSE,
              help="Rank-normalize the phenotype after all QC: [default: %default]"),
  make_option("--save_hsq", action="store_true", default=FALSE,
              help="Save heritability results even if weights are not computed [default: %default]"),			  
  make_option("--models", action="store", default="blup,lasso,top1,enet", type='character',
              help="Comma-separated list of prediction models [default: %default]\n
					Available models:\n
					top1:\tTop eQTL (standard marginal eQTL Z-scores always computed and stored)\n
					blup:\t Best Unbiased Linear Predictor (dual of ridge regression)\n
					bslmm:\t Bayesian Sparse Linear Model (spike/slab MCMC)\n
					lasso:\t LASSO regression (with heritability used as lambda)\n
					enet:\t Elastic-net regression (with mixing parameter of 0.5)\n")			  
		  
)

opt = parse_args(OptionParser(option_list=option_list))
models = unique( c(unlist(strsplit(opt$models,',')),"top1") )
M = length(models)

if ( ! all(models %in% c('blup', 'lasso', 'top1', 'enet', 'bslmm')) ) {
	cat( "ERROR: --models flag included invalid models\n" , sep='' , file=stderr() )
	q()
}

if ( opt$verbose == 2 ) {
  SYS_PRINT = F
} else {
  SYS_PRINT = T
}

# --- PREDICTION MODELS

# BSLMM
weights.bslmm = function( input , bv_type , snp , out=NA ) {
	if ( is.na(out) ) out = paste(input,".BSLMM",sep='')

	arg = paste( opt$PATH_gemma , " -miss 1 -maf 0 -r2 1 -rpace 1000 -wpace 1000 -bfile " , input , " -bslmm " , bv_type , " -o " , out , sep='' )
	system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)
	eff = read.table( paste(out,".param.txt",sep=''),head=T,as.is=T)
	eff.wgt = rep(NA,length(snp))
	m = match( snp , eff$rs )
	m.keep = !is.na(m)
	m = m[m.keep]
	eff.wgt[m.keep] = (eff$alpha + eff$beta * eff$gamma)[m]
	return( eff.wgt )
}

# PLINK: LASSO
weights.lasso = function( input , hsq , snp , out=NA ) {
	if ( is.na(out) ) out = paste(input,".LASSO",sep='')

	arg = paste( opt$PATH_plink , " --allow-no-sex --bfile " , input , " --lasso " , hsq , " --out " , out , sep='' )
	system( arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT )
	if ( !file.exists(paste(out,".lasso",sep='')) ) {
	cat( paste(out,".lasso",sep='') , " LASSO output did not exist\n" )
	eff.wgt = rep(NA,length(snp))
	} else {
	eff = read.table( paste(out,".lasso",sep=''),head=T,as.is=T)
	eff.wgt = rep(0,length(snp))
	m = match( snp , eff$SNP )
	m.keep = !is.na(m)	
	m = m[m.keep]
	eff.wgt[m.keep] = eff$EFFECT[m]
	}
	return( eff.wgt )
}

# Marginal Z-scores (used for top1)
weights.marginal = function( genos , pheno , beta=F ) {
	if ( beta ) eff.wgt = t( genos ) %*% (pheno) / ( nrow(pheno) - 1)
	else eff.wgt = t( genos ) %*% (pheno) / sqrt( nrow(pheno) - 1 )
	return( eff.wgt )
}

# Elastic Net
weights.enet = function( genos , pheno , alpha=0.5 ) {
	eff.wgt = matrix( 0 , ncol=1 , nrow=ncol(genos) )
	# remove monomorphics
	sds = apply( genos  , 2 , sd )
	keep = sds != 0 & !is.na(sds)
	enet = cv.glmnet( x=genos[,keep] , y=pheno , alpha=alpha , nfold=5 , intercept=T , standardize=F )
	eff.wgt[ keep ] = coef( enet , s = "lambda.min")[2:(sum(keep)+1)]
	return( eff.wgt )
}

# --- CLEANUP
cleanup = function() {
	if ( ! opt$noclean ) {
		arg = paste("rm -f " , opt$tmp , "*", sep='')
		system(arg)
	}
}

# Perform i/o checks here:
files = paste(opt$bfile,c(".bed",".bim",".fam"),sep='')
if ( !is.na(opt$pheno) ) files = c(files,opt$pheno)
if ( !is.na(opt$covar) ) files = c(files,opt$covar)

for ( f in files ) {
	if ( !file.exists(f) ){
		cat( "ERROR: ", f , " input file does not exist\n" , sep='', file=stderr() )
		cleanup()
		q()
	}
}

if ( system( paste(opt$PATH_plink,"--help") , ignore.stdout=T,ignore.stderr=T ) != 0 ) {
	cat( "ERROR: plink could not be executed, set with --PATH_plink\n" , sep='', file=stderr() )
	cleanup()
	q()
}

if ( !is.na(opt$hsq_set) && system( opt$PATH_gcta , ignore.stdout=T,ignore.stderr=T ) != 0 ){
	cat( "ERROR: gcta could not be executed, set with --PATH_gcta\n" , sep='', file=stderr() )
	cleanup()
	q()
}

if ( sum(models=="bslmm" | models=="blup") != 0 && system( paste(opt$PATH_gemma,"-h") , ignore.stdout=T,ignore.stderr=T ) != 0 ){
	cat( "ERROR: gemma could not be executed, set with --PATH_gemma or remove 'bslmm' and 'blup' from models\n" , sep='', file=stderr() )
	cleanup()
	q()
}
# ---

fam = read.table(paste(opt$bfile,".fam",sep=''),as.is=T)

# Make/fetch the phenotype file
if ( !is.na(opt$pheno) ) {
	pheno.file = opt$pheno
	pheno = read.table(pheno.file,as.is=T)
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(pheno[,1],pheno[,2]) )
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	m = m[m.keep]
	pheno = pheno[m,]
} else {
	pheno.file = paste(opt$tmp,".pheno",sep='')
	pheno = fam[,c(1,2,6)]
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}

if ( opt$rn ) {
	library('GenABEL')
	library(preprocessCore)
	pheno[,3] = rntransform( pheno[,3] )
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}

# Load in the covariates if needed
if ( !is.na(opt$covar) ) {
	covar = ( read.table(opt$covar,as.is=T,head=T) )
	if ( opt$verbose >= 1 ) cat( "Loaded",ncol(covar)-2,"covariates\n")
	# Match up data
	m = match( paste(fam[,1],fam[,2]) , paste(covar[,1],covar[,2]) )
	m.keep = !is.na(m)
	fam = fam[m.keep,]
	pheno = pheno[m.keep,]
	m = m[m.keep]
	covar = covar[m,]
	reg = summary(lm( pheno[,3] ~ as.matrix(covar[,3:ncol(covar)]) ))
	if ( opt$verbose >= 1 ) cat( reg$r.sq , "variance in phenotype explained by covariates\n" )
	pheno[,3] = scale(reg$resid)
	raw.pheno.file = pheno.file
	pheno.file = paste(pheno.file,".resid",sep='')
	write.table(pheno,quote=F,row.names=F,col.names=F,file=pheno.file)
}

geno.file = opt$tmp
# recode to the intersection of samples and new phenotype
arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$bfile," --pheno ",pheno.file," --keep ",pheno.file," --make-bed --out ",geno.file,sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# --- HERITABILITY ANALYSIS
if ( is.na(opt$hsq_set) ) {
if ( opt$verbose >= 1 ) cat("### Estimating heritability\n")

# 1. generate GRM
arg = paste( opt$PATH_plink," --allow-no-sex --bfile ",opt$tmp," --make-grm-bin --out ",opt$tmp,sep='')
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# 2. estimate heritability
if ( !is.na(opt$covar) ) {
arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",raw.pheno.file," --qcovar ",opt$covar," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
} else {
arg = paste( opt$PATH_gcta ," --grm ",opt$tmp," --pheno ",pheno.file," --out ",opt$tmp," --reml --reml-no-constrain --reml-lrt 1",sep='')
}
system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

# 3. evaluate LRT and V(G)/Vp
if ( !file.exists( paste(opt$tmp,".hsq",sep='') ) ) {
	cat(opt$tmp,"does not exist, likely GCTA could not converge, skipping gene\n",file=stderr())
	cleanup()
	q()
}

hsq.file = read.table(file=paste(opt$tmp,".hsq",sep=''),as.is=T,fill=T)
hsq = as.numeric(unlist(hsq.file[hsq.file[,1] == "V(G)/Vp",2:3]))
hsq.pv = as.numeric(unlist(hsq.file[hsq.file[,1] == "Pval",2]))

if ( opt$verbose >= 1 ) cat("Heritability (se):",hsq,"LRT P-value:",hsq.pv,'\n')
if ( opt$save_hsq ) cat( opt$out , hsq , hsq.pv , '\n' , file=paste(opt$out,".hsq",sep='') )

# 4. stop if insufficient
if ( hsq[1] < 0 || hsq.pv > opt$hsq_p ) {
	cat(opt$tmp," : heritability ",hsq[1],"; LRT P-value ",hsq.pv," : skipping gene\n",sep='',file=stderr())
	cleanup()
	q()
}

} else {
if ( opt$verbose >= 1 ) cat("### Skipping heritability estimate\n")
hsq = opt$hsq_set
hsq.pv = NA
}

# read in genotypes
genos = read_plink(geno.file,impute="avg")
mafs = apply(genos$bed,2,mean)/2
sds = apply(genos$bed,2,sd)
# important : genotypes are standardized and scaled here:
genos$bed = scale(genos$bed)
pheno = genos$fam[,c(1,2,6)]
pheno[,3] = scale(pheno[,3])

# check if any genotypes are NA
nasnps = apply( is.na(genos$bed) , 2 , sum )
if ( sum(nasnps) != 0 ) {
	cat( "WARNING :",sum(nasnps != 0),"SNPs could not be scaled and were zeroed out, make sure all SNPs are polymorphic\n" , file=stderr())
	genos$bed[,nasnps != 0] = 0
}

# regress covariates out of the genotypes as well (this is more accurate but slower)
if ( !is.na(opt$covar) && opt$resid ) {
	if ( opt$verbose >= 1 ) cat("regressing covariates out of the genotypes\n")
	for ( i in 1:ncol(genos$bed) ) {
		genos$bed[,i] = summary(lm( genos$bed[,i] ~ as.matrix(covar[,3:ncol(covar)]) ))$resid
	}
	genos$bed = scale(genos$bed)
}

N.tot = nrow(genos$bed)
if ( opt$verbose >= 1 ) cat(nrow(pheno),"phenotyped samples, ",nrow(genos$bed),"genotyped samples, ",ncol(genos$bed)," markers\n")

# --- CROSSVALIDATION ANALYSES
set.seed(1)
cv.performance = matrix(NA,nrow=2,ncol=M)
rownames(cv.performance) = c("rsq","pval")
colnames(cv.performance) = models

if ( opt$crossval <= 1 ) {
if ( opt$verbose >= 1 ) cat("### Skipping cross-validation\n")
} else {
if ( opt$verbose >= 1 ) cat("### Performing",opt$crossval,"fold cross-validation\n")
cv.all = pheno
N = nrow(cv.all)
cv.sample = sample(N)
cv.all = cv.all[ cv.sample , ]
folds = cut(seq(1,N),breaks=opt$crossval,labels=FALSE)

cv.calls = matrix(NA,nrow=N,ncol=M)

for ( i in 1:opt$crossval ) {
	if ( opt$verbose >= 1 ) cat("- Crossval fold",i,"\n")
	indx = which(folds==i,arr.ind=TRUE)
	cv.train = cv.all[-indx,]
	# store intercept
	intercept = mean( cv.train[,3] )
	cv.train[,3] = scale(cv.train[,3])
	
	# hide current fold
	cv.file = paste(opt$tmp,".cv",sep='')
	write.table( cv.train , quote=F , row.names=F , col.names=F , file=paste(cv.file,".keep",sep=''))	
	arg = paste( opt$PATH_plink ," --allow-no-sex --bfile ",opt$tmp," --keep ",cv.file,".keep --out ",cv.file," --make-bed",sep='')
	system(arg , ignore.stdout=SYS_PRINT,ignore.stderr=SYS_PRINT)

	for ( mod in 1:M ) {
		if ( models[mod] == "blup" ) {
			pred.wgt = weights.bslmm( cv.file , bv_type=2 , snp=genos$bim[,2] )
		}
		else if ( models[mod] == "bslmm" ) {
			pred.wgt = weights.bslmm( cv.file , bv_type=1 , snp=genos$bim[,2] )
		}		
		else if ( models[mod] == "lasso" ) {
			pred.wgt = weights.lasso( cv.file , hsq[1] , snp=genos$bim[,2] )
		}
		else if ( models[mod] == "enet" ) {
			pred.wgt = weights.enet( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3]) , alpha=0.5 )
		}		
		else if ( models[mod] == "top1" ) {
			pred.wgt = weights.marginal( genos$bed[ cv.sample[ -indx ],] , as.matrix(cv.train[,3,drop=F]) , beta=T )
			pred.wgt[ - which.max( pred.wgt^2 ) ] = 0
		}

		# predict from weights into sample
		pred.wgt[ is.na(pred.wgt) ] = 0
		cv.calls[ indx , mod ] = genos$bed[ cv.sample[ indx ] , ] %*% pred.wgt
	}
}

# compute rsq + P-value for each model
for ( mod in 1:M ) {
	if ( !is.na(sd(cv.calls[,mod])) && sd(cv.calls[,mod]) != 0 ) {
		reg = summary(lm( cv.all[,3] ~ cv.calls[,mod] ))
		cv.performance[ 1, mod ] = reg$adj.r.sq
		cv.performance[ 2, mod ] = reg$coef[2,4]
	} else {
		cv.performance[ 1, mod ] = NA
		cv.performance[ 2, mod ] = NA
	}
}
if ( opt$verbose >= 1 ) write.table(cv.performance,quote=F,sep='\t')
}

# --- FULL ANALYSES
if ( opt$verbose >= 1 ) cat("Computing full-sample weights\n")

# call models to get weights
wgt.matrix = matrix(0,nrow=nrow(genos$bim),ncol=M)
colnames(wgt.matrix) = models
rownames(wgt.matrix) = genos$bim[,2]
for ( mod in 1:M ) {
	if ( models[mod] == "blup" ) {
		wgt.matrix[,mod] = weights.bslmm( geno.file , bv_type=2 , snp=genos$bim[,2] , out=opt$tmp )
	}
	else if ( models[mod] == "bslmm" ) {
		wgt.matrix[,mod] = weights.bslmm( geno.file , bv_type=1 , snp=genos$bim[,2] , out=opt$tmp )
	}		
	else if ( models[mod] == "lasso" ) {
		wgt.matrix[,mod] = weights.lasso( geno.file , hsq[1] , snp=genos$bim[,2] , out=opt$tmp )
	}
	else if ( models[mod] == "enet" ) {
		wgt.matrix[,mod] = weights.enet( genos$bed , as.matrix(pheno[,3]) , alpha=0.5 )
	}	
	else if ( models[mod] == "top1" ) {
		wgt.matrix[,mod] = weights.marginal( genos$bed , as.matrix(pheno[,3]) , beta=F )
	}
}

# save weights, rsq, p-value for each model, and hsq to output
snps = genos$bim
save( wgt.matrix , snps , cv.performance , hsq, hsq.pv, N.tot , file = paste( opt$out , ".wgt.RDat" , sep='' ) )
# --- CLEAN-UP
if ( opt$verbose >= 1 ) cat("Cleaning up\n")
cleanup()
