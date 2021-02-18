#!/bin/sh
# MAKE SURE FUSION.compute_weights.R IS IN YOUR PATH
# FILL IN THESE PATHS
GCTA="<PATH TO GCTA>"
PLINK="<PATH TO PLINK>"
GEMMA="<PATH TO GEMMA>"
# ALTERNATIVELY: ENSURE THAT plink, gcta, gemma CAN BE CALLED FROM PATH AND REMOVE --PATH_* FLAGS BELOW
# PATH TO DIRECTORY CONTAINING LDREF DATA (FROM FUSION WEBSITE or https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2)
LDREF="<PATH TO LDREF>"
# THIS IS USED TO RESTRICT INPUT SNPS TO REFERENCE IDS ONLY

# PATH TO GEUVADIS GENE EXPRESSION MATRIX:
PRE_GEXP="GD462.GeneQuantRPKM.50FN.samplename.resk10.txt"
# GEUVADIS DATA WAS DOWNLOADED FROM https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/files/analysis_results/
# GENERATE ID FILE
cat $PRE_GEXP | head -n1 | tr '\t' '\n' | tail -n+5 | awk '{ print $1,$1 }' > ${PRE_GEXP}.ID

# PATH TO PREFIX FOR GEUVADIS GENOTYPES SPLIT BY CHROMOSOME
# SUBSAMPLE THESE TO THE LDREF SNPS FOR EFFICIENCY
PRE_GENO="1000G"

# PATH TO OUTPUT DIRECTORY (population-specific subdirs will be made)
OUT_DIR="./WEIGHTS"

# ROWS IN THE MATRIX TO ANALYZE (FOR BATCHED RUNS)
BATCH_START=1
BATCH_END=10

# --- BEGIN SCRIPT:

NR="${BATCH_START}_${BATCH_END}"
mkdir --parents tmp/$NR
mkdir --parents hsq/$NR
mkdir --parents out/$NR
# THIS IS DIRECTORY WHERE THE OUTPUT WILL GO:
mkdir $OUT_DIR

# Loop through each gene expression phenotype in the batch
cat $PRE_GEXP | awk -vs=$BATCH_START -ve=$BATCH_END 'NR > s && NR <= e' |  while read PARAM; do

# Get the gene positions +/- 500kb
CHR=`echo $PARAM | awk '{ print $3 }'`
P0=`echo $PARAM | awk '{ print $4 - 0.5e6 }'`
P1=`echo $PARAM | awk '{ print $4 + 0.5e6 }'`
GNAME=`echo $PARAM | awk '{ print $1 }'`

OUT="tmp/$NR/$PRE_GEXP.$GNAME"

echo $GNAME $CHR $P0 $P1

# Pull out the current gene expression phenotype
echo $PARAM | tr ' ' '\n' | tail -n+5 | paste $PRE_GEXP.ID - > $OUT.pheno

# Get the locus genotypes for all samples and set current gene expression as the phenotype
$PLINK --bfile $PRE_GENO.$CHR --pheno $OUT.pheno --make-bed --out $OUT --keep $OUT.pheno --chr $CHR --from-bp $P0 --to-bp $P1 --extract $LDREF/1000G.EUR.$CHR.bim

# Process all samples together (for reference purposes only since this is mult-ethnic data)
mkdir $OUT_DIR/ALL
FINAL_OUT="$OUT_DIR/ALL/ALL.$GNAME"

Rscript FUSION.compute_weights.R --bfile $OUT --tmp $OUT.$pop.tmp --out $FINAL_OUT --verbose 0 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA --models blup,lasso,top1,enet
# ALTERNATIVELY ADD COVARIATES HERE USING THE --covar FLAG
# MINIMAL COMMAND IS: `Rscript FUSION.compute_weights.R --bfile $OUT --tmp $OUT.$pop.tmp --out $FINAL_OUT`

# Append heritability output to hsq file
cat $FINAL_OUT.hsq >> hsq/$NR.hsq

# Clean-up just in case
rm -f $FINAL_OUT.hsq $OUT.tmp.*

# Process each population
# THIS REQUIRES HAVING A $pop.keep FILE LISTING IDS FOR EACH POPULATION
for pop in EUR YRI TSI CEU FIN GBR; do
	mkdir --parents WEIGHTS/$pop
	FINAL_OUT="WEIGHTS/$pop/$pop.$GNAME"
	# EXTRACT THE POPULATION AND 1% SNPS
	$PLINK --bfile $OUT --keep $pop.keep --make-bed --out $OUT.$pop --maf 0.01

	# MAKE SURE FUSION.compute_weights.R IS IN YOUR PATH
	Rscript FUSION.compute_weights.R --bfile $OUT.$pop --tmp $OUT.$pop.tmp --out $FINAL_OUT --verbose 0 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA --models blup,lasso,top1,enet
	
	# Append heritability output to hsq file
	cat $FINAL_OUT.hsq >> hsq/$NR.hsq

	# Clean-up just in case
	rm -f $FINAL_OUT.hsq $OUT.tmp.*
done

# Remove all intermediate files
rm $OUT.*

# GO TO THE NEXT GENE
done
