#!/bin/bash

module load plink2
GCTA="./gcta_nr_robust"
PLINK="plink --allow-no-sex"
GEMMA="./gemma"
FUSION="./FUSION/"

PRE="Adipose_Subcutaneous"
PRE="Prostate"
PRE="Whole_Blood"
NR=10

PRE=$1
NR=${SLURM_ARRAY_TASK_ID}

# ---

mkdir tmp/$PRE.$NR
rm -f HSQ/$PRE.$NR.hsq HSQ/TSPEC.$PRE.$NR.hsq

COVAR="$PRE.v7.covariates.txt.covar"
zcat $PRE.v7.normalized_expression.bed.gz | tail -n+2 | awk -v i=$NR 'NR > (i-1)*100 && NR <= i*100' |  while read PARAM; do

GNAME=`echo $PARAM | awk '{ print $4 }'`
CHR=`echo $PARAM | awk '{ print $1 }'`
P0=`echo $PARAM | awk '{ p=$2 - 500e3; if(p<0) p=0; print p; }'`
P1=`echo $PARAM | awk '{ print $3 + 500e3 }'`

OUT="tmp/$PRE.$NR/$PRE.$GNAME"

echo $GNAME $CHR $P0 $P1
echo $PARAM | tr ' ' '\n' | paste $PRE.v7.normalized_expression.bed.gz.HEADER $PRE.v7.normalized_expression.bed.gz.HEADER - | tail -n+5 > $OUT.pheno

rm -f $OUT.bed
$PLINK --silent --bfile GTEx.v7.HM3 --chr $CHR --from-bp $P0 --to-bp $P1 --make-bed --out $OUT --pheno $OUT.pheno --keep $OUT.pheno --maf 0.0001
if [ ! -f $OUT.bed ]; then
continue
fi

FINAL_OUT="WEIGHTS/$PRE/$PRE.$GNAME"
Rscript $FUSION/FUSION.compute_weights.R \
--bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --verbose 0 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA --models blup,lasso,top1,enet --covar $PRE.v7.covariates.txt.covar --hsq_p 1.0
cat $FINAL_OUT.hsq | awk -vw="$PRE $w" '{ print w,$0 }' >> HSQ/$PRE.$NR.hsq
rm -f $FINAL_OUT.hsq

# Get best eQTL for this gene
META=`cat GTEx.v7.meta.id | awk -vg=$GNAME '$3 == g { print $2 }'`
if [ "$META" == "" ]; then
continue
fi
cat GTEx.v7.meta.covar | awk -vm=$META '$1 == m' | tr '\t' '\n' | paste GTEx.v7.meta.covar.HEADER - | tail -n+2 | sort -k1,1 \
| join -1 1 -2 1 - <(cut -f2- $COVAR | sort -k1,1) | awk '$2 != "NA" { print $1,$0 }' | tr -s ' ' '\t' > $OUT.meta.covar

FINAL_OUT="WEIGHTS/$PRE.TSPEC/$PRE.$GNAME"

Rscript $FUSION/FUSION.compute_weights.R \
--bfile $OUT --tmp $OUT.tmp --out $FINAL_OUT --verbose 0 --save_hsq --PATH_gcta $GCTA --PATH_gemma $GEMMA --models blup,lasso,top1,enet --covar $OUT.meta.covar --resid
cat $FINAL_OUT.hsq | awk -vw="$PRE $w" '{ print w,$0 }' >> HSQ/TSPEC.$PRE.$NR.hsq
rm -f $FINAL_OUT.hsq

rm $OUT.*
done

rm -fr tmp/$PRE.$NR
touch HSQ/$PRE.$NR.done
