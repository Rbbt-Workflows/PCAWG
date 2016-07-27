#!/bin/bash

BASE_DIR='/home/mvazquezg/git/workflows/PCAWG/'
PROJECT_DIR="$BASE_DIR/project/preliminary_final_release"
DATA_DIR="$BASE_DIR/data/preliminary_final_release"

# PREPARE GENOTYPES

cd $PROJECT_DIR
tar xvfz $DATA_DIR/preliminary_final_release.snvs.tgz

rm preliminary_final_release/snv_mnv/*.tbi

mkdir genotypes
for f in  preliminary_final_release/snv_mnv/*.gz; do
  clean_name=$(basename $f .annotated.snv_mnv.vcf.gz).vcf
  zcat $f |grep -v LOWSUPPORT > "genotypes/$clean_name"
  gzip "genotypes/$clean_name"
done

rm -Rf preliminary_final_release

# GET SAMPLE INFO


