#!/bin/bash

cd /rds/general/user/art4017/home/snRNAseq_workflow/

. ~/.bashrc
conda activate nfcore
module load gcc/8.2.0
NXF_OPTS='-Xms1g -Xmx4g'

nextflow run . \
  -c config/oesophageal_10X.config \
  -c config/imperial.config \
  -profile imperial \
  -with-singularity /rds/general/user/art4017/home/nf-core/scflow/work/singularity/snRNAseq_workflow \
  --input.manifest_file input/oesophageal_10X/test/sample_manifest.tsv \
  --output.dir output/test/ \
  -w work_test/
