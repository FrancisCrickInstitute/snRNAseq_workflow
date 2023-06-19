#!/bin/bash
#SBATCH --job-name=snRNAseq_workflow
#SBATCH --time=1-00:00:0
#SBATCH --mem=10GB
#SBATCH -o /nemo/project/proj-tracerX/working/VHL_GERMLINE/tidda/snRNAseq_workflow/sbatch_log/sbatch.out
#SBATCH -e /nemo/project/proj-tracerX/working/VHL_GERMLINE/tidda/snRNAseq_workflow/sbatch_log/sbatch.err

cd /nemo/project/proj-tracerX/working/VHL_GERMLINE/tidda/snRNAseq_workflow/

# env
ml purge
. ~/.bashrc
conda activate vhl2

# run
nextflow run . \
  -c config/VHL_ParseBio.config \
  -c config/crick.config \
  -resume
