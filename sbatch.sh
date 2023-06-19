#!/bin/bash
#SBATCH --job-name=snRNAseq_workflow
#SBATCH --time=1-00:00:0
#SBATCH --mem=10GB

cd /nemo/project/proj-tracerX/working/VHL_GERMLINE/tidda/snRNAseq_workflow/

# env
ml purge
. ~/.bashrc
conda activate vhl2

# run
run . \
  -c config/VHL_ParseBio.config \
  -c config/crick.config \
  -resume