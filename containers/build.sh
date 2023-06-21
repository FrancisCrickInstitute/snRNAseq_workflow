#!/bin/bash
#SBATCH --job-name=snRNAseq_workflow
#SBATCH --time=1-00:00:0
#SBATCH --mem=10GB
#SBATCH -o /nemo/project/proj-tracerX/working/VHL_GERMLINE/tidda/snRNAseq_workflow/sbatch_log/sbatch.out
#SBATCH -e /nemo/project/proj-tracerX/working/VHL_GERMLINE/tidda/snRNAseq_workflow/sbatch_log/sbatch.err

cd /nemo/project/proj-tracerX/working/VHL_GERMLINE/tidda/snRNAseq_workflow/containers/
ml purge
. ~/.bashrc
ml Singularity/3.6.4
rm conda.sif
singularity build --fakeroot conda.sif singularity.def


# ml purge
# module load anaconda3/personal
# mamba create -n snRNAseq_workflow r-essentials r-base=4.2.0
# mamba activate snRNAseq_workflow
# mamba install \
#   -c bioconda \
#   -c bu_cnio \
#   -c bioconda \
#   -c conda-forge \
#   nextflow
#   r-rjson
#   r-seuratwrappers \
#   r-seurat \
#   bioconductor-summarizedexperiment \
#   bioconductor-singler \
#   r-dplyr \
#   r-tidyr \
#   r-purrr \
#   r-knitr \
#   r-readr \
#   r-rlang \
#   r-rmarkdown \
#   r-magrittr \
#   r-patchwork \
#   r-pheatmap \
#   r-ggrepel \
#   r-intrinsicdimension \
#   r-lemon \
#   r-clustree \
#   r-ggforce \
#   r-ggplot2 \
#   r-gridextra \
#   r-xfun \
#   bioconductor-celldex \
#   r-matrix \
#   bioconductor-scdblfinder \
#   r-stringr \
#   r-reshape2 \
#   r-rcolorbrewer \
#   bioconductor-dittoseq \
#   r-monocle3 \
#   bioconductor-infercnv
