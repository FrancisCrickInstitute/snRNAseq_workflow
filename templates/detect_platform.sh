#!/bin/bash

set -euo pipefail

# check directory structure
if [ -f "test/matrix.mtx" -a -f "test/genes.tsv" -a -f "test/barcodes.tsv" ] ; then
    echo "10X Genomics platform detected"
    platform="10X Genomics"
elif [ -f "test/DGE.mtx" -a -f "test/all_genes.csv" -a -f "test/cell_metadata.csv" ] ; then
    echo "Parse Biosciences platform detected"
    platform="Parse Biosciences"
fi