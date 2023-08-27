# snRNAseq Nextflow workflow

The essential components of the workflow repository are:
- `[##]_*.nf`: Scripts containing the primary workflow code for each step
- `modules/`: Contains all of the sub-workflows which are used to organize large chunks of analysis
- `templates/`: Contains all of the code which is executed in each individual step of the workflow, as R and bash scripts
- `config/`: Contains configuration files for different datasets and institutional cluster configurations
- `input/`: Contains patient and sample manifests for different input datasets.

PBSPro scripts for running each of the steps from this workflow can be found in the [snRNAseq_analysis](github.com/alextidd/snRNAseq_analysis) repository.

Before running Nextflow scripts on the Imperial clusters, the following must be run:

```
. ~/.bashrc
conda activate nfcore
module load gcc/8.2.0
NXF_OPTS='-Xms1g -Xmx4g'
```
