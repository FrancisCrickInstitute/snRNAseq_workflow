# snRNAseq nextflow workflow

The essential components of the workflow repository are:
- `main.nf`: Contains the primary workflow code which pulls in all additional code from the repository
- `modules/`: Contains all of the sub-workflows which are used to organize large chunks of analysis
- `templates/`: Contains all of the code which is executed in each individual step of the workflow

```
ml purge
ml Nextflow/22.10.3
ml Singularity/3.6.4
ml CAMP_proxy

cd /nemo/project/proj-tracerX/working/VHL_GERMLINE/tidda/snRNAseq_workflow/

nextflow run . -c config/VHL_ParseBio.config -resume
```