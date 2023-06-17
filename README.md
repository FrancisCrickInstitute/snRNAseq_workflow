# snRNAseq nextflow workflow

The essential components of the workflow repository are:
- `main.nf`: Contains the primary workflow code which pulls in all additional code from the repository
- `modules/`: Contains all of the sub-workflows which are used to organize large chunks of analysis
- `templates/`: Contains all of the code which is executed in each individual step of the workflow