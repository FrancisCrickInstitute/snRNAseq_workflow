snRNAseq_workflow image built as a modification of 'almurphy-scfdev-dev.img'.

Build a sandbox of the image, enter the sandbox and open R:

```
singularity build --sandbox snRNAseq_workflow/ docker://almurphy/scfdev:dev
singularity shell --writable snRNAseq_workflow/
Apptainer *:~> R
```

Then in R inside the container, run install.packages() or BiocManager::install() to download the following packages in this order:

* dplyr
* Seurat
* SeuratWrapper
* inferCNV
* SingleR
* scDblFinder
* celldex
* Matrix
* stringr
* reshape2
* RColorBrewer
* ggplot2
* pheatmap
* purrr
* rmarkdown
* xfun
* Monocle3
* rjson
* clustree
* dittoSeq
* ggforce
* ggrepel
* magrittr
* readr
* SummarizedExperiment
* janitor
* gridExtra
* intrinsicDimension
* knitr
* lemon
* patchwork
* rlang

Then exit the container (`exit`) and convert it to an image:

```
singularity build snRNAseq_workflow.img snRNAseq_workflow/
```


