snRNAseq_workflow image built as a modification of 'almurphy-scfdev-dev.img'.

Build a sandbox of the image, enter the sandbox and open R:

```
singularity build --sandbox snRNAseq_workflow/ docker://almurphy/scfdev:dev
$ singularity shell -f --writable snRNAseq_workflow/
INFO:    User not listed in /etc/subuid, trying root-mapped namespace
INFO:    Using fakeroot command combined with root-mapped namespace
Apptainer> R
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

Install a previous version of EWCE from GitHub:

```
devtools::install_github("neurogenomics/EWCE")
```

The install PyTorch and scGPT using pip:

```
pip install torchvision
pip install cuda-python==11.7.0
pip install scgpt
```

Then exit the container (`exit`) and convert it to an image:

```
singularity build snRNAseq_workflow.img snRNAseq_workflow/
```


