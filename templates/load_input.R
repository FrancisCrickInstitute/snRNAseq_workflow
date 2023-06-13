# pipe
library(magrittr)

# args
args <-
  commandArgs(trailingOnly = T) %>% 
  setNames(c("sample", "dir", "sample_metadata_file")) %>%
  strsplit(",") 

# input files
input_files <- list.files(args$dir)

# platform specs
platform_files <- list(
  "Parse Biosciences" = c("DGE.mtx", "all_genes.csv", "cell_metadata.csv"),
  "10X Genomics" = c("matrix.mtx.gz", "features.tsv.gz", "genes.tsv.gz", "barcodes.tsv.gz")
)

# detect platform
platform <-
  platform_files %>%
  purrr::map(~ (.x %>% intersect(input_files) %>% length) == 3) %>%
  unlist() %>%
  {platform_files[.]} %>%
  names()

# read
cat("Platform", platform, "detected. Reading...\n")
if (platform == "10X Genomics") {
  
  # read in expression matrix
  dge_mat <- Seurat::Read10X(args$dir)
  
  # create seurat object
  seu <- Seurat::CreateSeuratObject(counts = dge_mat)
  
} else if (platform == "Parse Biosciences") {
  
  # read in expression matrix
  dge_mat <- Seurat::ReadParseBio(args$dir)
  
  # read in cell metadata
  cell_metadata <-
    args$dir %>%
    paste0("/cell_metadata.csv") %>%
    readr::read_csv(show_col_types = F) %>%
    tibble::column_to_rownames("bc_wells") %>%
    as.data.frame()
  
  # create seurat object
  seu <-
    dge_mat %>%
    Seurat::CreateSeuratObject(
      names.field = 0,
      meta.data = cell_metadata, 
      min.cells = 1
    )
  
}

# add sample metadata
cat("Adding sample metadata...\n")
  
# read in sample metadata and add to misc slot
seu@misc$sample_metadata <-
  readr::read_tsv(args$sample_metadata_file, show_col_types = F)

# add to Seurat object meta.data
seu@meta.data <-
  dplyr::left_join(
    seu@meta.data,
    seu@misc$sample_metadata,
    by = "sample")
rownames(seu@meta.data) <- colnames(seu)

# remove genes with no counts
counts <- as.matrix(seu@assays$RNA@counts)
seu <- subset(seu, features = rownames(counts[rowSums(counts) > 0, ]))

# subset to sample
seu <- subset(x = seu, subset = sample == args$sample)
seu@misc$sample_metadata %>%
  dplyr::filter(sample == args$sample)

# save
saveRDS(seu, "seu.rds")
