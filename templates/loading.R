# pipe
library(magrittr)

# args
args <-
  commandArgs(trailingOnly = T) %>% 
  setNames(c("id", "id_col", "dir", "sample_metadata_file")) %>%
  strsplit(",") %>%
  as.list()

# platform specs
platform_files <- list(
  "Parse Biosciences" = c("DGE.mtx", "all_genes.csv", "cell_metadata.csv"),
  "10X Genomics" = c("matrix.mtx.gz", "features.tsv.gz", "genes.tsv.gz", "barcodes.tsv.gz")
)

seu_ls <-
  args$dir %>%
  purrr::map(function(dir) {
    
    # input files
    input_files <- list.files(dir)
    
    # detect platform
    platform <-
      platform_files %>%
      purrr::map(~ (.x %>% intersect(input_files) %>% length) == 3) %>%
      unlist() %>%
      {platform_files[.]} %>%
      names()
    
    if (platform == "10X Genomics" & args$id_col != "sample") {
      stop("10X Genomics implementation only supports sample IDs!")
    }
    
    # read
    cat("Platform", platform, "detected. Reading...\n")
    if (platform == "10X Genomics") {
      
      # read in expression matrix
      dge_mat <- Seurat::Read10X(dir)
      
      # create seurat object
      seu <- Seurat::CreateSeuratObject(counts = dge_mat)
      
      # add id column
      seu@meta.data[, args$id_col] <- args$id
      
    } else if (platform == "Parse Biosciences") {
      
      # read in expression matrix
      dge_mat <- Seurat::ReadParseBio(dir)
      
      # read in cell metadata
      cell_metadata <-
        dir %>%
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
    
    # add input path to metadata
    seu@meta.data$dir <- dir
    
    # return
    seu
    
  })

# merge
seu <- Reduce(merge, seu_ls)

# add sample metadata
cat("Adding sample metadata...\n")

# read in sample metadata and add to Seurat misc slot, subset to sample, remove columns that are all NA
seu@misc$sample_metadata <-
  readr::read_tsv(args$sample_metadata_file, show_col_types = F) 

# add to Seurat meta.data slot
seu@meta.data <-
  dplyr::left_join(
    seu@meta.data,
    seu@misc$sample_metadata,
    by = "sample")
rownames(seu@meta.data) <- colnames(seu)

# optionally subset Seurat object to id
if (args$id != "") {
  cat("Subsetting", args$id_col, "to", args$id, "...\n")
  seu <- seu[, seu@meta.data[, args$id_col] == args$id] 
  seu@misc$sample_metadata <- 
    seu@misc$sample_metadata %>%
    dplyr::filter(get(args$id_col) == args$id)
}
cat("Done subsetting!\n")

# remove columns that are all NA
seu@misc$sample_metadata <-
  seu@misc$sample_metadata[, colSums(is.na(seu@misc$sample_metadata)) != nrow(seu@misc$sample_metadata)]

# remove genes with no counts
counts <- as.matrix(seu@assays$RNA@counts)
seu <- subset(seu, features = rownames(counts[rowSums(counts) > 0, ]))

# save
saveRDS(seu, "seu.rds")
