# pipe
library(magrittr)

# args
args <-
  commandArgs(trailingOnly = T) %>% 
  setNames(c("input_dir", "sample_metadata_file", "sample_subset", "cell_subset")) %>%
  strsplit(",") 


seu_ls <-
  purrr::map2(
    args$input_dir, args$sample_metadata_file, 
    function(input_dir, sample_metadata_file) {
      
      # input dir
      input_files <- list.files(input_dir)
      
      # detect platform
      platform_files <- list(
        "Parse Biosciences" = c("DGE.mtx", "all_genes.csv", "cell_metadata.csv"),
        "10X Genomics" = c("matrix.mtx.gz", "features.tsv.gz", "genes.tsv.gz", "barcodes.tsv.gz")
      )
      platform <-
        platform_files %>%
        purrr::map(~ (.x %>% intersect(input_files) %>% length) == 3) %>%
        unlist() %>%
        {platform_files[.]} %>%
        names
      
      # read
      cat("Platform", platform, "detected. Reading...\n")
      if (platform == "10X Genomics") {
        
        # read in expression matrix
        dge_mat <- Seurat::Read10X(input_dir)
        
        # create seurat object
        seu <- Seurat::CreateSeuratObject(counts = dge_mat)
        
      } else if (platform == "Parse Biosciences") {
        
        # read in expression matrix
        dge_mat <- Seurat::ReadParseBio(input_dir)
        
        # read in cell metadata
        cell_metadata <-
          input_dir %>%
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
        
      # conditionally add sample metadata
      cat("Adding sample metadata...\n")
      if (sample_metadata_file != "false") {
        
        # read in sample metadata and add to misc slot
        seu@misc$sample_metadata <-
          readr::read_tsv(sample_metadata_file, show_col_types = F)
        
        # add to Seurat object meta.data
        seu@meta.data <-
          dplyr::left_join(
            seu@meta.data,
            seu@misc$sample_metadata,
            by = "sample")
        rownames(seu@meta.data) <- colnames(seu)
        
      }
        
      # remove genes with no counts
      counts <- as.matrix(seu@assays$RNA@counts)
      seu <- subset(seu, features = rownames(counts[rowSums(counts) > 0, ]))
      
      # remove NA samples
      seu <- subset(seu, subset = sample %in% seu$sample[!is.na(seu$sample)])
      
      # subset to cell subset
      if (all(args$cell_subset != "false")) {
        seu <- seu[, intersect(args$cell_subset, colnames(seu))]
      }
      
      # subset to sample subset
      if (all(args$sample_subset != "false")) {
        
        if (length(intersect(args$sample_subset, seu$sample)) == 0) {
          
          seu <- NULL
          
        } else {
          
          seu <- subset(x = seu, subset = sample %in% args$sample_subset)
          seu@misc$sample_metadata <-
            seu@misc$sample_metadata %>%
            dplyr::filter(sample %in% args$sample_subset)
          
        }
        
      }
        
      # return
      seu
      
    }
  ) 

# drop empty elements
seu_ls <- Filter(Negate(is.null), seu_ls)

# merge if multiple runs, combine misc slot
if (length(seu_ls) > 1) {
  seu <- Reduce(merge, seu_ls)
  seu@misc <- do.call(Map, c(rbind, purrr::map(seu_ls, ~ .x@misc))) %>% purrr::map(dplyr::distinct)
} else {
  seu <- seu_ls[[1]]
}

# save
saveRDS(seu, "seu.rds")
