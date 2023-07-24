# pipe
library(magrittr)

# # args
# args <-
#   commandArgs(trailingOnly = T) %>% 
#   setNames(c("rds_files")) %>%
#   strsplit(",") 
# 
# # read files
# seu_ls <-
#   args$rds_files %>%
#   lapply(readRDS)
# 
# if (length(seu_ls) > 1) {
#   # normalize and identify variable features for each dataset independently
#   seu_ls <- lapply(X = seu_ls, FUN = function(x) {
#     x <- Seurat::NormalizeData(x)
#     x <- Seurat::FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#     x
#   })
#   
#   # select features that are repeatedly variable across datasets for integration
#   features <- Seurat::SelectIntegrationFeatures(object.list = seu_ls)
#   
#   # identify integration anchors
#   seu_anchors <- 
#     Seurat::FindIntegrationAnchors(
#       object.list = seu_ls, 
#       anchor.features = features
#     )
#   
#   if (nrow(seu_anchors@anchors) < 100) {
#     # if n anchor cells < 100, merge
#     seu <- Reduce(merge, seu_ls)
#   } else {
#     # create 'integrated' data assay
#     seu <- Seurat::IntegrateData(anchorset = seu_anchors)
#   }
#   # merge metadata
#   seu@misc$sample_metadata <-
#     seu_ls %>% 
#     purrr::map(~ .x@misc$sample_metadata) %>% 
#     dplyr::bind_rows() %>% 
#     dplyr::distinct()
# } else {
#   seu <- Seurat::NormalizeData(seu_ls[[1]])
#   seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
# }

# # save
# saveRDS(seu, "seu.rds")

seu <- readRDS("seu.rds")
seu@misc$sample_metadata <-
  seu@meta.data %>%
  dplyr::select(4:23) %>% 
  tibble::as_tibble() %>% 
  dplyr::distinct() 
saveRDS(seu, "seu.rds")

