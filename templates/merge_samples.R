# pipe
library(magrittr)

# args
args <-
  commandArgs(trailingOnly = T) %>% 
  setNames(c("rds_files")) %>%
  strsplit(",") 

# read files
seu_ls <-
  args$rds_files %>%
  lapply(readRDS)

# merge
seu <- Reduce(merge, seu_ls)
seu@misc$sample_metadata <-
  seu_ls %>% 
  purrr::map(~ .x@misc$sample_metadata) %>% 
  dplyr::bind_rows() %>% 
  dplyr::distinct()

# save
saveRDS(seu, "seu.rds")