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
seu@misc <- 
  do.call(Map, c(dplyr::bind_rows, 
    seu_ls %>%
      purrr::map(function(i) {
        i@misc %>%
          purrr::map(function(slot) {
            if (is.data.frame(slot)) {
              slot
            }
          })
      }) 
    )
  )

# save
saveRDS(seu, "seu.rds")