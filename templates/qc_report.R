# define output directory
dir.create(out$base, showWarnings = F, recursive = T)

# set sample groupings to check at the clustering stage
groupings <-
  c("sample", "percent_mito", "percent_ribo", "percent_globin",
    "date_prep", "nih_pid", "rin", "lesion_type", "tumour_size", "fuhrman_grade") %>%
  # if genome is human, do cell cycle scoring (doesn't work with other genomes)
  { if (grepl("hg38", genome)) c(., "Phase") else . } %>%
  # if checking for doublets, add to the groupings
  { if (remove_doublets) c(., "doublet") else . } %>%
  # if multiple experiments, add to the groupings
  { if (length(experiment) > 1) c(., "experiment") else . }

# statistics of interest in the Parse analysis
statistics <-
  c("median_tscp_per_cell",
    "median_genes_per_cell",
    "fraction_tscp_in_cells",
    "cell_tscp_cutoff",
    "valid_barcode_fraction"
  )

# save args to output directory
dput(args, paste0(out$base, "/args_for_generate_qc_report.R"))

# pass args to rmd
rmarkdown::render(system.file("rmd", "generate_qc_report.rmd", package = "vhl"),
                  knit_root_dir = rprojroot::find_rstudio_root_file(),
                  output_dir = out$base,
                  output_file = "qc_report",
                  params = list(args, out = out,
                                groupings = groupings, statistics = statistics))
