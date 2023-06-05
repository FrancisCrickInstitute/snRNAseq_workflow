#!/usr/bin/env Rscript
rmarkdown::render("templates/generate_qc_report.rmd",
                  knit_root_dir = rprojroot::find_rstudio_root_file(),
                  output_dir = "output/",
                  output_file = "qc_report",
                  params = list(params, out = out,
                                groupings = groupings, statistics = statistics))