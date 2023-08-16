#!/usr/bin/env nextflow

// default params set in nextflow.config

// import json function
import java.io.BufferedReader

// using DSL-2
nextflow.enable.dsl=2

log.info """\

snRNAseq workflow
=================
manifest file:        ${params.input.manifest_file}
sample metadata file: ${params.input.sample_metadata_file}
output dir:           ${params.output.dir}
"""

// merge infercnv samples
process merge_infercnv {
  tag "${patient}"
  label 'process_medium'
  publishDir "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/",
    mode: 'copy',
    pattern: "seu.rds"

  input:
    tuple val(patient), path(rds_file)

  output:
    tuple val(patient), path('seu.rds'), optional: true

  script:
    """
    #!/usr/bin/env Rscript
    library(magrittr)
    
    # read files
    seu <- readRDS("seu_annotated.rds")
    
    # get cnv files
    infercnv_cache_dir <- 
      "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/infercnv_cache/"
    cnv_files <- 
      list.files(
        infercnv_cache_dir,
        pattern = "HMM_CNV_predictions.HMMi3.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat",
        recursive = T
      )
      
    # continue if cnv files
    if (length(cnv_files) > 0) {
      aneuploid_cells <- 
        cnv_files %>%
        purrr::map(function(cnv_file_i) {
          sample_i <- gsub("/.*", "", cnv_file_i)
          readr::read_tsv(paste0(infercnv_cache_dir, cnv_file_i)) %>%
            dplyr::mutate(sample = sample_i)
        }) %>%
        dplyr::bind_rows() %>%
        dplyr::group_by(sample, cell_group_name) %>%
        dplyr::summarise(aneuploid = any(c(1, 3) %in% state)) %>%
        dplyr::filter(aneuploid)
      malignant_cells <-
        aneuploid_cells %>%
        dplyr::filter(cell_group_name %in% colnames(seu[, seu\$label.reduced == "Epithelial_cells"]))
      
      # exclude epithelial samples in samples for which infercnv run was unsuccessful 
      # (malignant annotations not present)
      seu <- seu[, (seu\$sample %in% gsub("/.*", "", cnv_files) |
                    seu\$label.reduced != "Epithelial_cells")]
                    
      # add aneuploid and malignant labels 
      seu\$aneuploid <- 0
      seu@meta.data[aneuploid_cells\$cell_group_name, "aneuploid"] <- 1
      seu\$label.infercnv <- seu\$label.reduced
      seu@meta.data[malignant_cells\$cell_group_name, "label.infercnv"] <- "Malignant_cells"
      
      # save
      saveRDS(seu, "seu.rds")
    }
    """
}

workflow {

  // get files
  Channel
    .fromPath("${params.output.dir}/by_patient_wo_organoids/*/integrating/seurat_clustering/seu_annotated.rds")
    .map { file -> tuple(file.toString().tokenize('/')[-4], file) }
    .groupTuple()
    .set { ch_run }
    
  // merge infercnv
  merge_infercnv(
    ch_run
  )
  
}

