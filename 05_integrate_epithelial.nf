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

// import modules
include { save_params } from './modules/save_params'

// integrating
// integrate malignant samples
process malignant_integrating {
  tag "${id}"
  label 'process_high_memory'
  publishDir "${output_dir}",
    mode: 'copy',
    pattern: "*.rds"

  input:
    path rmd_file
    path rds_files, stageAs: "seu??.rds"
    path params_file
    val output_dir

  output:
    path 'seu.rds', emit: ch_integrated
    path 'integrating.html'

  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        params_file = "${params_file}",
        rds_file = "${rds_files.join(',')}",
        cache_dir = "${output_dir}/integrating_cache/"),
      output_file = "integrating.html",
      output_dir = getwd()
    )
    """
}

// cluster malignant samples
process malignant_clustering {
  tag "${id}"
  label 'process_high_memory'
  cpus = { check_max( 12 * task.attempt, 'cpus' ) }
  publishDir "${output_dir}",
    mode: 'copy',
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"

  input:
    path rmd_file
    path rds_file
    path params_file
    val output_dir

  output:
    path 'clustering.html'
    path 'clustering_files/figure-html/*.png'

  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        params_file = "${params_file}",
        rds_file = "${rds_file}",
        cache_dir = "${output_dir}/clustering_cache/"),
      output_file = "clustering.html",
      output_dir = getwd()
    )
    """
}

workflow {

  // save params
  save_params()

  // integrate all preloaded samples
  Channel
    .fromPath("${params.output.dir}/by_patient_wo_organoids/*/integrating/infercnv/seu_infercnv_malignant.rds")
    .set { ch_run }

  // integrate all malignant cells
  malignant_integrating(
    "${baseDir}/templates/integrating.rmd",
    ch_run,
    save_params.out.ch_params,
    "${params.output.dir}/malignant/integrating/"
  )

  // cluster all malignant cells
  malignant_clustering(
    "${baseDir}/templates/seurat_clustering.rmd",
    malignant_integrating.out.ch_integrated
    save_params.out.ch_params,
    "${params.output.dir}/malignant/integrating/clustering/"    
  )
}
