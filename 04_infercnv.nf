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
include { infercnv    } from './modules/infercnv'

// integrate malignant samples
process malignant_integrating {
  tag "${id}"
  label 'process_medium'
  publishDir "${params.output.dir}/${dir}/${id}/integrating/infercnv/malignant_integrating/",
    mode: 'copy',
    pattern: "*.rds"

  input:
    path rmd_file
    tuple val(id), val(dir), val(subdir), path(rds_files, stageAs: "seu??.rds")
    path params_file

  output:
    tuple val(id), val(dir), val(subdir), path('seu.rds'), emit: ch_integrated
    path 'integrating.html'

  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        params_file = "${params_file}",
        rds_file = "${rds_files.join(',')}",
        cache_dir = "${params.output.dir}/${dir}/${id}/integrating/infercnv/malignant_integrating/malignant_integrating_cache/"),
      output_file = "integrating.html",
      output_dir = getwd()
    )
    """
}

// cluster malignant samples
process malignant_clustering {
  tag "${id}"
  label 'process_medium'
  cpus = { check_max( 12 * task.attempt, 'cpus' ) }
  publishDir "${params.output.dir}/${dir}/${id}/integrating/infercnv/malignant_clustering/",
    mode: 'copy',
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"

  input:
    path rmd_file
    tuple val(id), val(dir), val(subdir), path(rds_files, stageAs: "seu??.rds")
    path params_file

  output:
    tuple val(id), val(dir), val(subdir), path('seu_annotated.rds'), emit: ch_clustered
    path 'malignant_clustering.html'
    path 'malignant_clustering_files/figure-html/*.png'

  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        params_file = "${params_file}",
        rds_file = "${rds_file}",
        cache_dir = "${params.output.dir}/${dir}/${id}/integrating/infercnv/malignant_clustering/malignant_clustering_cache/"),
      output_file = "malignant_clustering.html",
      output_dir = getwd()
    )
    """
}

workflow {

  // save params
  save_params()

  // integrate all preloaded samples
  Channel
    .fromPath("${params.output.dir}/by_patient_wo_organoids/*/integrating/seurat_clustering/seu_annotated.rds")
    .set { ch_run }

  // get id, dir, subdir and filepath
  ch_run
    .map { file ->
            tuple(file.toString().tokenize('/')[-4],
                  file.toString().tokenize('/')[-5],
                  file.toString().tokenize('/')[-3],
                  file) }
    .set { ch_run }

  // perform infercnv
  infercnv(
    "${baseDir}/templates/infercnv.rmd",
    ch_run,
    save_params.out.ch_params
  )
  
  // perform malignant integration
  malignant_integrating(
    "${baseDir}/templates/integrating.rmd",
    infercnv.out.ch_malignant,
    save_params.out.ch_params
  )
  
  // perform malignant clustering
  malignant_clustering(
    "${baseDir}/templates/integrating.rmd",
    malignant_integrating.out.ch_integrated,
    save_params.out.ch_params
  )
  
}
