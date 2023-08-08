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

// integrate malignant samples
process integrating {
  tag "${patient}"
  label 'process_medium'
  publishDir "${output_dir}/integrating/",
    mode: 'copy',
    pattern: "*.rds"

  input:
    path rmd_file
    tuple val(patient), path(rds_files, stageAs: "seu??.rds"), val(output_dir)
    path params_file

  output:
    tuple val(patient), path('seu.rds'), val(output_dir), emit: ch_integrated
    path 'integrating.html'

  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        params_file = "${params_file}",
        rds_file = "${rds_files.join(',')}",
        cache_dir = "${output_dir}/integrating/integrating_cache/"),
      output_file = "integrating.html",
      output_dir = getwd()
    )
    """
}

// cluster malignant samples
process clustering {
  tag "${patient}"
  label 'process_medium'
  container '/rds/general/user/art4017/home/snRNAseq_analysis/singularity/snRNAseq_workflow.img'
  cpus = { check_max( 12 * task.attempt, 'cpus' ) }
  publishDir "${output_dir}/clustering/",
    mode: 'copy',
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"

  input:
    path rmd_file
    tuple val(patient), path(rds_file), val(output_dir)
    path params_file

  output:
    tuple val(patient), path('seu_annotated.rds'), val(output_dir), emit: ch_clustered
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
        cache_dir = "${output_dir}/clustering/clustering_cache/"),
      output_file = "clustering.html",
      output_dir = getwd()
    )
    """
}

workflow {

  // save params
  save_params()

  // get malignant files
  Channel
    .fromPath("${params.output.dir}/by_patient_wo_organoids/A/integrating/infercnv/infercnv_cache/*/seu_infercnv_malignant.rds")
    .set { ch_malignant_files }
  
  // one channel of all patients
  ch_malignant_files 
    .collect()
    .map { file -> tuple('malignant', file, "${params.output.dir}/by_patient_wo_organoids/malignant/") }
    .set { ch_run }
  
  // one channel per patient 
  ch_malignant_files
    .map { file -> 
            tuple(file.toString().tokenize('/')[-6], file) }
    .groupTuple()
    .map { patient, file -> tuple(
            patient, 
            file, 
            "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/malignant/")}
    .concat(ch_run)
    .set { ch_run }
  
  ch_run.view()
  
  // integrate
  integrating(
    "${baseDir}/templates/integrating.rmd",
    ch_run,
    save_params.out.ch_params
  )
  
  // cluster
  clustering(
    "${baseDir}/templates/seurat_clustering.rmd",
    integrating.out.ch_integrated,
    save_params.out.ch_params
  )
  
}

