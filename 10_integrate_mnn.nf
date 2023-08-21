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

// integrate with mnn
process integrating_mnn {
  tag "${patient}"
  label 'process_medium'
  
  publishDir "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating_mnn/",
    mode: 'copy',
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"

  input:
    path rmd_file
    tuple val(patient), path(rds_file)

  output:
    path 'seu.rds'
    path 'integrating_mnn.html'
    path 'integrating_mnn_files/figure-html/*.png', optional: true

  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        cache_dir = "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating_mnn/integrating_mnn_cache/",
        infercnv_cache_dir = "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/infercnv_cache/",
        malignancy_score_file = "${params.annotate.malignancy_score_file}"),
      output_file = "integrating_mnn.html",
      output_dir = getwd()
    )
    """
}

workflow {

  // by patient, wo organoids
  Channel
    .fromPath("${params.output.dir}/by_patient_wo_organoids/*/integrating/seurat_clustering/seu_annotated.rds")
    .map { rds_file -> tuple(rds_file.toString().tokenize('/')[-4].split("_")[0], rds_file) }
    .set { ch_run }
    
  // all malignant cells, wo organoids
  Channel
    .fromPath("${params.output.dir}/by_patient_wo_organoids/malignant/clustering/seu_annotated.rds")
    .map { rds_file -> tuple('malignant', rds_file) }
    .concat(ch_run)
    .set { ch_run }
    
  // all cells
  Channel
    .fromPath("${params.output.dir}/by_patient_wo_organoids/*/integrating/seurat_clustering/seu_annotated.rds")
    .map { rds_file -> tuple('all', rds_file) }
    .groupTuple()
    .concat(ch_run)
    .set { ch_run }
  
  // perform MNN integration
  integrating_mnn(
    "${baseDir}/templates/integrating_mnn.rmd",
    ch_run
  )

}
