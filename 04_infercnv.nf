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

// perform infercnv on each sample
process infercnv {
  tag "${sample}"
  label 'process_medium'
  publishDir "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/infercnv_cache/${sample}/",
    mode: 'copy',
    pattern: "*.rds"

  input:
    path rmd_file
    tuple val(sample), val(patient), val(rds_file)

  output:
    tuple val(sample), val(patient), path('seu_infercnv.rds'), emit: ch_infercnv, optional: true
    path 'infercnv.html'

  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        sample = "${sample}",
        rds_file = "${rds_file}",
        gene_order_file = "${params.infercnv.gene_order_file}",
        cache_dir = "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/infercnv_cache/${sample}/"),
      output_file = "infercnv.html",
      output_dir = getwd()
    )
    """
}

// integrate malignant samples
process malignant_integrating {
  tag "${id}"
  label 'process_medium'
  publishDir "${output_dir}",
    mode: 'copy',
    pattern: "*.rds"

  input:
    path rmd_file
    tuple val(id), val(dir), val(subdir), path(rds_files, stageAs: "seu??.rds")
    path params_file
    val output_dir

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
        cache_dir = "${output_dir}/malignant_integrating_cache/"),
      output_file = "integrating.html",
      output_dir = getwd()
    )
    """
}

// cluster malignant samples
process malignant_clustering {
  tag "${id}"
  label 'process_medium'
  container '/rds/general/user/art4017/home/snRNAseq_analysis/singularity/snRNAseq_workflow.img'
  cpus = { check_max( 12 * task.attempt, 'cpus' ) }
  publishDir "${output_dir}",
    mode: 'copy',
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"

  input:
    path rmd_file
    tuple val(id), val(dir), val(subdir), path(rds_file)
    path params_file
    val output_dir

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
        cache_dir = "${output_dir}/malignant_clustering_cache/"),
      output_file = "malignant_clustering.html",
      output_dir = getwd()
    )
    """
}

workflow {

  // generate one channel per id
  Channel
    .fromPath(params.input.manifest_file)
    .splitCsv(header:true, sep:'\t')
    .map { row -> tuple(
              row.id, 
              row.id.split("_")[0], 
              "${params.output.dir}/by_patient_wo_organoids/"+row.id.split("_")[0]+"/integrating/seurat_clustering/seu_annotated.rds") }
    .set { ch_run }

  // run infercnv
  infercnv(
    "${baseDir}/templates/infercnv.rmd",
    ch_run
  )

  //// integrate all preloaded samples
  //Channel
  //  .fromPath("${params.output.dir}/by_patient_wo_organoids/*/integrating/seurat_clustering/seu_annotated.rds")
  //  .set { ch_run }

  //// get id, dir, subdir and filepath
  //ch_run
  //  .map { file ->
  //          tuple(file.toString().tokenize('/')[-4],
  //                file.toString().tokenize('/')[-5],
  //                file.toString().tokenize('/')[-3],
  //                file) }
  //  .set { ch_run }

  //// perform infercnv
  //infercnv(
  //  "${baseDir}/templates/infercnv.rmd",
  //  ch_run,
  //  save_params.out.ch_params
  //)
  
  //// perform malignant integration
  //malignant_integrating(
  //  "${baseDir}/templates/integrating.rmd",
  //  infercnv.out.ch_malignant,
  //  save_params.out.ch_params,
  //  "${params.output.dir}/${dir}/${id}/integrating/infercnv/malignant_integrating/"
  //)
  
  //// perform malignant clustering
  //malignant_clustering(
  //  "${baseDir}/templates/seurat_clustering.rmd",
  //  malignant_integrating.out.ch_integrated,
  //  save_params.out.ch_params,
  //  "${params.output.dir}/${dir}/${id}/integrating/infercnv/malignant_clustering/"
  //)
  
}
