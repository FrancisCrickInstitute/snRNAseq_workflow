#!/usr/bin/env nextflow

// default params set in nextflow.config

// import json function
import groovy.json.JsonOutput

// using DSL-2
nextflow.enable.dsl=2

log.info """\

snRNAseq workflow
=================
input dir    : ${params.input.dir}
output dir   : ${params.output.dir}
"""

// function which prints help message text
def help_message() {
    log.info"""
usage:

nextflow run <ARGUMENTS>

  input:
  --input.dir   directory containing sparse matrix of aligned snRNA-seq data
  
  output:
  --output.dir  output directory
  
    """.stripIndent()
}

process save_params {
  publishDir params.output.dir, mode: 'copy'
  
  output:
    path 'params.json'

  script:
    "echo '${JsonOutput.toJson(params)}'  > params.json"
      
}

// quality control
process quality_control {
  publishDir "$params.output.dir/quality_control/", mode: 'copy', pattern: "{*.html,*.rds}"
  publishDir "$params.output.dir/quality_control/", mode: 'copy', pattern: "*_files/figure-html/*.png"
  conda "environment.yml"
  cpus 2
  memory 20.GB
  time 2.hour

  input:
    path rmd_file
    path params_file

  output: 
    path 'seu_quality_controlled.rds', emit: ch_filtered
    path 'quality_control_report.html' 
    path 'quality_control_files/figure-html/*.png'
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}"), output_file = "quality_control_report.html", output_dir = getwd())'
    """
}

// clustering
process clustering {
  publishDir "$params.output.dir/clustering/", mode: 'copy', pattern: "{*.html,*.rds}"
  publishDir "$params.output.dir/clustering/", mode: 'copy', pattern: "clustering_files/figure-html/*.png"
  conda "environment.yml"
  cpus 2
  memory 20.GB
  time 2.hour

  input:
    path rmd_file
    path params_file
    path rds_file

  output: 
    path 'cds_clustered.rds', emit: ch_clustered
    path 'clustering_report.html' 
    path 'clustering_files/figure-html/*.png'
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}", rds_file = "${rds_file}"), output_file = "clustering_report.html", output_dir = getwd())'
    """
}

// celltype_annotation
process celltype_annotation {
  publishDir "$params.output.dir/celltype_annotation/", mode: 'copy', pattern: "{*.html,*.rds}"
  publishDir "$params.output.dir/celltype_annotation/", mode: 'copy', pattern: "*_files/figure-html/*.png"
  conda "environment.yml"
  cpus 2
  memory 20.GB
  time 2.hour

  input:
    path rmd_file
    path params_file
    path rds_file

  output: 
    path 'cds_celltype_annotated.rds', emit: ch_annotated
    path 'celltype_annotation_report.html' 
    path 'celltype_annotation_files/figure-html/*.png'
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}", rds_file = "${rds_file}"), output_file = "celltype_annotation_report.html", output_dir = getwd())'
    """
}

// infercnv
process infercnv {
  publishDir "$params.output.dir/infercnv/", mode: 'copy', pattern: "{*.html,*.rds}"
  publishDir "$params.output.dir/infercnv/", mode: 'copy', pattern: "*_files/figure-html/*.png"
  publishDir "$params.output.dir/infercnv/", mode: 'copy', pattern: "infercnv/*"
  conda "environment.yml"
  cpus 2
  memory 20.GB
  time 2.hour

  input:
    path rmd_file
    path params_file
    path rds_file

  output: 
    path 'cds_infercnv.rds', emit: ch_infercnv
    path 'infercnv_report.html' 
    path 'infercnv_files/figure-html/*.png'
    path 'infercnv/*'
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}", rds_file = "${rds_file}"), output_file = "infercnv_report.html", output_dir = getwd())'
    """
}

// main workflow
workflow {

  // show help message if the user specifies the --help flag at runtime
  // or if any required params are not provided
  if ( params.help || ! params.input.dir ){
      // invoke the function above which prints the help message
      help_message()
      // exit out and do not run anything else
      exit 1
  }
  
  // save params
  params_ch = save_params()
  
  // quality control and filtering
  quality_control("${baseDir}/templates/quality_control.rmd", params_ch)
  
  // dimensionality reduction and clustering
  clustering(
    "${baseDir}/templates/clustering.rmd", 
    params_ch, 
    quality_control.out.ch_filtered
  ) 
  
  // cell type annotation
  celltype_annotation(
    "${baseDir}/templates/celltype_annotation.rmd", 
    params_ch, 
    clustering.out.ch_clustered
  ) 
  
  // infercnv - run if reference celltypes provided
  if ( params.infercnv.reference_celltypes != false & params.infercnv.gene_order_file != false ) {
    infercnv(
      "${baseDir}/templates/infercnv.rmd",
      params_ch,
      celltype_annotation.out.ch_annotated
    )
  }
}