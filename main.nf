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
genome       : ${params.genome}
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

// render qc report
process render_qc_report {
  publishDir params.output.dir, mode: 'copy', pattern: '{*.html,*.tsv,qc_report_files/}'
  conda "environment.yml"

  input:
    path rmd_file
    path params_file

  output: 
    path 'qc_report.html' 

  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}"), output_file = "qc_report.html", output_dir = getwd())'
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
  
  // render QC report
  rmd_ch = file("${baseDir}/templates/qc_report.rmd")
  params_ch = save_params()
  render_qc_report(rmd_ch, params_ch) 

}