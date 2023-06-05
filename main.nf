#!/usr/bin/env nextflow

// default params set in nextflow.config

// import json function
import groovy.json.JsonOutput

// using DSL-2
nextflow.enable.dsl=2

log.info """\

snRNAseq workflow
=================
data dir     : ${params.data_dir}
output dir   : ${params.output_dir}
genome       : ${params.genome}
"""

// function which prints help message text
def help_message() {
    log.info"""
usage:

nextflow run <ARGUMENTS>

  data:
  --data_dir            directory containing sparse matrix of aligned snRNA-seq data
  
  output:
  --output_dir          output directory
  
    """.stripIndent()
}

process save_params {
  publishDir params.output_dir, mode: 'copy'
  
  output:
    path 'params.json'

  script:
    "echo '${JsonOutput.toJson(params)}' |  > params.json"
      
}

// check snRNAseq platform
process detect_platform {
  
  input:
    path data_dir
  
  output:
    stdout
    
  script:
    """
    (
      cd '${data_dir}'
      if [ -f "matrix.mtx" -a -f "genes.tsv" -a -f "barcodes.tsv" ] ; then
          echo ; echo "10X Genomics platform detected"
          platform="10X Genomics"
      elif [ -f "DGE.mtx" -a -f "all_genes.csv" -a -f "cell_metadata.csv" ] ; then
          echo ; echo "Parse Biosciences platform detected"
          platform="Parse Biosciences"
      fi
    )
    """
}

// render qc report
process render_qc_report {
  publishDir params.output_dir, mode: 'copy', pattern: '*.html'
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
  if ( params.help || ! params.data_dir ){
      // invoke the function above which prints the help message
      help_message()
      // exit out and do not run anything else
      exit 1
  }
  
  // detect platform (Parse or 10X)
  data_dir_ch = Channel.fromPath(params.data_dir)
  platform_ch = detect_platform(data_dir_ch) 
  platform_ch.view { it }
  
  // render QC report
  rmd_ch = file("${baseDir}/templates/qc_report.rmd")
  params_ch = save_params()
  render_qc_report(rmd_ch, params_ch)

}