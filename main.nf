#!/usr/bin/env nextflow

// default params set in nextflow.config

// import json function
import groovy.json.JsonOutput
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

// function which prints help message text
def help_message() {
    log.info"""
usage:

nextflow run <ARGUMENTS>
  
  (required)
  
  --input.manifest_file            directory containing sparse matrix of aligned snRNA-seq data
  --input.sample_metadata_file    file containing sample metadata, must contain a 'sample' column
  
  (optional)
  
  input:
  --input.sample_subset           comma-delimited list of samples to use
  --input.cell_subset             comma-delimited list of cells to use

  output:
  --output.dir                    optional. output directory, default is './output/'
  
  filters:
  --filter.do_filtering           logical, default is true
  --filter.doublets               logical, default is true
  --filter.n_nuclei_per_gene.max  integer, defailt is 5
  --filter.nCount_RNA.max         integer, default is 10000
  --filter.nCount_RNA.min         integer, default is 300
  --filter.nFeature_RNA.max       integer, default is 10000
  --filter.nFeature_RNA.min       integer, default is 200
  --filter.percent_mito.max       proportion (0-1), default is 0.1
  
  dimensionality reduction and clustering:
  --n_dims                        integer, number of PCs
  --vars_to_regress               vector of variables to regress out of the SCTransform residuals
  
  celltype_annotation:
  --markers_file                  file containing markers for celltype annotation, must contain 'gene' and 'population' columns
  --celltype_annotations_file     file containing final celltype annotations for the dataset, must contain an 'annotation' column and a column matching one of the metadata columns (e.g., 'cell', 'cluster' or 'partition') of the object
  
  infercnv:
  --infercnv.reference_celltypes  comma-delimited list of annotations to use as a reference
  --infercnv.gene_order_file      file containing the positions of each gene along each chromosome in the genome
 
    """.stripIndent()
}

process save_params {
  publishDir "${params.output.dir}", mode: 'copy'
  
  output:
    path 'params.json', emit: ch_params

  script:
    "echo '${JsonOutput.toJson(params)}'  > params.json"
      
}

// load input
process load_input {
  tag "${sample}"
  publishDir "${params.output.dir}/${sample}/", 
    mode: 'copy'
  cpus 1
  memory 64.GB
  time 1.hour
    
  input:
    tuple val(sample), path(dir)
    
  output:
    path 'seu.rds', emit: ch_input
    
  script:
    """
    Rscript ${baseDir}/templates/load_input.R \
      "${sample}" \
      "${dir}" \
      "${params.input.sample_metadata_file}"
    """
}

// filtering
process filtering {
  tag "${sample}"
  publishDir "${params.output.dir}/${sample}/filtering/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"
  conda "environment.yml"
  cpus 2
  memory 40.GB
  time 2.hour

  input:
    tuple val(sample), path(dir)
    path rmd_file
    path params_file
    path rds_file

  output: 
    path 'seu_filtered.rds', emit: ch_filtered
    path 'filtering.html' 
    path 'filtering_files/figure-html/*.png'
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}", rds_file = "${rds_file}"), output_file = "filtering.html", output_dir = getwd())'
    """
}

// clustering
process clustering {
  tag "${sample}"
  publishDir "${params.output.dir}/${sample}/clustering/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"
  conda "environment.yml"
  cpus 2
  time 2.hour
  memory { 50.GB * task.attempt }
  maxRetries 4
  errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'finish' }  

  input:
    tuple val(sample), path(dir)
    path rmd_file
    path params_file
    path rds_file

  output: 
    path 'cds_clustered.rds', emit: ch_clustered
    path 'clustering.html' 
    path 'clustering_files/figure-html/*.png'
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}", rds_file = "${rds_file}"), output_file = "clustering.html", output_dir = getwd())'
    """
}

// celltype annotation
process celltype_annotation {
  tag "${sample}"
  publishDir "${params.output.dir}/${sample}/celltype_annotation/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"
  conda "environment.yml"
  cpus 2
  memory 20.GB
  time 2.hour

  input:
    tuple val(sample), path(dir)
    path rmd_file
    path params_file
    path rds_file

  output: 
    path 'cds_celltype_annotated.rds', emit: ch_annotated, optional: true
    path 'celltype_annotation.html' 
    path 'celltype_annotation_files/figure-html/*.png'
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}", rds_file = "${rds_file}"), output_file = "celltype_annotation.html", output_dir = getwd())'
    """
}

// infercnv
process infercnv {
  tag "${sample}"
  publishDir "${params.output.dir}/${sample}/infercnv/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png,infercnv/*}"
  conda "environment.yml"
  cpus 2
  memory 60.GB
  time 10.hour

  input:
    tuple val(sample), path(dir)
    path rmd_file
    path params_file
    path rds_file

  output: 
    path 'infercnv.rds', emit: ch_infercnv
    path 'infercnv.html' 
    path 'infercnv_files/figure-html/*.png'
    path 'infercnv/*'
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}", rds_file = "${rds_file}"), output_file = "infercnv.html", output_dir = getwd())'
    """
}

// merge samples
process merge_samples {
  publishDir "${params.output.dir}/merge/",
    mode: 'copy',
    pattern: "*.rds"
  
  input:
    path rmd_files, stageAs: "?/*"
    
  output:
    path 'seu.rds', emit: ch_merge
    
  script:
    """
    Rscript ${baseDir}/templates/merge_samples.R \
      ${rmd_files.join(',')}
    """
}

// workflow
workflow {
  
  // show help message if the user specifies the --help flag at runtime
  // or if any required params are not provided
  if ( params.help || ! params.input.manifest_file || !params.input.sample_metadata_file ){
      // invoke the function above which prints the help message
      help_message()
      // exit out and do not run anything else
      exit 1
  }
  
  // generate one channel per sample
  Channel
    .fromPath(params.input.manifest_file)
    .splitCsv(header:true, sep:'\t')
    .map { row -> tuple( row.sample, row.dir ) }
    .set { ch_input }

  // save params
  save_params()
  
  // load input
  load_input(
    ch_input
  )
  
  // quality control and filtering
  filtering(
    ch_input,
    "${baseDir}/templates/filtering.rmd", 
    save_params.out.ch_params,
    load_input.out.ch_input
  )
  
  // dimensionality reduction and clustering
  clustering(
    ch_input,
    "${baseDir}/templates/clustering.rmd", 
    save_params.out.ch_params, 
    filtering.out.ch_filtered
  ) 
  
  // cell type annotation
  celltype_annotation(
    ch_input,
    "${baseDir}/templates/celltype_annotation.rmd", 
    save_params.out.ch_params, 
    clustering.out.ch_clustered
  ) 
  
  // infercnv - run if reference celltypes provided
  if ( params.infercnv.reference_celltypes != false & params.infercnv.gene_order_file != false ) {
    infercnv(
      ch_input,
      "${baseDir}/templates/infercnv.rmd",
      save_params.out.ch_params,
      celltype_annotation.out.ch_annotated
    )
  }
  
  // merge samples
  merge_samples (
    filtering.out.ch_filtered.collect()
  )
  
}
