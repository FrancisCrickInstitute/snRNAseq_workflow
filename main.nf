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
  
  annotating celltypes:
  --markers_file                  file containing markers for celltype annotation, must contain 'gene' and 'population' columns
  --annotations_file              file containing final celltype annotations for the dataset, must contain an 'annotation' column and a column matching one of the metadata columns (e.g., 'cell', 'cluster' or 'partition') of the object
  
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
  tag "${id}"
  label 'process_low'
  publishDir "${params.output.dir}/${id}/", 
    mode: 'copy'
  conda "environment.yml"

    
  input:
    tuple val(id), val(id_col), val(dir)
    
  output:
    tuple val(id), path('seu.rds'), emit: ch_loaded
    
  script:
    """
    Rscript ${baseDir}/templates/load_input.R \
      "${id}" \
      "${id_col}" \
      "${dir}" \
      "${params.input.sample_metadata_file}"
    """
}

// merge samples
process merge_samples {
  tag "merged"
  label 'process_medium'
  publishDir "${params.output.dir}/merged/",
    mode: 'copy',
    pattern: "*.rds"
  
  input:
    path rds_files, stageAs: "seu??.rds"
    
  output:
    tuple val("merged"), path('seu.rds'), emit: ch_merged
    
  script:
    """
    Rscript ${baseDir}/templates/merge_samples.R \
      ${rds_files.join(',')}
    """
}

// filtering
process filtering {
  tag "${id}"
  label 'process_medium'
  publishDir "${params.output.dir}/${id}/filtering/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"
  conda "environment.yml"

  input:
    tuple val(id), path(rds_file)
    path rmd_file
    path params_file

  output: 
    tuple val(id), path('seu_filtered.rds'), emit: ch_filtered
    path 'filtering.html' 
    path 'filtering_files/figure-html/*.png'
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}", rds_file = "${rds_file}"), output_file = "filtering.html", output_dir = getwd())'
    """
}

// clustering
process clustering {
  tag "${id}"
  label "${ id = 'merged' ? 'process_high_memory' : 'process_medium' }"
  publishDir "${params.output.dir}/${id}/clustering/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"
  conda "environment.yml"
  
  input:
    tuple val(id), path(rds_file)
    path rmd_file
    path params_file

  output: 
    tuple val(id), path('cds_clustered.rds'), emit: ch_clustered
    path 'clustering.html' 
    path 'clustering_files/figure-html/*.png'
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}", rds_file = "${rds_file}"), output_file = "clustering.html", output_dir = getwd())'
    """
}

// celltype annotation
process annotating {
  tag "${id}"  
  label 'process_medium'
  publishDir "${params.output.dir}/${id}/annotating/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"
  conda "environment.yml"

  input:
    tuple val(id), path(rds_file)
    path rmd_file
    path params_file

  output: 
    tuple val(id), path('cds_celltype_annotated.rds'), emit: ch_annotated, optional: true
    path 'annotating.html' 
    path 'annotating_files/figure-html/*.png'
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}", rds_file = "${rds_file}"), output_file = "annotating.html", output_dir = getwd())'
    """
}

// infercnv
process infercnv {
  tag "${id}"
  label 'process_medium'
  publishDir "${params.output.dir}/${id}/infercnv/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png,infercnv/*}"
  conda "environment.yml"

  input:
    tuple val(id), path(rds_file)
    path rmd_file
    path params_file

  output: 
    path 'infercnv.html' 
    path 'infercnv_files/figure-html/*.png', optional: true
    path 'infercnv.rds', emit: ch_infercnv, optional: true
    path 'infercnv/*', optional: true
    
  script:
    """
    Rscript -e 'rmarkdown::render("${rmd_file}", params = list(params_file = "${params_file}", rds_file = "${rds_file}"), output_file = "infercnv.html", output_dir = getwd())'
    """
}

// pipeline workflow
workflow snRNAseq_workflow {
  take:
    ch_loaded
    ch_params
    
  main:
    
    // quality control and filtering
    filtering(
      ch_loaded,
      "${baseDir}/templates/filtering.rmd", 
      ch_params,
    )
    
    // dimensionality reduction and clustering
    clustering(
      filtering.out.ch_filtered,
      "${baseDir}/templates/clustering.rmd", 
      ch_params
    ) 
    
    // cell type annotation
    annotating(
      clustering.out.ch_clustered,
      "${baseDir}/templates/annotating.rmd", 
      ch_params
    ) 
    
    // infercnv - run if reference celltypes provided
    if ( params.annotating.annotations_file != false &
         params.infercnv.reference_celltypes != false & 
         params.infercnv.gene_order_file != false ) {
      infercnv(
        annotating.out.ch_annotated,
        "${baseDir}/templates/infercnv.rmd",
        ch_params
      )
    }

}

workflow {
  
  // show help message if the user specifies the --help flag at runtime
  // or if any required params are not provided
  if ( params.help || ! params.input.manifest_file || !params.input.sample_metadata_file ) {
      // invoke the function above which prints the help message
      help_message()
      // exit out and do not run anything else
      exit 1
  }
  
  // save params
  save_params()
    
  // generate one channel per id
  Channel
    .fromPath(params.input.manifest_file)
    .splitCsv(header:true, sep:'\t')
    .map { row -> tuple( row.id, row.id_col, row.dir ) }
    .set { ch_ids }
  
  // load each id
  load_input(
    ch_ids
  )
  
  // merge all ids
  merge_samples(
    load_input.out.ch_loaded
    .map { id, rds_file -> rds_file }
    .collect()
  )
  
  // run on each sample and on all samples
  snRNAseq_workflow(
    merge_samples.out.ch_merged
      .concat(load_input.out.ch_loaded),
    save_params.out.ch_params
  )
  
}
