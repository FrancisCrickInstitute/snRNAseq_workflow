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
  
  inputs:
  --input.manifest_file           directory containing sparse matrix of aligned snRNA-seq data
  --input.sample_metadata_file    file containing sample metadata, must contain a 'sample' column
  
  parallelisation:
  --input.run_each                run the pipeline on each sample
  --input.run_all                 run the pipeline on all samples
  --input.run_by_patient          run the pipeline on each patient (sample IDs must be in the form "patient_*")

  output:
  --output.dir                    optional. output directory, default is './output/'
  
  filters:
  --filter.do_filtering           logical, default is true
  --filter.n_mads                  integer, number of median absolute deviancies for outlier detection
  --filter.doublets               logical, default is true
  --filter.n_nuclei_per_gene.max  integer or 'adaptive', defailt is 5
  --filter.nCount_RNA.max         integer or 'adaptive', default is 10000
  --filter.nCount_RNA.min         integer or 'adaptive', default is 300
  --filter.nFeature_RNA.max       integer or 'adaptive', default is 10000
  --filter.nFeature_RNA.min       integer or 'adaptive', default is 200
  --filter.percent_mito.max       proportion (0-1) or 'adaptive', default is 0.1
  
  dimensionality reduction and clustering:
  --n_dims                        integer, number of PCs
  --vars_to_regress               vector of variables to regress out of the SCTransform residuals
  --clustering_resolution         resolution for clustering
  
  annotating celltypes:
  --annotate.markers_file         file containing markers for celltype annotation, must contain 'gene' and 'population' columns
  --annotate.annotations_file     file containing final celltype annotations for the dataset, must contain an 'annotation' column and a column matching one of the metadata columns (e.g., 'cell', 'cluster' or 'partition') of the object
  
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
  publishDir "${params.output.dir}/by_sample/${id}/", 
    mode: 'copy'
  
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

// filtering
process filtering {
  tag "${id}"
  label 'process_medium'
  publishDir "${params.output.dir}/by_sample/${id}/filtering/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png,nucleus_filtering.rds}"

  input:
    path rmd_file
    tuple val(id), path(rds_file)
    path params_file

  output: 
    tuple val(id), path('seu_filtered.rds'), emit: ch_filtered
    path 'seu_filters.rds', emit: ch_qc
    path 'filtering.html' 
    path 'filtering_files/figure-html/*.png'
    
  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}", 
      params = list(params_file = "${params_file}", rds_file = "${rds_file}"), 
      output_file = "filtering.html", 
      output_dir = getwd()
    )
    """
}

// merge samples
process merging {
  tag "${id}"
  label 'process_medium'
  publishDir "${params.output.dir}/${id}/",
    mode: 'copy',
    pattern: "*.rds"
  
  input:
    tuple val(id), path(rds_files, stageAs: "seu??.rds")
    
  output:
    tuple val(id), path('seu.rds'), emit: ch_merged
    
  script:
    """
    Rscript ${baseDir}/templates/merging.R \
      ${rds_files.join(',')}
    """
}

// integrate samples
process integrating {
  tag "${id}"
  label 'process_medium'
  publishDir "${params.output.dir}/${id}/",
    mode: 'copy',
    pattern: "*.rds"
  
  input:
    tuple val(id), path(rds_files, stageAs: "seu??.rds")
    
  output:
    tuple val(id), path('seu_integrated.rds'), emit: ch_integrated
    
  script:
    """
    Rscript ${baseDir}/templates/integrating.R \
      ${rds_files.join(',')}
    """
}

// clustering
process clustering {
  tag "${id}"
  label 'process_medium'
  publishDir "${params.output.dir}/${subdir}/${id}/clustering/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"

  input:
    path rmd_file
    tuple val(id), val(subdir), path(rds_file)
    path params_file

  output: 
    tuple val(id), val(subdir), path('cds_clustered.rds'), emit: ch_clustered
    path 'clustering.html' 
    path 'clustering_files/figure-html/*.png'
    
  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}", 
      params = list(params_file = "${params_file}", rds_file = "${rds_file}"), 
      output_file = "clustering.html",
      output_dir = getwd()
    )
    """
}

// celltype annotation
process annotating {
  tag "${id}"  
  label 'process_medium'
  publishDir "${params.output.dir}/${subdir}/${id}/annotating/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png,cell_groupings_summary.tsv,top_3_markers.tsv}"

  input:
    path rmd_file
    tuple val(id), val(subdir), path(rds_file)
    path params_file

  output: 
    tuple val(id), val(subdir), path('cds_celltype_annotated.rds'), emit: ch_annotated, optional: true
    path 'cds_singler_annotated.rds'
    path 'annotating.html' 
    path 'annotating_files/figure-html/*.png'
    path 'top_3_markers.tsv'
    
  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}", 
      params = list(params_file = "${params_file}", rds_file = "${rds_file}"), 
      output_file = "annotating.html", 
      output_dir = getwd()
    )
    """
}

// infercnv
process infercnv {
  tag "${id}"
  label 'process_medium'
  publishDir "${params.output.dir}/${id}/infercnv/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png,infercnv/*}"

  input:
    path rmd_file
    tuple val(id), path(rds_file)
    path params_file

  output: 
    path 'infercnv.html' 
    path 'infercnv_files/figure-html/*.png', optional: true
    path 'infercnv.rds', emit: ch_infercnv, optional: true
    path 'infercnv/*', optional: true
    
  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}", 
      params = list(params_file = "${params_file}", rds_file = "${rds_file}"), 
      output_file = "infercnv.html",
      output_dir = getwd()
    )
    """
} 

// pipeline workflow
workflow snRNAseq_workflow {
  take:
    ch_filtered
    ch_params
    
  main:
    
    // dimensionality reduction and clustering
    clustering(
      "${baseDir}/templates/clustering.rmd",
      ch_filtered,
      ch_params
    ) 
    
    // cell type annotation
    annotating(
      "${baseDir}/templates/annotating.rmd",
      clustering.out.ch_clustered,
      ch_params
    )

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
  
  // quality control and filtering
  filtering(
    "${baseDir}/templates/filtering.rmd",
    load_input.out.ch_loaded,
    save_params.out.ch_params,
  )
  
  // initiate ch_run, ch_integrate and ch_merge
  Channel
    .empty()
    .set { ch_run }
  
  if ( params.input.run_all ) {
    // integrate all samples
    filtering.out.ch_filtered
      .map { id, rds_file -> tuple("integrated", "integrated", rds_file) }
      .groupTuple()
      .set { ch_run_all }
      
    // perform integration
    integrating(
      ch_run_all
    )
    
    // add integrated samples to ch_run
    ch_run
      .concat(integrating.out.ch_integrated)
      .set { ch_run }
  }
  
  if ( params.input.run_by_patient ) {
    // merge all samples by patient
    filtering.out.ch_filtered
      .map { id, rds_file -> tuple(id.split("_")[0], "by_patient", rds_file) }
      .groupTuple()
      .set { ch_run_by_patient }
      
    // perform merge
    merging(
      ch_run_by_patient
    )
    
    // add merged samples to ch_run
    ch_run
      .concat(merging.out.ch_merged)
      .set { ch_run }
  }
  
  if ( params.input.run_each ) {
    // get each sample
    filtering.out.ch_filtered
      .map { id, rds_file -> tuple(id, "by_sample", rds_file) }
      .set { ch_run_each }
      
    // add each sample to ch_run
    ch_run
      .concat(ch_run_each)
      .set { ch_run }
  }
  
  // run on each/all samples/patients
  snRNAseq_workflow(
    ch_run,
    save_params.out.ch_params
  )
  
}



  
  // infercnv - run if reference celltypes or samples provided, on merged output
  // TODO: make this work on annotated output for infercnv.reference_samples or
  // clustered output for infercnv.reference_celltypes and work on all merged samples
  // or all patient samples or annotated individual samples!
  //if (  params.infercnv.gene_order_file != false &
  //      ((params.infercnv.reference_samples != false) ||
  //      (params.infercnv.reference_celltypes != false & params.annotate.annotations_file != false))
  //   ) {
  //  infercnv(
  //    // merged samples will be the first out
  //    ch_merged_annotated,
  //    "${baseDir}/templates/infercnv.rmd",
  //    ch_params
  //  )
  //}

