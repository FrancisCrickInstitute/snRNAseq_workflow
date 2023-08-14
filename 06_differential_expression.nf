#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// differential expression analysis
process diffexp {
  tag "${patient}"
  queue { 50 * task.attempt < 921 ? 'v1_medium72' : 'v1_largemem72' }
  memory { 50.GB * task.attempt }
  cpus { 6 * task.attempt }
  time { 24.hour * task.attempt }
  publishDir "${output_dir}",
    mode: 'copy',
    pattern: "{*.rds,*.html}"
    
  input:
    path(rmd_file)
    tuple val(patient), path(rds_file), val(infercnv_cache_dir), val(output_dir)
    
  output:
    path 'differential_expression.html'
    path 'aneuploid_cells.rds', optional: true
    path 'malignant_cells.rds', optional: true
    path 'degs.rds', optional: true
    path 'scpa.rds', optional: true
    
  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        cache_dir = "${output_dir}/differential_expression_cache/",
        infercnv_cache_dir = "${infercnv_cache_dir}"),
      output_file = "differential_expression.html",
      output_dir = getwd()
    )
    """
}

workflow {
  
  // get seu_annotated files
  Channel
    .fromPath("${params.output.dir}/by_patient_wo_organoids/*/integrating/seurat_clustering/seu_annotated.rds")
    .map { rds_file -> tuple(rds_file.toString().tokenize('/')[-4], rds_file)}
    .map { patient, rds_file -> tuple(
      patient, rds_file, 
      "${params.output.dir}/by_patient_wo_organoids/"+patient+"/integrating/infercnv/infercnv_cache/",
      "${params.output.dir}/by_patient_wo_organoids/"+patient+"/integrating/differential_expression/")}
    .set { ch_run }
  
  // differential expression
  diffexp(
    "${baseDir}/templates/differential_expression.rmd",
    ch_run
  )
  
}

