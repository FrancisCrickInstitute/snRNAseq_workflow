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