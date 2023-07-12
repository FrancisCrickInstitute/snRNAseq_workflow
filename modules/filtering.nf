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
    tuple val(id), path('seu.rds'), emit: ch_filtered
    path 'filters.rds', emit: ch_qc
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