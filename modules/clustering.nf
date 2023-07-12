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
    tuple val(id), val(subdir), path('cds.rds'), emit: ch_clustered
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