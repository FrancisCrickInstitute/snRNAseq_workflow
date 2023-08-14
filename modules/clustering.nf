// clustering
process clustering {
  tag "${id}"
  label 'process_medium'
  publishDir "${params.output.dir}/${dir}/${id}/${subdir}/clustering/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png,*.tsv}"

  input:
    path rmd_file
    tuple val(id), val(dir), val(subdir), path(rds_file)
    path params_file

  output: 
    tuple val(id), val(dir), val(subdir), path('cds.rds'), emit: ch_clustered
    path 'clustering.html' 
    path 'clustering_files/figure-html/*.png'
    path 'cluster_markers', optional: true
    
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