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