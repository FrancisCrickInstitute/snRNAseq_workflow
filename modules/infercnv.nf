// infercnv
process infercnv {
  tag "${id}"
  label 'process_high_memory'
  publishDir "${params.output.dir}/${dir}/${id}/${subdir}/infercnv/",
    mode: 'copy',
    pattern: "{*.rds,*.html,*/*}"

  input:
    path rmd_file
    tuple val(id), val(dir), val(subdir), path(rds_file)
    path params_file

  output:
    path 'infercnv.html'
    path 'infercnv.rds', optional: true
    tuple val(id), val(dir), val(subdir), path('seu_infercnv_malignant.rds'), emit: ch_malignant, optional: true
    
  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        params_file = "${params_file}",
        rds_file = "${rds_file}",
        cache_dir = "${params.output.dir}/${dir}/${id}/${subdir}/infercnv/infercnv_cache/"),
      output_file = "infercnv.html",
      output_dir = getwd()
    )
    """
}
