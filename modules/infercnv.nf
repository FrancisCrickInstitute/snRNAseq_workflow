// infercnv
process infercnv {
  tag "${id}"
  label 'process_high_memory'
  publishDir "${params.output.dir}/${dir}/${id}/${subdir}/infercnv/",
    mode: 'copy',
    pattern: "{*.html,*.tsv,*.rds,*_files/figure-html/*.png,infercnv/*}"

  input:
    path rmd_file
    tuple val(id), val(dir), val(subdir), path(rds_file)
    path params_file

  output:
    path 'infercnv.html'
    path 'infercnv_annots.tsv'
    path 'infercnv_files/figure-html/*.png', optional: true
    path 'infercnv.rds', emit: ch_infercnv, optional: true
    path 'infercnv/*', optional: true

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
