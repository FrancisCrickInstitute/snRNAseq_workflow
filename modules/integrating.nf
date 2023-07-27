// integrate samples
process integrating {
  tag "${id}"
  label 'process_high_memory'
  publishDir "${params.output.dir}/${dir}/${id}/integrating/",
    mode: 'copy',
    pattern: "*.rds"

  input:
    path rmd_file
    tuple val(id), val(dir), path(rds_files, stageAs: "seu??.rds")
    path params_file

  output:
    tuple val(id), val(dir), val("integrating"), path('seu.rds'), emit: ch_integrated
    path 'integrating.html'
    path 'integrating_files/figure-html/*.png'

  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        params_file = "${params_file}",
        rds_file = "${rds_files.join(',')}",
        cache_dir = "${params.output.dir}/${dir}/${id}/integrating/integrating_cache/"),
      output_file = "integrating.html",
      output_dir = getwd()
    )
    """
}
