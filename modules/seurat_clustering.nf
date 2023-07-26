// seurat clustering
process seurat_clustering {
  tag "${id}"
  label 'process_high'
  cpus = { check_max( 12 * task.attempt, 'cpus' ) }
  publishDir "${params.output.dir}/${dir}/${id}/${subdir}/seurat_clustering/", 
    mode: 'copy', 
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"

  input:
    path rmd_file
    tuple val(id), val(dir), val(subdir), path(rds_file)
    path params_file

  output: 
    tuple val(id), val(dir), val(subdir), path('seu.rds'), emit: ch_clustered
    path 'seurat_clustering.html' 
    path 'seurat_clustering_files/figure-html/*.png'
    path 'cds_singler.rds'
    
  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}", 
      params = list(
        params_file = "${params_file}", 
        rds_file = "${rds_file}",
        cache_dir = "${baseDir}/${params.output.dir}/${dir}/${id}/${subdir}/seurat_clustering/seurat_clustering_cache/"), 
      output_file = "seurat_clustering.html",
      output_dir = getwd()
    )
    """
}