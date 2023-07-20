// merge samples
process merging {
  tag "${id}"
  label 'process_medium'
  publishDir "${params.output.dir}/${dir}/${id}/merging/",
    mode: 'copy',
    pattern: "*.rds"
  
  input:
    tuple val(id), val(dir), path(rds_files, stageAs: "seu??.rds")
    
  output:
    tuple val(id), val(dir), val("merging"), path('seu.rds'), emit: ch_merged
    
  script:
    """
    Rscript ${baseDir}/templates/merging.R \
      ${rds_files.join(',')}
    """
}