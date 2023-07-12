// merge samples
process merging {
  tag "${id}"
  label 'process_medium'
  publishDir "${params.output.dir}/${subdir}/${id}/merging/",
    mode: 'copy',
    pattern: "*.rds"
  
  input:
    tuple val(id), val(subdir), path(rds_files, stageAs: "seu??.rds")
    
  output:
    tuple val(id), val(subdir), path('seu.rds'), emit: ch_merged
    
  script:
    """
    Rscript ${baseDir}/templates/merging.R \
      ${rds_files.join(',')}
    """
}