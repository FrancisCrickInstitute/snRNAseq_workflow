// integrate samples
process integrating {
  tag "${id}"
  label 'process_high_memory'
  publishDir "${params.output.dir}/${subdir}/integrating/",
    mode: 'copy',
    pattern: "*.rds"
  
  input:
    tuple val(id), val(subdir), path(rds_files, stageAs: "seu??.rds")
    
  output:
    tuple val(id), val(subdir), path('seu.rds'), emit: ch_integrated
    
  script:
    """
    Rscript ${baseDir}/templates/integrating.R \
      ${rds_files.join(',')}
    """
}