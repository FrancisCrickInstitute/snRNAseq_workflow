// load input
process loading {
  tag "${id}"
  label 'process_low'
  publishDir "${params.output.dir}/by_sample/${id}/loading/", 
    mode: 'copy'
  
  input:
    tuple val(id), val(id_col), val(dir)
    
  output:
    tuple val(id), path('seu.rds'), emit: ch_loaded
    
  script:
    """
    Rscript ${baseDir}/templates/loading.R \
      "${id}" \
      "${id_col}" \
      "${dir}" \
      "${params.input.sample_metadata_file}"
    """
}