#!/usr/bin/env nextflow

// default params set in nextflow.config

// import json function
import java.io.BufferedReader

// using DSL-2
nextflow.enable.dsl=2

log.info """\

snRNAseq workflow
=================
manifest file:        ${params.input.manifest_file}
sample metadata file: ${params.input.sample_metadata_file}
output dir:           ${params.output.dir}
"""

// import modules
include { save_params } from './modules/save_params'

// perform infercnv on each sample
process infercnv {
  tag "${sample}"
  memory { 200.GB * task.attempt }
  cpus { 12 * task.attempt }
  time { 24.hour * task.attempt }
  queue { 200 * task.attempt < 921 ? 'v1_medium72' : 'v1_largemem72' }

  publishDir "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/infercnv_cache/${sample}/",
    mode: 'copy',
    pattern: "{*.rds,*.html}"

  input:
    path rmd_file
    tuple val(sample), val(patient), val(rds_file)

  output:
    tuple val(patient), path('seu_infercnv.rds'), val(''), emit: ch_infercnv, optional: true
    tuple val(patient), path('seu_infercnv_malignant.rds'), val('malignant'), emit: ch_infercnv_malig, optional:true
    path 'infercnv.html'

  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        sample = "${sample}",
        rds_file = "${rds_file}",
        gene_order_file = "${params.infercnv.gene_order_file}",
        cache_dir = "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/infercnv_cache/${sample}/"),
      output_file = "infercnv.html",
      output_dir = getwd()
    )
    """
}

// integrate malignant samples
process integrating {
  tag "${patient}_${subdir}"
  label 'process_medium'
  publishDir "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/${subdir}/integrating/",
    mode: 'copy',
    pattern: "*.rds"

  input:
    path rmd_file
    tuple val(patient), path(rds_files, stageAs: "seu??.rds"), val(subdir)
    path params_file

  output:
    tuple val(patient), path('seu.rds'), val(subdir), emit: ch_integrated
    path 'integrating.html'

  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        params_file = "${params_file}",
        rds_file = "${rds_files.join(',')}",
        cache_dir = "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/${subdir}/integrating/integrating_cache/"),
      output_file = "integrating.html",
      output_dir = getwd()
    )
    """
}

// cluster malignant samples
process clustering {
  tag "${patient}"
  label 'process_medium'
  container '/rds/general/user/art4017/home/snRNAseq_analysis/singularity/snRNAseq_workflow.img'
  cpus { check_max( 12 * task.attempt, 'cpus' ) }
  publishDir "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/malignant/clustering/",
    mode: 'copy',
    pattern: "{*.html,*.rds,*_files/figure-html/*.png}"

  input:
    path rmd_file
    tuple val(patient), path(rds_file)
    path params_file

  output:
    tuple val(patient), path('seu_annotated.rds'), emit: ch_clustered
    path 'clustering.html'
    path 'clustering_files/figure-html/*.png'

  script:
    """
    #!/usr/bin/env Rscript
    rmarkdown::render(
      "${rmd_file}",
      params = list(
        params_file = "${params_file}",
        rds_file = "${rds_file}",
        cache_dir = "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/infercnv/malignant/clustering/clustering_cache/"),
      output_file = "clustering.html",
      output_dir = getwd()
    )
    """
}

workflow {

  // save params
  save_params()

  // generate one channel per id (exclude organoids + J_post_T1_biopsy, throwing unknown error)
  //Channel
    //.fromPath(params.input.manifest_file)
    //.splitCsv(header:true, sep:'\t')
    //.map { row -> tuple(
    //        row.id, row.id.split("_")[0],
    //         "${params.output.dir}/by_patient_wo_organoids/"+row.id.split("_")[0]+"/integrating/seurat_clustering/seu_annotated.rds") } 
    //.filter { !it[0].contains("J_post_T1_biopsy") }    
    //.filter { !it[0].contains("organoid") } 
    //.set { ch_run }

  Channel
    .fromList([ 'B_post_T1_biopsy',
                'F_post_N2_biopsy',
                'H_post_T1_biopsy',
                'I_pre_T0_biopsy',
                'N_post_T1_biopsy',
                'P_post_T1_biopsy'])
    .map { id -> tuple(
	    id, id.split("_")[0],
	     "${params.output.dir}/by_patient_wo_organoids/"+id.split("_")[0]+"/integrating/seurat_clustering/seu_annotated.rds") }     
    .set { ch_run }

  // run infercnv
  infercnv(
    "${baseDir}/templates/infercnv.rmd",
    ch_run
  )
  
}


//// integrate all preloaded samples
//Channel
//  .fromPath("${params.output.dir}/by_patient_wo_organoids/*/integrating/seurat_clustering/seu_annotated.rds")
//  .set { ch_run }

//// get id, dir, subdir and filepath
//ch_run
//  .map { file ->
//          tuple(file.toString().tokenize('/')[-4],
//                file.toString().tokenize('/')[-5],
//                file.toString().tokenize('/')[-3],
//                file) }
//  .set { ch_run }

//// perform infercnv
//infercnv(
//  "${baseDir}/templates/infercnv.rmd",
//  ch_run,
//  save_params.out.ch_params
//)

//// perform malignant integration
//malignant_integrating(
//  "${baseDir}/templates/integrating.rmd",
//  infercnv.out.ch_malignant,
//  save_params.out.ch_params,
//  "${params.output.dir}/${dir}/${id}/integrating/infercnv/malignant_integrating/"
//)

//// perform malignant clustering
//malignant_clustering(
//  "${baseDir}/templates/seurat_clustering.rmd",
//  malignant_integrating.out.ch_integrated,
//  save_params.out.ch_params,
//  "${params.output.dir}/${dir}/${id}/integrating/infercnv/malignant_clustering/"
//)
