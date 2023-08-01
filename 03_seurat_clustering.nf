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
include { save_params        } from './modules/save_params'
include { seurat_clustering  } from './modules/seurat_clustering'

workflow {
  
  // save params
  save_params()
  
  // initiate ch_run
  Channel
    .empty()
    .set { ch_run }
  
  // get subdir
  subdir = params.input.integrate ? "integrating" : "merging"
  
  // get ch_run(s)
  if ( params.input.run_all ) {
    Channel
      .fromPath("${params.output.dir}/all/*/${subdir}/seu.rds")
      .map { file -> tuple("${subdir}", file)}
      .concat(ch_run)
      .set { ch_run }
  }
  if ( params.input.run_by_patient ) {
    Channel
    .fromPath("${params.output.dir}/by_patient/*/${subdir}/seu.rds")
    .map { file -> tuple("${subdir}", file)}
    .concat(ch_run)
    .set { ch_run }
  }
  if ( params.input.run_by_patient_wo_organoids ) {
    Channel
      .fromPath("${params.output.dir}/by_patient_wo_organoids/*/${subdir}/seu.rds")
      .map { file -> tuple("${subdir}", file)}
      .concat(ch_run)
      .set { ch_run }
  }
  
  // by sample
  if ( params.input.run_by_sample ) {
    Channel
      .fromPath("${params.output.dir}/by_sample/*/filtering/seu.rds")
      .map { file -> tuple("", file)}
      .concat(ch_run)
      .set { ch_run }
  }
  
  // get id, dir, subdir and filepath
  ch_run
    .map { subdir, file -> 
            tuple(file.toString().tokenize('/')[-3], file.toString().tokenize('/')[-4], subdir, file) }
    .set { ch_run }
  
  // dimensionality reduction and clustering
    seurat_clustering(
      "${baseDir}/templates/seurat_clustering.rmd",
      ch_run,
      save_params.out.ch_params
    ) 
  
}

