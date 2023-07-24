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
include { loading     } from './modules/loading'
include { filtering   } from './modules/filtering'
include { merging     } from './modules/merging'
include { integrating } from './modules/integrating'

workflow {
  
  // save params
  save_params()

  // generate one channel per id
  Channel
    .fromPath(params.input.manifest_file)
    .splitCsv(header:true, sep:'\t')
    .map { row -> tuple( row.id, row.id_col, row.dir ) }
    .set { ch_ids }
  
  // load each id
  loading(
    ch_ids
  )
  
  // quality control and filtering
  filtering(
    "${baseDir}/templates/filtering.rmd",
    loading.out.ch_loaded,
    save_params.out.ch_params,
  )
  
  // initiate ch_run
  Channel
    .empty()
    .set { ch_run }
  
  //get runs
  
  if ( params.input.run_by_sample ) {
    // by sample
    filtering.out.ch_filtered
      .map { id, rds_file -> tuple(id, "by_sample", rds_file)}
      .concat(ch_run)
      .set { ch_run }
  }
  
  if ( params.input.run_by_patient ) {
    // by patient
    filtering.out.ch_filtered
      .map { id, rds_file -> tuple(id.split("_")[0], rds_file) }
      .groupTuple()
      .map { id, rds_file -> tuple(id, "by_patient", rds_file) }
      .concat(ch_run)
      .set { ch_run }
  }
  
  if ( params.input.run_by_patient_wo_organoids ) {
    // by patient without organoids
    filtering.out.ch_filtered
      .filter { !it[0].contains("organoid") }
      .map { id, rds_file -> tuple(id.split("_")[0], rds_file) }
      .groupTuple()
      .map { id, rds_file -> tuple(id, "by_patient_wo_organoids", rds_file) }
      .concat(ch_run)
      .set { ch_run }
  }
    
  if ( params.input.run_all ) {
    // all
    filtering.out.ch_filtered
      .map { id, rds_file -> tuple("all", "all", rds_file) }
      .groupTuple()
      .concat(ch_run)
      .set { ch_run }
  }
  
  if ( ! params.input.integrate ) {
    // perform merge
    merging(
      ch_run
    )
  } else {
    // perform integration
    integrating(
      ch_run
    )
  }
  
}

