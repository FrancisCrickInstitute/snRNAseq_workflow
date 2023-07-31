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
include { integrating } from './modules/integrating'


workflow {

  // save params
  save_params()

  // initialise
  Channel
    .fromPath("${params.output.dir}/by_sample/*/filtering/seu.rds")
    .map { rds_file -> tuple(rds_file.toString().tokenize('/')[-3], rds_file) }
    .set { ch_files }
  Channel
    .empty()
    .set { ch_run }

  if ( params.input.run_by_patient ) {
    // by patient
    ch_files
      .map { id, rds_file -> tuple(id.split("_")[0], rds_file) }
      .groupTuple()
      .map { id, rds_file -> tuple(id, "by_patient", rds_file) }
      .concat(ch_run)
      .set { ch_run }
  }

  if ( params.input.run_by_patient_wo_organoids ) {
    ch_files
      .filter { !it[0].contains("organoid") }
      .map { id, rds_file -> tuple(id.split("_")[0], rds_file) }
      .groupTuple()
      .map { id, rds_file -> tuple(id, "by_patient_wo_organoids", rds_file) }
      .concat(ch_run)
      .set { ch_run }
  }

  if ( params.input.run_all ) {
    ch_files
      .map { rds_file -> tuple("all", rds_file) }
      .groupTuple()
      .map { dir, rds_file -> tuple(dir, "", rds_file)}
      .set { ch_run }
  }

  // perform integration
  integrating(
    "${baseDir}/templates/integrating.rmd",
    ch_run,
    save_params.out.ch_params
  )

}
