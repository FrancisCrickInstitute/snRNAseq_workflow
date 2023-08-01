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

}

