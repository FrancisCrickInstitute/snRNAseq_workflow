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
include { clustering  } from './modules/clustering'
include { annotating  } from './modules/annotating'
include { infercnv    } from './modules/infercnv'


workflow {

  // save params
  save_params()

  // integrate all preloaded samples
  Channel
    .fromPath("${params.output.dir}/by_sample/*/filtering/seu.rds")
    .map { file -> tuple("all", file) }
    .groupTuple()
    .map { dir, file -> tuple(dir, "", file)}
    .set { ch_run }

  // perform integration
  integrating(
    "${baseDir}/templates/integrating.rmd",
    ch_run,
    save_params.out.ch_params
  )

}
