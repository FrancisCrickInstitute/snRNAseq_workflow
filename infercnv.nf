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
include { infercnv    } from './modules/infercnv'


workflow {

  // save params
  save_params()

  // integrate all preloaded samples
  Channel
    .fromPath("${params.output.dir}/by_patient_wo_organoids/*/integrating/seurat_clustering/seu_annotated.rds")
    .set { ch_run }

  // get id, dir, subdir and filepath
  ch_run
    .map { file ->
            tuple(file.toString().tokenize('/')[-4],
                  file.toString().tokenize('/')[-5],
                  file.toString().tokenize('/')[-3],
                  file) }
    .set { ch_run }

  // perform integration
  infercnv(
    "${baseDir}/templates/infercnv.rmd",
    ch_run,
    save_params.out.ch_params
  )

}
