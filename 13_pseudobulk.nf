#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// differential expression analysis
process pseudobulk {
  tag "${patient}"
  label 'process_medium'
  publishDir "${params.output.dir}/by_patient_wo_organoids/${patient}/pseudobulk/",
    mode: 'copy',
    pattern: "{*.rds}"
    
  input:
    tuple val(patient), path(rds_file)
    
  output:
    path 'pse.rds', optional: true
    
  script:
    """
    #!/usr/bin/env Rscript
    # load seu, convert to sce
    seu <- readRDS("seu.rds")
    Seurat::DefaultAssay(seu) <- "RNA"
    sce <- Seurat::as.SingleCellExperiment(seu)
    
    # take highly expressed genes and proper cells
    sce <- sce[
      rowSums(as.matrix(BiocGenerics::counts(sce))) > 5,
      sce\$label.infercnv != "Others"]
    sce\$cycling_score <- rowMeans(cbind(sce\$S.Score, sce\$G2M.Score))
    sce\$response <- factor(
      sce\$response, levels = c("Progressor", "Non-responder", "Partial Responder", "Responder")
    )
    
    # convert to dense matrix
    BiocGenerics::counts(sce) <- as.matrix(BiocGenerics::counts(sce))
    
    # pseudobulk
    pse <- glmGamPoi::pseudobulk(sce, group_by = glmGamPoi::vars(
      patient_id, label.infercnv, timepoint, sample_site),
      malignancy = mean(malignancy), age = unique(age), 
      treatment = unique(treatment), 
      response = unique(response),
      cycling_score = mean(cycling_score),
      gender = unique(gender))
    
    # save
    saveRDS(pse, "pse.rds")
    """
}

workflow {
  
  // get seu_annotated files
  Channel
    .fromPath("${params.output.dir}/by_patient_wo_organoids/*/integrating/infercnv/seu.rds")
    .map { rds_file -> tuple(rds_file.toString().tokenize('/')[-4], rds_file) }
    .set { ch_run }
  
  // differential expression
  pseudobulk(
    ch_run
  )
  
}

