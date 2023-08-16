#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// differential expression analysis
process pathway_analysis {
  tag "${patient}"
  label 'process_medium'
  publishDir "${params.output.dir}/by_patient_wo_organoids/${patient}/integrating/differential_expression/",
    mode: 'copy',
    pattern: "{*.rds}"
    
  input:
    tuple val(patient), path(rds_file)
    
  output:
    path 'scpa.rds'
    
  script:
    """
    #!/usr/bin/env Rscript
    library(magrittr)
    
    print("load seu")
    seu <- readRDS("seu.rds")
    
    print("convert to cds")
    cds <- SeuratWrappers::as.cell_data_set(seu, assay = "RNA")
    SummarizedExperiment::rowData(cds)[, "gene_short_name"] <- 
      row.names(SummarizedExperiment::rowData(cds))
    
    print("preprocess")
    cds <- monocle3::preprocess_cds(cds, num_dim = 30)
    cds <- monocle3::reduce_dimension(cds)

    print("get pathways")
    pathways <-
      msigdbr::msigdbr(species = "Homo sapiens", category = "H") %>%
      SCPA::format_pathways()
      
    print("SCPA pathway analysis")
    scpa <-
      SCPA::compare_sce(cds,
                        group1 = "label.infercnv",
                        group1_population = c("Malignant_cells", "Epithelial_cells"),
                        pathways = pathways)
    
    print("save rds")
    saveRDS(scpa, "scpa.rds")
    """
}

workflow {
  
  // get rds files
  Channel
    .fromPath("${params.output.dir}/by_patient_wo_organoids/*/integrating/infercnv/seu.rds")
    .map { rds_file -> tuple(rds_file.toString().tokenize('/')[-4], rds_file) }
    .set { ch_run }
  
  // differential expression
  pathway_analysis(
    ch_run
  )
  
}

