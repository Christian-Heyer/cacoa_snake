if(!require("cacoa")) {
    install.packages("coda.base", repos =  "https://cloud.r-project.org")
    #devtools::install_github("kharchenkolab/sccore", ref = "dev")
    devtools::install_github("kharchenkolab/sccore")
    devtools::install_github("kharchenkolab/cacoa")

}
library(cacoa)
library(igraph)
library(magrittr)
library(Seurat)
library(future)
library(org.Mm.eg.db)
if (exists("snakemake")) {
    cao_path <- snakemake@input[["cacoa_obj"]]
    cao_output <- snakemake@output[["cacoa_processed"]]
    xlsx_path <- snakemake@output[["cao_xlsx"]]
    threads <- snakemake@threads
    organism <- snakemake@config[["organism"]]
    is_raw <- snakemake@config[["is_raw"]]
       
    
} else {
    threads <- 2
    
    base_fp = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/TabularMuris"
    cao_path <- file.path(base_fp, "cao_obj.RDS.gz")
    cao_output <- file.path(base_fp, "processed_cao.RDS.gz")
    is_raw <- FALSE
}
#plan("multicore", workers = threads)
#options(error=function() traceback(2))
        
cao <- readRDS(cao_path)
        
cao$data.object<- FindNeighbors(cao$data.object,features = VariableFeatures(cao$data.object))
if(is_raw) {
    cao$data.object@misc$graph.name <- "SCT_snn"
} else {
    cao$data.object@misc$graph.name <- "RNA_snn"
}
cao$n.cores <- threads

cao$estimateCellLoadings()
cao$estimateCellDensity()
#cao$estimateCellDensity(method='kde', name='cell.density.kde',)

cao$estimateDiffCellDensity()
#cao$estimateDiffCellDensity(type='permutation', verbose=FALSE, name='cell.density.kde')
cao$estimateExpressionShiftMagnitudes(min.cells.per.sample = 5,top.n.genes = 1500, n.permutations = 2500)
#cao$estimateExpressionShiftMagnitudes()


cao$estimateClusterFreeDE()
exc.genes <- cao$test.results$cluster.free.z %>%  colnames() %>%
  .[grepl("Mt-",x = ., ignore.case = TRUE)] %>% c("Malat1")
#cao$smoothClusterFreeZScores(n.top.genes=1000, excluded.genes=exc.genes)

#cao$estimateClusterFreeExpressionShifts(n.top.genes=1500, n.permutations = 2500)
        
#cao$estimateGenePrograms(n.programs=10, n.top.genes = 1500)
## Estimate DE



cao$estimateDEPerCellType(independent.filtering=FALSE, name = "de.Wald", 
                          test='DESeq2.Wald', 
                          resampling.method=NULL, n.resamplings=10,
                         min.cell.count = 40)


cao$estimateDEPerCellType(n.cells.subsample= 50, name='deFixed_LRT',
                          n.resamplings=30)

cao$estimateDEPerCellType(name = 'de.fix',  
                            n.cells.subsample = 100,
                             min.cell.count = 100, n.resamplings = 30)
        
cao$estimateDEPerCellType(name='de.loo', resampling.method='loo', n.resamplings = 30)



estimateAllStabs <- function(cao_obj, de_n, org_db = org.Mm.eg.db::org.Mm.eg.db) {
    

    
    cao_obj$estimateDEStabilityPerGene(de.name = de_n,
                             top.n.genes = 500)
    cao_obj$estimateOntology(type="GSEA", org.db=org_db, verbose=T, n.cores=1,name = paste0(de_n, "GO"), de.name = de_n)

    
    return(cao_obj)
}

if (organism == "hs") {
    org_db  <- org.Hs.eg.db::org.Hs.eg.db
} else if (organism == "mm") {
    org_db <- org.Mm.eg.db::org.Mm.eg.db
}

for(de_n in list("de.Wald", "deFixed_LRT", "de.fix", "de.loo")) {
    cao <- estimateAllStabs(de_n = de_n, cao_obj = cao, org_db = org_db)
} 

#from cao_obj$test.results$de.loo extract the results from [[cell.type]]$res for each cell type

extract_res <- purrr::map(cao$test.results$de.loo, function(x) x$res)

writexl::write_xlsx(extract_res, path = xlsx_path)
saveRDS(cao,  cao_output ) 

