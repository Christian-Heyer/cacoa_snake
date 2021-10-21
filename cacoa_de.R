if(!require("cacoa")) {
    install.packages("coda.base", repos =  "https://cloud.r-project.org")
    devtools::install_github("kharchenkolab/cacoa", ref = "dev", upgrade = "never")
}
library(cacoa)
library(igraph)
library(magrittr)
library(Seurat)
library(future)
if (exists("snakemake")) {
    cao_path <- snakemake@input[["cacoa_obj"]]
    cao_output <- snakemake@output[["cacoa_processed"]]
    threads <- snakemake@threads
} else {
    threads <- 2

    base_fp = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/adams_et_al//"
    cao_path <- file.path(base_fp, "cao_obj.RDS.gz")
    cao_output <- file.path(base_fp, "processed_cao.RDS.gz")
}
plan("multicore", workers = threads)
#options(error=function() traceback(2))
        
cao <- readRDS(cao_path)
        
cao$data.object<- FindNeighbors(cao$data.object,features = VariableFeatures(cao$data.object))
cao$data.object@misc$graph.name <- "SCT_snn"
cao$n.cores <- threads

cao$estimateCellLoadings()
cao$estimateCellDensity()
cao_obj$estimateCellDensity(method='graph', name='cell.density.graph',)

cao$estimateDiffCellDensity(type='permutation', verbose=FALSE)
cao_obj$estimateDiffCellDensity(type='permutation', verbose=FALSE, name='cell.density.graph')
cao$estimateExpressionShiftMagnitudes(min.cells.per.sample=10, dist="cor")

cao$estimateCommonExpressionShiftMagnitudes(n.cores = threads)

cao$estimateClusterFreeDE(min.expr.frac=0.01)
exc.genes <- cao_obj$test.results$cluster.free.z %>%  colnames() %>%
  .[grepl("Mt-", .)] %>% c("Malat1")
cao_obj$smoothClusterFreeZScores(n.top.genes=1000, progress.chunks=10, excluded.genes=exc.genes)
cao$estimateGenePrograms(n.programs=10)

cao$estimateClusterFreeExpressionShifts(n.top.genes=3000)
        
## Estimate DE

cao$estimateDEPerCellType(independent.filtering=TRUE, name = "de.Wald", 
                          test='DESeq2.Wald', 
                          resampling.method='bootstrap', n.resamplings=30)

cao$estimateDEPerCellType(n.cells.subsample= 50, name='deFixed_LRT',
                          resampling.method='bootstrap', n.resamplings=30)

cao$estimateDEPerCellType(resampling.method = 'fix.cells', name = 'de.fix',  
                            n.cells.subsample = 100,
                             min.cell.count = 100)
        
cao$estimateDEPerCellType(name='de.loo', resampling.method='loo', n.resamplings = 30)



estimateAllStabs <- function(cao_obj, de_n) {
    cao_obj$estimateDEStabilityPerCellType(top.n.genes = 300, de.name = de_n, 
                                     name = paste0(de_n,'_stab.fix'))

    cao_obj$estimateDEStabilityTrend(de.name = de_n, name = paste0(de_n, "_trend"),
                             top.n.genes = seq(50,500,50))  
    
    return(cao_obj)
}

for(de_n in list("de.Wald", "deFixed_LRT", "de.fix", "de.loo")) {
    cao <- estimateAllStabs(de_n = de_n, cao_obj = cao)
} 

saveRDS(cao,  cao_output ) 

