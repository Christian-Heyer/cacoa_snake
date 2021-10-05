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
    threads <- 8

    base_fp = "/omics/odcf/analysis/OE0228_projects/VascularAging/rna_sequencing/public_scrnaseq/TabularMuris/"
    cao_path <- file.path(base_fp, "cao_obj.RDS.gz")
    cao_output <- file.path(base_fp, "processed_cao.RDS.gz")
}
plan("multicore", workers = threads)
#options(error=function() traceback(2))
        
cao <- readRDS(cao_path)
        
cao$data.object<- FindNeighbors(cao$data.object,features = VariableFeatures(cao$data.object))
cao$data.object@misc$graph.name <- "SCT_snn"
cao$n.cores <- 8

cao$estimateCellLoadings()
cao$estimateCellDensity()

cao$estimateDiffCellDensity(type='permutation', verbose=FALSE)
cao$estimateExpressionShiftMagnitudes(min.cells=10, n.cells=1e3, dist="cor", n.subsamples=50)

cao$estimateCommonExpressionShiftMagnitudes(n.cores = 8)

cao$estimateClusterFreeDE(min.expr.frac=0.01)

cao$smoothClusterFreeZScores()
cao$estimateGenePrograms(n.programmes=10, n.sampled.cells=10000, cyc=1500)

cao$estimateClusterFreeExpressionShifts(n.top.genes=3000)
        
## Estimate DE

cao$estimateDEPerCellType(independent.filtering=TRUE, name = "de.Wald", 
                          test='DESeq2.Wald', 
                          resampling.method='bootstrap', max.resamplings=30)

cao$estimateDEPerCellType(max.cell.count = 50, name='deFixed_LRT',
                          resampling.method='bootstrap', max.resamplings=30)

cao$estimateDEPerCellType(resampling.method = 'fix.count', name = 'de.fix', max.cell.count = 100, min.cell.count = 100)
        
cao$estimateDEPerCellType(name='de.loo', resampling.method='loo')



estimateAllStabs <- function(cao_obj, de_n) {
    cao_obj$estimateDEStabilityPerCellType(top.n.genes = 300, de.name = de_n, 
                                     name = paste0(de_n,'_stab.fix'))

    cao_obj$estimateDEStabilityTrend(de.name = de_n, name = paste0(de_n, "_trend"),
                             top.n.genes = seq(50,500,50))  
    
    return(cao_obj)
}

for(de_n in list("de.Wald", "deFixed_LRT", "de.fix", "de.loo")) {
    cao_obj <- estimateAllStabs(de_n = de_n, cao_obj = cao_obj)
} 

saveRDS(cao,  cao_output ) 

