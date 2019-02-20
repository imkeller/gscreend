# pool Screen experiment class

setClass("PoolScreenExp", slots = c(
                   sgRNAData = "SummarizedExperiment",
                   FittingIntervals = "vector",
                   LFCModelParameters = "list",
                   GeneData = "SummarizedExperiment",
                   FittingOptions = "list")
         )
