# sgRNA slot
setGeneric("sgRNAData", function(x) standardGeneric("sgRNAData"))
setGeneric("sgRNAData<-", function(x, value) standardGeneric("sgRNAData<-"))

# Gene slot
setGeneric("GeneData", function(x) standardGeneric("GeneData"))
setGeneric("GeneData<-", function(x, value) standardGeneric("GeneData<-"))

# set and get normalized counts
setGeneric("normcounts", function(x) standardGeneric("normcounts"))
setGeneric("normcounts<-", function(x, value) standardGeneric("normcounts<-"))

# get reference(library) counts
setGeneric("refcounts", function(x) standardGeneric("refcounts"))

# get sample counts
setGeneric("samplecounts", function(x) standardGeneric("samplecounts"))

# get sample lfc
setGeneric("samplelfc", function(x) standardGeneric("samplelfc"))

# get sample pvalues
setGeneric("samplepval", function(x) standardGeneric("samplepval"))
