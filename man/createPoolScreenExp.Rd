% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DataInput.R
\name{createPoolScreenExp}
\alias{createPoolScreenExp}
\title{Create PoolScreenExp Experiment}
\usage{
createPoolScreenExp(data)
}
\arguments{
\item{data}{Input data object containing gRNA level data
(SummarizedExperiment)}
}
\value{
object  PoolScreenExp object
}
\description{
Create PoolScreenExp Experiment
}
\examples{
raw_counts <- read.table(
                        system.file('extdata', 'simulated_counts.txt',
                        package = 'gscreend'),
                        header=TRUE)

counts_matrix <- cbind(raw_counts$library0, raw_counts$R0_0, raw_counts$R1_0)

rowData <- data.frame(sgRNA_id = raw_counts$sgrna_id,
gene = raw_counts$Gene)

colData <- data.frame(samplename = c('library', 'R1', 'R2'),
timepoint = c('T0', 'T1', 'T1'))

library(SummarizedExperiment)
se <- SummarizedExperiment(assays=list(counts=counts_matrix),
rowData=rowData, colData=colData)

# create a PoolScreenExp experiment
pse <- createPoolScreenExp(se)

}
