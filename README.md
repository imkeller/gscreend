# gscreend - analysis of pooled CRIPR screens


### Run gscreend, as also explained in the vignette section

Set up a SummarizedExperiment object containing a matrix of raw count data, rowData on gRNAs and genes and colData on the sample type.

```{r}
counts_matrix <- cbind(raw_counts$library0, raw_counts$R0_0, raw_counts$R1_0)

rowData <- data.frame(sgRNA_id = raw_counts$sgrna_id,
                           gene = raw_counts$Gene)

colData <- data.frame(samplename = c("library", "R1", "R2"),
                      # timepoint naming convention: 
                      # T0 -> reference, 
                      # T1 -> selected
                      timepoint = c("T0", "T1", "T1"))

se <- SummarizedExperiment(assays=list(counts=counts_matrix),
                        rowData=rowData, colData=colData)
```

Create a PoolScreenExp object

```{r}
pse <- createPoolScreenExp(se)
```

Run gscreend

```{r}
pse <- RunGscreend(pse)

```
