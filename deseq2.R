print("Hello, world!")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager", lib="~/R/library")

BiocManager::install("DESeq2", lib="~/R/library")
BiocManager::install("GenomicRanges", lib="~/R/library")

library("DESeq2", lib="~/R/library")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition)
dds