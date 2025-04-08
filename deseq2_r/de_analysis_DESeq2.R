# prepare environment
rm(list=ls()) # clear envoronment

# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install required Bioconductor packages if not already installed
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

# Install required Bioconductor packages if not already installed
if (!requireNamespace("fgsea", quietly = TRUE))
  BiocManager::install("fgsea")

if (!requireNamespace("msigdbr", quietly = TRUE))
  BiocManager::install("msigdbr")

# Install CRAN packages if not already installed
if (!requireNamespace("dplyr", quietly = TRUE))
  install.packages("dplyr")

if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

if (!requireNamespace("tibble", quietly = TRUE))
  install.packages("tibble")


# Load required libraries
library(DESeq2)
library(ggplot2)
library(fgsea)
library(msigdbr)
library(dplyr)
library(tibble)  # for deframe()

# Set up your working directory that contains all the working files. Verify that you have set it correctly. See
# whether you see all the working files (.map .ped .list)
setwd("C:\\Monikos\\gDrive\\Osmianskiene\\VU Sistemų biologija\\Transkriptomika\\Practice - Gediminas Alzbutas\\")
getwd()
list.files()

print("Hello, world!")

##############################################################
# 1. Differential Expression for Collibri
##############################################################

# Read in the Collibri count matrix
collibri_counts <- read.table("data_output_matrix/Collibri.txt", header = TRUE, row.names = 1, sep="\t")

# Create sample information by extracting condition from column names
# (Samples with “UHRR” are cancer-derived and those with “HBR” are normal)
collibri_samples <- colnames(collibri_counts)
collibri_condition <- ifelse(grepl("UHRR", collibri_samples), "UHRR", "HBR")
collibri_colData <- data.frame(row.names = collibri_samples,
                               condition = factor(collibri_condition))

# Build DESeqDataSet for Collibri and run the DESeq2 pipeline
dds_collibri <- DESeqDataSetFromMatrix(countData = collibri_counts,
                                       colData = collibri_colData,
                                       design = ~ condition)

# Full differential expression analysis on DESeqDataSet object
# Normalization, Dispersion Estimation, Model Fitting, Hypothesis Testing
dds_collibri <- DESeq(dds_collibri)

# Extract the differential expression results from the DESeqDataSet object
res_collibri <- results(dds_collibri)

# Order results by adjusted p-value and extract DE genes (padj < 0.05)
res_collibri_ordered <- res_collibri[order(res_collibri$padj), ]
de_genes_collibri <- subset(as.data.frame(res_collibri_ordered), padj < 0.05)

# A listof DE genes wit Padj values
df_c <- as.data.frame(res_collibri)
df_c$gene <- rownames(df_c)
write.csv(df_c[, c("gene", "padj")], "r output/DE_collibri_results.txt", row.names = FALSE)

write.csv(as.data.frame(res_collibri_ordered), "r output/DE_collibri_ordered_results.txt")
write.csv(as.data.frame(de_genes_collibri), "r output/DE_collibri_genes_results.txt")

# Define a function to generate a volcano plot with colored points
make_volcano <- function(res, title) {
  res_df <- as.data.frame(res)
  # Mark genes significant if padj < 0.05 and absolute log2FoldChange > 1
  res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Yes", "No")
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("green", "red")) +
    labs(title = title, x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
    theme_minimal()
}

# Print the volcano plot to the R console
print(make_volcano(res_collibri, "Volcano Plot - Collibri"))

# Save the volcano plot to a PNG file
png("r output/volcano_collibri.png", width = 800, height = 600)
print(make_volcano(res_collibri, "Volcano Plot - Collibri"))
dev.off()

##############################################################
# 2. Differential Expression for KAPA
##############################################################

# Read in the KAPA count matrix
kapa_counts <- read.table("data_output_matrix/KAPA.txt", header = TRUE, row.names = 1, sep="\t")

# Create sample information for KAPA by extracting the condition from sample names
kapa_samples <- colnames(kapa_counts)
kapa_condition <- ifelse(grepl("UHRR", kapa_samples), "UHRR", "HBR")
kapa_colData <- data.frame(row.names = kapa_samples,
                           condition = factor(kapa_condition))

# Build DESeqDataSet for KAPA and run the analysis
dds_kapa <- DESeqDataSetFromMatrix(countData = kapa_counts,
                                   colData = kapa_colData,
                                   design = ~ condition)
# Full differential expression analysis on DESeqDataSet object
# Normalization, Dispersion Estimation, Model Fitting, Hypothesis Testing
dds_kapa <- DESeq(dds_kapa)

# Extract the differential expression results from the DESeqDataSet object
res_kapa <- results(dds_kapa)

# Order results and extract DE genes (padj < 0.05)
res_kapa_ordered <- res_kapa[order(res_kapa$padj), ]
de_genes_kapa <- subset(as.data.frame(res_kapa_ordered), padj < 0.05)

# A listof DE genes with Padj values
df_k <- as.data.frame(res_kapa)
df_k$gene <- rownames(df_k)
write.csv(df_k[, c("gene", "padj")], "r output/DE_KAPA_results.txt", row.names = FALSE)

write.csv(as.data.frame(res_kapa_ordered), "r output/DE_KAPA_ordered_results.txt")
write.csv(as.data.frame(de_genes_kapa), "r output/DE_KAPA_genes_results.txt")

# Define a function to generate a volcano plot with colored points
make_volcano_KAPA <- function(res, title) {
  res_df <- as.data.frame(res)
  # Mark genes significant if padj < 0.05 and absolute log2FoldChange > 1
  res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, "Yes", "No")
  ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.7) +
    scale_color_manual(values = c("blue", "red")) +
    labs(title = title, x = "log2 Fold Change", y = "-log10 Adjusted P-value") +
    theme_minimal()
}

# Print the volcano plot to the R console
print(make_volcano_KAPA(res_kapa, "Volcano Plot - KAPA"))

# Save the volcano plot to a PNG file
png("r output/volcano_KAPA.png", width = 800, height = 600)
print(make_volcano_KAPA(res_kapa, "Volcano Plot - KAPA"))
dev.off()

##############################################################
# 3. PCA Plot Based on Common DE Genes
##############################################################

# Identify genes that are DE (padj < 0.05) in both methods
common_de_genes <- intersect(rownames(de_genes_collibri), rownames(de_genes_kapa))
cat("Number of common DE genes:", length(common_de_genes), "\n")

# Transform counts (rlog transformation) for each dataset
rld_collibri <- rlog(dds_collibri, blind = TRUE)
rld_kapa <- rlog(dds_kapa, blind = TRUE)

# Subset the rlog-transformed data for common DE genes
rld_collibri_common <- assay(rld_collibri)[common_de_genes, ]
rld_kapa_common <- assay(rld_kapa)[common_de_genes, ]

# Combine the transformed data from both sample preparations
combined_counts <- cbind(rld_collibri_common, rld_kapa_common)

# Create a metadata table for the combined samples.
combined_samples <- colnames(combined_counts)
combined_method <- c(rep("Collibri", ncol(rld_collibri_common)),
                     rep("KAPA", ncol(rld_kapa_common)))
combined_condition <- c(as.character(collibri_colData$condition),
                        as.character(kapa_colData$condition))
combined_colData <- data.frame(row.names = combined_samples,
                               method = factor(combined_method),
                               condition = factor(combined_condition))

# Perform PCA using the combined rlog data (transpose so samples are rows)
pca <- prcomp(t(combined_counts))
pca_data <- data.frame(PC1 = pca$x[,1],
                       PC2 = pca$x[,2],
                       method = combined_colData$method,
                       condition = combined_colData$condition)

# Create PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = method, shape = condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle("PCA Plot Based on Common DE Genes") +
  xlab("PC1") +
  ylab("PC2")
print(pca_plot)

##############################################################
# 4. Saving DE Gene Lists (Optional)
##############################################################

write.table(de_genes_collibri, file = "DE_genes_Collibri.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)
write.table(de_genes_kapa, file = "DE_genes_KAPA.txt", sep = "\t",
            quote = FALSE, row.names = TRUE)

##############################################################
# 5. Which Fold Changes to Use?
##############################################################
# For enrichment analysis using fgsea, you should use the unshrunken fold changes.
# Although DESeq2 provides shrunken fold changes via lfcShrink for interpretation,
# these shrink the magnitude of changes and can affect gene ranking.
# Therefore, I used the original (unshrunken) log2 fold changes for the ranking.

##############################################################
# 6. Preparing Gene Rankings for fgsea
##############################################################

# Prepare ranking for Collibri
# Note: Our gene IDs contain version numbers (e.g., ENSG00000011304.12) so we remove them.
rank_collibri <- res_collibri$log2FoldChange
names(rank_collibri) <- sapply(strsplit(rownames(res_collibri), "\\."), `[`, 1)
rank_collibri <- sort(rank_collibri, decreasing = TRUE)

# Prepare ranking for KAPA
rank_kapa <- res_kapa$log2FoldChange
names(rank_kapa) <- sapply(strsplit(rownames(res_kapa), "\\."), `[`, 1)
rank_kapa <- sort(rank_kapa, decreasing = TRUE)

##############################################################
# 7. Obtaining Reactome Pathways
##############################################################
# We use msigdbr to get Reactome gene sets.
# Here we extract Ensembl gene IDs so they match our ranking names.
reactome_pathways <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(ensembl_gene, gs_name)

# Build a list of pathways (each element is a vector of Ensembl gene IDs)
reactome_list <- reactome_pathways %>% 
  group_by(gs_name) %>% 
  summarise(genes = list(unique(ensembl_gene))) %>% 
  deframe()

##############################################################
# 8. Running fgsea for Both Kits
##############################################################
set.seed(42)  # for reproducibility
fgsea_collibri <- fgsea(pathways = reactome_list,
                        stats = rank_collibri,
                        nperm = 1000)
fgsea_kapa <- fgsea(pathways = reactome_list,
                    stats = rank_kapa,
                    nperm = 1000)

# Order results by adjusted p-value
fgsea_collibri <- fgsea_collibri[order(fgsea_collibri$padj), ]
fgsea_kapa <- fgsea_kapa[order(fgsea_kapa$padj), ]
fgsea_collibri
fgsea_kapa

##############################################################
# 9. Examining Enriched Pathways
##############################################################
# Print top 10 enriched pathways for each kit
print("Top enriched pathways in Collibri:")
print(head(fgsea_collibri, 10))

print("Top enriched pathways in KAPA:")
print(head(fgsea_kapa, 10))

# Identify cancer-related pathways by searching for the term "cancer" (case-insensitive)
cancer_pathways_collibri <- fgsea_collibri[grep("cancer", fgsea_collibri$pathway, ignore.case = TRUE), ]
cancer_pathways_kapa <- fgsea_kapa[grep("cancer", fgsea_kapa$pathway, ignore.case = TRUE), ]

print("Cancer-related pathways in Collibri:")
print(cancer_pathways_collibri)

print("Cancer-related pathways in KAPA:")
print(cancer_pathways_kapa)

##############################################################
# 10. Comparing Results Between Kits
##############################################################
# Merge the fgsea results from both kits for common pathways
merged_fgsea <- merge(fgsea_collibri, fgsea_kapa, by = "pathway", suffixes = c("_collibri", "_kapa"))

# Plot comparison of normalized enrichment scores (NES) between kits
comparison_plot <- ggplot(merged_fgsea, aes(x = NES_collibri, y = NES_kapa, label = pathway)) +
  geom_point() +
  geom_text(size = 3, vjust = 1.5, hjust = 1.5) +
  theme_minimal() +
  ggtitle("Comparison of NES: Collibri vs KAPA") +
  xlab("NES Collibri") +
  ylab("NES KAPA")
print(comparison_plot)

##############################################################
# 11. Summary Comments
##############################################################
# - I used the unshrunken log2 fold changes to rank genes because these better capture the differences for enrichment analysis.
# - The fgsea analysis on Reactome pathways produced a set of enriched pathways.
# - The printed output will show if any cancer-related pathways are enriched.
# - The comparison plot can help determine if there are systematic differences in enrichment results
#   between the two sample preparation kits.



