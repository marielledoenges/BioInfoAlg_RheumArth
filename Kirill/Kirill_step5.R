# Install and load required libraries if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install clusterProfiler for Gene Ontology Enrichment
if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

# Install org.Hs.eg.db for human gene annotations
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Hs.eg.db")
}

# Load libraries
library(clusterProfiler)
library(org.Hs.eg.db)

# Step 1: Load your gene expression data
new_data <- read.csv("all_genes.csv", header = TRUE)

# Step 2: Filter genes with adjusted p-value < 0.05
significant_genes <- subset(new_data, padj < 0.05)

# Extract the list of significant gene names (assuming 'Gene' column has gene symbols)
gene_list <- significant_genes$Gene

# Step 3: Convert gene symbols to Entrez IDs (required for clusterProfiler)
# This uses the org.Hs.eg.db package to map gene symbols to Entrez IDs
entrez_ids <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Step 4: Run GO Enrichment Analysis using clusterProfiler, focusing on Molecular Function (MF)
go_enrichment <- enrichGO(gene = entrez_ids$ENTREZID, 
                          OrgDb = org.Hs.eg.db, 
                          keyType = "ENTREZID", 
                          ont = "CC",       # Use "MF" for Molecular Function ontology
                          pAdjustMethod = "BH", 
                          pvalueCutoff = 0.2, 
                          readable = TRUE)    # Convert Entrez IDs back to gene symbols

# Step 5: View the top enriched GO terms
head(go_enrichment)

# Step 6: Save the GO enrichment results to a CSV file
write.csv(as.data.frame(go_enrichment), "go_enrichment_molecular_function.csv", row.names = FALSE)

# Step 7: Visualize the enrichment results (optional)
# Bar plot of the top enriched GO terms
barplot(go_enrichment, showCategory = 10, title = "Top 10 Enriched GO Terms (Molecular Function)")

# Dot plot of the GO terms
dotplot(go_enrichment, showCategory = 10, title = "Top 10 Enriched GO Terms (Molecular Function)")
