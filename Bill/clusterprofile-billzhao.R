#Loading Libraries

library(clusterProfiler)

library(org.Hs.eg.db)

# Reading dataframe of most significant genes and their associated log2FoldChanges
gene_summary = read.csv("C:/Users/aquat/OneDrive/Documents/GitHub/BioInfoAlg_RheumArth/most_significant.csv")
gene_list <- gene_summary$log2FoldChange
gene_symbols <- gene_summary$Gene

# Conversion to EntrezID
gene_entrez <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Attaching log2FoldChange scores to each ID
names(gene_list) <- gene_entrez$ENTREZID

# Sorting for appropriate use in gseGO
gene_list <- sort(gene_list, decreasing=TRUE)
gene_list <- gene_list[!duplicated(names(gene_list))]
length(gene_list)

# Running Gene Enrichment
gsea_go <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, ont = "BP",
                 + minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05, verbose=FALSE)

# Writing File Out
write.csv(gsea_go, "C:/Users/aquat/OneDrive/Documents/GitHub/BioInfoAlg_RheumArth/enrichment.csv")