install.packages("BiocManager")
BiocManager::install("topGO")
BiocManager::install("org.Hs.eg.db")  # Assuming you are working with human data

# Load the libraries
library(topGO)
library(org.Hs.eg.db)

# Load the CSV file
significant_genes <- read.delim("all_genes.csv")

head(significant_genes)  # check out first few rows

# create the gene vector of p-values with names attached
geneList <- significant_genes$pvalue
names(geneList) <- significant_genes$Gene

# Convert p-values into a binary factor to show significance
geneFactor <- factor(as.integer(geneList < 0.05))
levels(geneFactor)
names(geneFactor) <- names(geneList)

# Create the topGO data object
GOdata <- new("topGOdata",
              ontology = "BP", 
              allGenes = geneFactor,
              annot = annFUN.org,  
              mapping = "org.Hs.eg.db", 
              ID = "symbol")  # Use gene symbols as IDs

# Run Fisher's exact test for the actual enrichment analysis
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Get the top significant GO terms
topGOres <- GenTable(GOdata, 
                     classicFisher = resultFisher, 
                     orderBy = "classicFisher", 
                     topNodes = 10)  # keeping top 10

# Print the results just to check
print(topGOres)

# Save the results to a CSV file
write.csv(topGOres, "GSEA_topGO_results.csv")

