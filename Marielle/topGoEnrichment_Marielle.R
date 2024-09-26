install.packages("BiocManager")
BiocManager::install("topGO")
BiocManager::install("org.Hs.eg.db")  # Assuming you are working with human data

# Load the libraries
library(topGO)
library(org.Hs.eg.db)

# Load the CSV file
significant_genes <- read.delim("all_genes.csv")

# Assuming your CSV has two columns: 'Gene' and 'p_value'
head(significant_genes)  # Inspect the first few rows of your data

# Prepare the gene list: a named vector of p-values
geneList <- significant_genes$pvalue
names(geneList) <- significant_genes$Gene

# Convert p-values into a binary factor indicating significance (e.g., threshold at 0.05)
geneFactor <- factor(as.integer(geneList < 0.05))
levels(geneFactor)
names(geneFactor) <- names(geneList)

# Create the topGOdata object
GOdata <- new("topGOdata",
              ontology = "BP",  # You can choose "MF" or "CC" based on your focus
              allGenes = geneFactor,
              annot = annFUN.org,  # Annotation function using the org package
              mapping = "org.Hs.eg.db",  # Adjust if using a different species
              ID = "symbol")  # Use gene symbols as IDs (adjust as needed)
# Run Fisher's exact test for enrichment analysis
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Get the top significant GO terms
topGOres <- GenTable(GOdata, 
                     classicFisher = resultFisher, 
                     orderBy = "classicFisher", 
                     topNodes = 10)  # Adjust topNodes for more/less results

# Print the results
print(topGOres)

# Save the results to a CSV file
write.csv(topGOres, "GSEA_topGO_results.csv")

# Optional: Visualize the most significant GO terms
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')


