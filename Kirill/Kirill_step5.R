# Install and load required libraries if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Install GenomicSuperSignature if needed
if (!requireNamespace("GenomicSuperSignature", quietly = TRUE)) {
  BiocManager::install("GenomicSuperSignature")
}
# Install bcellViper package if needed (required for demo data)
if (!requireNamespace("bcellViper", quietly = TRUE)) {
  BiocManager::install("bcellViper")
}

# Load libraries
library(GenomicSuperSignature)
library(bcellViper)

# Step 1: Download and load the pre-trained PLIERpriors model
RAVmodel <- getModel("PLIERpriors", load = TRUE)
print(RAVmodel)

# Step 2: Load your gene expression data
new_data <- read.csv("all_genes.csv", header = TRUE)
head(new_data)

# Filter genes with adjusted p-value < 0.05
significant_genes <- subset(new_data, padj < 0.05)
nrow(significant_genes)
head(significant_genes)

# Extract the list of significant gene names
gene_list <- significant_genes$Gene
head(gene_list)

# Step 3: Perform validation on your dataset using the pre-trained model
# Use 'validate()' function to validate your data against the RAV model
# Here, for demonstration, the dset is used, replace it with your data as needed.
data(bcellViper)
val_all <- validate(dset, RAVmodel)
head(val_all)

# Step 4: Visualize the validation results
heatmapTable(val_all, RAVmodel, num.out = 5, swCutoff = 0)

# Step 5: Plot the validation results
plotValidate(val_all, interactive = FALSE)

# Step 6: Explore validated signatures and pathways
validated_ind <- validatedSignatures(val_all, RAVmodel, num.out = 3, swCutoff = 0, indexOnly = TRUE)
validated_ind

# Optional: Draw word clouds of the top signatures
drawWordcloud(RAVmodel, validated_ind[1])
drawWordcloud(RAVmodel, validated_ind[2])

# Step 7: Annotate the top enriched pathways
annotateRAV(RAVmodel, validated_ind[2])

# Step 8: Subset and display enriched pathways
subsetEnrichedPathways(RAVmodel, validated_ind[2], include_nes = TRUE)