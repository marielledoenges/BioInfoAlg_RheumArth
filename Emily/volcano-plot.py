import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from gseapy import enrichr, barplot

# Load the differential expression results from the CSV file
de_results = pd.read_csv('de_results_optimized.csv')

# Create the Volcano Plot
print(f"Making Volcano plot...")
plt.figure(figsize=(10, 6))

# Scatter plot with coloring based on thresholds for significance
plt.scatter(de_results['Log2FoldChange'], -np.log10(de_results['PValue']),
            c=np.where((de_results['PValue'] < 0.05) & (abs(de_results['Log2FoldChange']) > 1), 'red', 'gray'),
            alpha=0.75)

plt.title('Volcano Plot of Differential Gene Expression')
plt.xlabel('Log2 Fold Change')
plt.ylabel('-log10 P-value')

# Add threshold lines
plt.axhline(y=-np.log10(0.05), color='black', linestyle='--')  # P-value threshold
plt.axvline(x=1, color='black', linestyle='--')  # Log2 fold change threshold (upregulated)
plt.axvline(x=-1, color='black', linestyle='--')  # Log2 fold change threshold (downregulated)

# Show the plot
plt.tight_layout()
plt.show()

#significant_genes = de_results[(de_results['PValue'] < 0.05) & (abs(de_results['Log2FoldChange']) > 0.5)]

significant_gene_names = significant_genes['Gene'].tolist()
print("Systematically significant genes:")
print(significant_gene_names)

print("Top 50 significant genes:")
#top_50_genes = significant_genes.sort_values(by=['PValue', 'Log2FoldChange'], descending=[False, True]).head(50)
#print(top_50_genes)

# GENE EXPRESSION (Emily) INDV PORTION
#significant_gene_list = top_50_genes['Gene'].tolist()

significant_genes.to_csv('significant_genes.csv', index=False)

# Perform Gene Set Enrichment Analysis (GSEA) using gseapy
#print("Performing GSEA...")

#enrichr_results = enrichr(gene_list=top_50_gene_list,
#                          gene_sets='GO_Biological_Process_2021',
#                          organism='Human') 


#enrichr_results.res2d.to_csv('gsea_results.csv', index=False)
#print("GSEA results saved to 'gsea_results.csv'.")


#print("Generating barplot of top enriched terms...")
#barplot(enrichr_results.res2d, title='Top Enriched Terms (GO Biological Process)', cutoff=0.05)

#print("GSEA completed.")


