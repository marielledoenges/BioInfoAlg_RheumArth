import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from umap import UMAP
from scipy.stats import ttest_ind


df = pd.read_csv("SRP158491_converted.tsv", sep='\t')

patient_df = pd.read_csv("SRP158491/metadata_SRP158491.tsv", sep='\t')

healthy_rows = patient_df[patient_df['refinebio_disease'] == 'healthy']
unhealthy_rows = patient_df[patient_df['refinebio_disease'] == 'ra non treatment']

healthy_columns = healthy_rows['refinebio_accession_code'].tolist()
unhealthy_columns = unhealthy_rows['refinebio_accession_code'].tolist()

healthy_gene_expression = df[['Gene'] + healthy_columns]
unhealthy_gene_expression = df[['Gene'] + unhealthy_columns]

combined_gene_expression = pd.concat([healthy_gene_expression.iloc[:, 1:], unhealthy_gene_expression.iloc[:, 1:]], axis=1)

# Prepare data for t-SNE and UMAP
combined_gene_expression_transpose = combined_gene_expression.transpose()

# Apply t-SNE
tsne = TSNE(n_components=2, random_state=42)
tsne_results = tsne.fit_transform(combined_gene_expression_transpose)

# Apply UMAP
umap = UMAP(n_components=2, random_state=42)
umap_results = umap.fit_transform(combined_gene_expression_transpose)

labels = ['Healthy'] * len(healthy_columns) + ['Unhealthy'] * len(unhealthy_columns)

# Plot t-SNE
plt.figure(figsize=(14, 6))

plt.subplot(1, 2, 1)
plt.scatter(tsne_results[:len(healthy_columns), 0], tsne_results[:len(healthy_columns), 1], c='blue', label='Healthy')
plt.scatter(tsne_results[len(healthy_columns):, 0], tsne_results[len(healthy_columns):, 1], c='red', label='Unhealthy')
plt.title('t-SNE of Healthy and Unhealthy Samples')
plt.xlabel('Component 1')
plt.ylabel('Component 2')
plt.legend()

# Plot UMAP
plt.subplot(1, 2, 2)
plt.scatter(umap_results[:len(healthy_columns), 0], umap_results[:len(healthy_columns), 1], c='blue', label='Healthy')
plt.scatter(umap_results[len(healthy_columns):, 0], umap_results[len(healthy_columns):, 1], c='red', label='Unhealthy')
plt.title('UMAP of Healthy and Unhealthy Samples')
plt.xlabel('Component 1')
plt.ylabel('Component 2')
plt.legend()

plt.tight_layout()
plt.show()

nan_genes_count = df['Gene'].isna().sum()

total_genes_count = len(df['Gene'])


# --- Differential Expression Analysis ---
# Print the healthy gene expression data
print("Healthy Gene Expression Data:")
print(healthy_gene_expression)

# Print the unhealthy gene expression data
print("Unhealthy Gene Expression Data:")
print(unhealthy_gene_expression)

results = []
for gene in df['Gene']:
    #print(f"Processing gene: {gene}")
    
    # Check if gene exists in both datasets (healthy and unhealthy)
    if gene in healthy_gene_expression['Gene'].values and gene in unhealthy_gene_expression['Gene'].values:
        print(f"Gene {gene} found in both datasets")
        healthy_expr = healthy_gene_expression.loc[healthy_gene_expression['Gene'] == gene, healthy_columns].values.flatten()
        unhealthy_expr = unhealthy_gene_expression.loc[unhealthy_gene_expression['Gene'] == gene, unhealthy_columns].values.flatten()

        # Perform t-test between healthy and unhealthy samples
        t_stat, p_val = ttest_ind(healthy_expr, unhealthy_expr, equal_var=False)

        # Log2 fold change (mean of unhealthy samples - mean of healthy samples)
        log_fold_change = np.log2(unhealthy_expr.mean() + 1) - np.log2(healthy_expr.mean() + 1)

        results.append([gene, log_fold_change, p_val])
    else:
        print(f"Gene {gene} not found in both datasets")

print(f"Healthy mean expression: {healthy_expr.mean()}, Unhealthy mean expression: {unhealthy_expr.mean()}")

# Results
de_results = pd.DataFrame(results, columns=['Gene', 'Log2FoldChange', 'PValue'])
de_results.to_csv('de_results_optimized.csv', index=False)
