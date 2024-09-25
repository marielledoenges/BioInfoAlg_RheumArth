import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load the data into a DataFrame
df = pd.read_csv("SRP158491_converted.tsv", sep='\t')

# Count NaN genes after conversion
 # This should be the same as before

nan_genes_after = df['Gene'].isna().sum()
total_genes_after = len(df['Gene']) 

# Display the size of the expression matrix
matrix_size = df.shape
print(f"Expression matrix size: {matrix_size}")

# Calculate the number of genes
num_genes = df['Gene'].nunique()
print(f"Number of genes: {num_genes}")

# Calculate per-gene median expression ranges
median_expression = df.iloc[:, 1:].median(axis=1)

# Create a density plot of the median expression ranges
plt.figure(figsize=(10, 6))
sns.kdeplot(median_expression, fill=True)
plt.title("Density Plot of Per-Gene Median Expression Ranges")
plt.xlabel("Median Expression")
plt.ylabel("Density")
plt.savefig("median_expression_density_plot.png")
plt.show()

# Summarize findings
summary = f"""
Expression matrix size: {matrix_size}
Number of genes: {num_genes}
Median expression ranges:{median_expression.max()} - {median_expression.min()}
Mean: {median_expression.mean()}
Standard Deviation: {median_expression.std()}
"""

print(summary)

# Save the summary to a text file
with open("summary.txt", "w") as file:
    file.write(summary)

# PART 2 #

#csplit dataset into health vs unhealthy
patient_df = pd.read_csv("SRP158491/metadata_SRP158491.tsv", sep='\t')

healthy_rows = patient_df[patient_df['refinebio_disease'] == 'healthy']
unhealthy_rows = patient_df[patient_df['refinebio_disease'] == 'ra non treatment']

healthy_columns = healthy_rows['refinebio_accession_code'].tolist()
unhealthy_columns = unhealthy_rows['refinebio_accession_code'].tolist()

healthy_gene_expression = df[healthy_columns]
unhealthy_gene_expression = df[unhealthy_columns]


#treaded unhealthy samples
#ra treatment
#mtx_treated = patient_df[patient_df['refinebio_disease'] == 'ra mtx treatment']
#tcz_treated = patient_df[patient_df['refinebio_disease'] == 'ra tcz treatment']
#ifx_treated = patient_df[patient_df['refinebio_disease'] == 'ra ifx treatment']
#all_treated = pd.concat([mtx_treated, tcz_treated, ifx_treated], axis=1)

#PCA Plot
from sklearn.decomposition import PCA

pca = PCA(n_components=2)
pca.fit(df.iloc[:, 1:].transpose())

healthy_pca = pca.transform(healthy_gene_expression.iloc[:, 1:].transpose())
unhealthy_pca = pca.transform(unhealthy_gene_expression.iloc[:, 1:].transpose())
#treated_pca = pca_unhealthy.fit_transform(all_treated.iloc[:, 1:])
#plotting
plt.figure(figsize=(12, 8))
plt.scatter(healthy_pca[:, 0], healthy_pca[:, 1], c='blue', label='Healthy')
plt.scatter(unhealthy_pca[:, 0], unhealthy_pca[:, 1], c='red', label='Unhealthy')
plt.xlabel('PCA Component 1')
plt.ylabel('PCA Component 2')
plt.title('PCA of Healthy and Unhealthy Individuals')
plt.legend()
plt.show()

