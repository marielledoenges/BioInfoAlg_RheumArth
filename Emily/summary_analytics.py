import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

#STEP 1 - BASIC DATA INFO
df = pd.read_csv("SRP158491_converted.tsv", sep='\t')

df = df.dropna(subset=['Gene'])
df = df.drop_duplicates(subset='Gene', keep='first')

matrix_size = df.shape
print(f"Expression matrix size: {matrix_size}")

num_genes = df['Gene'].nunique()
print(f"Number of genes: {num_genes}")

median_expression = df.iloc[:, 1:].median(axis=1)

plt.figure(figsize=(10, 6))
sns.kdeplot(median_expression.clip(0, 12), fill=True)
plt.title("Density Plot of Per-Gene Median Expression Ranges")
plt.xlabel("Median Expression")
plt.ylabel("Density")
plt.savefig("median_expression_density_plot.png")
plt.show()

plt.figure(figsize=(10, 6))
sns.kdeplot(median_expression[median_expression > 0], color='blue', fill=True)
plt.title('Histogram of log counts ignoring 0 values')
plt.ylabel('Density')
plt.xlabel('Log Counts')
plt.show()

summary = f"""
Expression matrix size: {matrix_size}
Number of genes: {num_genes}
Median expression ranges:{median_expression.max()} - {median_expression.min()}
Mean: {median_expression.mean()}
Standard Deviation: {median_expression.std()}
"""

print(summary)

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

print(healthy_gene_expression.head)
print(unhealthy_gene_expression.head)




#PCA Plot
from sklearn.decomposition import PCA

pca = PCA(n_components=2)
pca.fit(df.iloc[:, 1:].transpose())

healthy_pca = pca.transform(healthy_gene_expression.iloc[:, 1:].transpose())
unhealthy_pca = pca.transform(unhealthy_gene_expression.iloc[:, 1:].transpose())


print("ra non treatment vs healthy")
plt.scatter(healthy_pca[:, 0], healthy_pca[:, 1], c='blue', label='Healthy')
plt.scatter(unhealthy_pca[:, 0], unhealthy_pca[:, 1], c='red', label='Unhealthy')
plt.xlabel('PCA Component 1')
plt.ylabel('PCA Component 2')
plt.title('PCA of Healthy and Unhealthy Individuals')
plt.legend()
plt.show()


