import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

de_results = pd.read_csv('de_results_optimized.csv')

df = pd.read_csv("SRP158491_converted.tsv", sep='\t')

df = df.dropna(subset=['Gene'])
df = df.drop_duplicates(subset='Gene', keep='first')

patient_df = pd.read_csv("SRP158491/metadata_SRP158491.tsv", sep='\t')

# healthy vs unhealthy samples based on metadata
healthy_rows = patient_df[patient_df['refinebio_disease'] == 'healthy']
unhealthy_rows = patient_df[patient_df['refinebio_disease'] == 'ra non treatment']

healthy_columns = healthy_rows['refinebio_accession_code'].tolist()
unhealthy_columns = unhealthy_rows['refinebio_accession_code'].tolist()

healthy_gene_expression = df[['Gene'] + healthy_columns]
unhealthy_gene_expression = df[['Gene'] + unhealthy_columns]

significant_genes = de_results[(de_results['PValue'] < 0.05) & (abs(de_results['Log2FoldChange']) > 0.5)]
print(f"Number of significantly differentially expressed genes: {len(significant_genes)}")

significant_gene_list = significant_genes['Gene'].tolist()

significant_healthy_expression = healthy_gene_expression[healthy_gene_expression['Gene'].isin(significant_gene_list)]
significant_unhealthy_expression = unhealthy_gene_expression[unhealthy_gene_expression['Gene'].isin(significant_gene_list)]

combined_expression = pd.concat([significant_healthy_expression.set_index('Gene'), 
                                 significant_unhealthy_expression.set_index('Gene')], axis=1)

plt.figure(figsize=(12, 8))
sns.heatmap(combined_expression, cmap='vlag', linewidths=0.5)
plt.title('Heatmap of Significantly Differentially Expressed Genes')
plt.tight_layout()
plt.show()
