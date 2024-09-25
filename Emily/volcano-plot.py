import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from gseapy import enrichr, barplot


de_results = pd.read_csv('de_results_optimized.csv')

print(f"Making Volcano plot...")
plt.figure(figsize=(10, 6))

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
plt.tight_layout()
plt.show()

significant_gene_names = significant_genes['Gene'].tolist()
print("Systematically significant genes:")
print(significant_gene_names)

#significant_gene_list = top_50_genes['Gene'].tolist()

significant_genes.to_csv('significant_genes.csv', index=False)



