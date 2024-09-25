import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference
import gseapy as gp

# Load your gene expression count data and metadata
metadata_df = pd.read_csv('SRP158491/metadata_SRP158491.tsv', sep='\t')
df = pd.read_csv('SRP158491_converted.tsv', sep='\t')

df = df.dropna()
df = df.drop_duplicates(subset='Gene', keep='first')

metadata_df_healthy = metadata_df.loc[metadata_df['refinebio_disease'] == 'healthy']
metadata_df_unhealthy = metadata_df.loc[metadata_df['refinebio_disease'] == 'ra non treatment']

metadata_df = pd.concat([metadata_df_healthy, metadata_df_unhealthy])

metadata_df

len(metadata_df['refinebio_accession_code'].unique())

deseq_df = df.set_index('Gene')

deseq_df = deseq_df.round(0)
deseq_df = deseq_df.astype(int)
deseq_df

deseq_df.transpose()
metadata_df.set_index('refinebio_accession_code', inplace=True)
metadata_df

inference = DefaultInference(n_cpus=4)
dds = DeseqDataSet(
    counts=deseq_df.transpose().loc[metadata_df.index],
    metadata=metadata_df,
    design_factors='refinebio_disease',
    refit_cooks=True,
    inference=inference
)

dds.deseq2()
stat_res = DeseqStats(dds, inference=inference)
stat_res.summary()

result_df = stat_res.results_df

results = result_df

results['significant'] = (results['padj'] < 0.05) & (abs(results['log2FoldChange']) > 1)

print("Significant genes:")
print(results.loc[results.significant])
results['nonsignificant'] = ~results['significant']

significant_genes_df = results[results['significant']]

significant_genes = significant_genes_df.index.tolist()

significant_genes_df.to_csv('significant_genes.csv')

rnk = results[['log2FoldChange']].copy()
rnk = rnk.dropna()
rnk = rnk[~rnk.index.duplicated(keep='first')]
rnk = rnk.sort_values('log2FoldChange', ascending=False)

rnk.reset_index(inplace=True)
rnk.columns = ['Gene', 'log2FoldChange']
rnk.to_csv('gene_rankings.rnk', sep='\t', index=False, header=False)
rnk_list = rnk.set_index('Gene')['log2FoldChange']
gene_sets = 'GO_Biological_Process_2023'

# Perform GSEA preranked analysis
pre_res = gp.prerank(
    rnk=rnk_list,
    gene_sets=gene_sets,
    processes=4,
    permutation_num=100,  # Adjust permutation number for accuracy; higher number increases computation time
    outdir='gsea_prerank_results',  # Output directory
    format='png',  # Output format for figures
    seed=6
)

# View the top enriched gene sets
print("Top enriched gene sets:")
print(pre_res.res2d.head())

# Save results to a CSV file
pre_res.res2d.to_csv('gsea_prerank_results.csv', index=False)