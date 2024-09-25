import pandas as pd
import numpy as np
import mygene as mygene
import matplotlib.pyplot as plt


# Load the TSV file into a DataFrame
df = pd.read_csv("SRP158491/SRP158491.tsv", sep='\t')

nan_genes_before = df['Gene'].isna().sum()
total_genes_before = len(df['Gene'])

mg = mygene.MyGeneInfo()
query_result = mg.querymany(df['Gene'].tolist(), scopes='ensembl.gene',
                           fields='symbol', species='human')
ensembl_to_hugo = {item['query']: item.get('symbol', '') for item in query_result}
df['Gene'] = df['Gene'].map(ensembl_to_hugo)

df = df.drop_duplicates(subset='Gene', keep='first')
df.to_csv('SRP158491_converted.tsv', sep='\t', index=False)