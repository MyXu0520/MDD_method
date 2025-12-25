import pandas as pd
bed_df = pd.read_csv('/student1/MDD/result/1/EAC_specificMDDs.bed', sep='\t', header=None, names=['chromosome', 'start', 'end'])
gene_df = pd.read_csv('gene_anno.csv', header=None, names=['chromosome', 'gene_start', 'gene_end', 'gene_id'])
results_df = pd.DataFrame(columns=['BED_chromosome', 'BED_start', 'BED_end', 'Gene_ID'])
# for _, bed_row in bed_df.iterrows():
#     bed_chromosome = bed_row['chromosome']
#     bed_start = bed_row['start']
#     bed_end = bed_row['end']
#     mask = (gene_df['chromosome'] == bed_chromosome) & \
#            (gene_df['gene_start'] <= bed_end) & \
#            (gene_df['gene_end'] >= bed_start)
#     covered_genes = gene_df.loc[mask]
#     for _, gene_row in covered_genes.iterrows():
#         results_df = results_df.append({
#             'BED_chromosome': bed_chromosome,
#             'BED_start': bed_start,
#             'BED_end': bed_end,
#             'Gene_ID': gene_row['gene_id']
#         }, ignore_index=True)
# results_df.to_csv('covered_genes.csv', index=False)
results_rows = []
for _, bed_row in bed_df.iterrows():
    bed_chromosome = bed_row['chromosome']
    bed_start = bed_row['start']
    bed_end = bed_row['end']
    mask = (gene_df['chromosome'] == bed_chromosome) & \
           (gene_df['gene_start'] <= bed_end) & \
           (gene_df['gene_end'] >= bed_start)
    covered_genes = gene_df.loc[mask]
    for _, gene_row in covered_genes.iterrows():
        results_rows.append({
            'BED_chromosome': bed_chromosome,
            'BED_start': bed_start,
            'BED_end': bed_end,
            'Gene_ID': gene_row['gene_id']
        })
results_df = pd.DataFrame(results_rows)
results_df.to_csv('covered_genes.csv', index=False)