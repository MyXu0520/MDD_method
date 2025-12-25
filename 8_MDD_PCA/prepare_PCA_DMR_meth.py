import os
import pandas as pd
input_dir = '/student1/MDD/ESCA_rmblackList_cov5'
result_dir = '/student1/MDD/ESCA_MDD_meth'
shared_regions_file = '/student1/MDD/result/1/ESCA_sharedMDDs.bed'
output_matrix_file = os.path.join(result_dir, 'methylation_matrix.csv')
os.makedirs(result_dir, exist_ok=True)
shared_regions = pd.read_csv(shared_regions_file, sep='\t', header=None, names=['chrom', 'start', 'end'])
sample_results = {}
for filename in os.listdir(input_dir):
    if filename.endswith('.bed'):
        sample_name = os.path.splitext(filename)[0]
        print(sample_name)
        sample_df = pd.read_csv(os.path.join(input_dir, filename), sep='\t', header=None, names=['chrom', 'start', 'end', 'meth_level', 'reads'])
        region_stats = {}
        for _, region in shared_regions.iterrows():
            chrom, start, end = region
            # print(region)
            region_df = sample_df[(sample_df['chrom'] == chrom) &
                                  (sample_df['start'].between(start, end, inclusive='left')) &
                                  (sample_df['end'].between(start, end, inclusive='right'))]
            num_sites = len(region_df)
            if num_sites > 0:
                avg_meth_level = region_df['meth_level'].mean()
            else:
                avg_meth_level = 0
            region_stats[(chrom, start, end)] = avg_meth_level
            # print(avg_meth_level)
        sample_results[sample_name] = region_stats
        print(sample_name)
for sample_name, results in sample_results.items():
    output_data = [[chrom, start, end, f"{avg_meth:.4f}"] for (chrom, start, end), avg_meth in results.items()]
    df = pd.DataFrame(output_data, columns=['chrom', 'start', 'end', 'avg_meth_level'])
    output_file = os.path.join(result_dir, f'{sample_name}_avg_meth.bed')
    df.to_csv(output_file, sep='\t', index=False, header=False)
    print(sample_name)
matrix_data = []
region_order = sorted(sample_results[next(iter(sample_results))].keys())
header = ['Region'] + list(sample_results.keys())
for region_key in region_order:
    row_data = [f"{region_key[0]}:{region_key[1]}-{region_key[2]}"]
    for sample_name in sample_results.keys():
        avg_meth = sample_results[sample_name].get(region_key, 0)
        row_data.append(f"{avg_meth:.4f}")
    matrix_data.append(row_data)
matrix_df = pd.DataFrame(matrix_data, columns=header)
matrix_df.to_csv(output_matrix_file, index=False)