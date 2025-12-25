import os
import pandas as pd
input_folder = '/student1/MDD/ESCA_rmblackList_cov5'
output_folder = '/student1/MDD/ESCA_10kb_methylation_2'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
for filename in os.listdir(input_folder):
    if filename.endswith('.bed'):
        file_path = os.path.join(input_folder, filename)
        df = pd.read_csv(file_path, sep='\t', header=None, names=['chromosome', 'start', 'end', 'methylation_level', 'reads'])
        df['start'] = df['start'].astype(int)
        df['end'] = df['end'].astype(int)
        chromosomes = df['chromosome'].unique()
        for chromosome in chromosomes:
            chromosome_df = df[df['chromosome'] == chromosome]
            chromosome_df.sort_values(by=['start'], inplace=True)
            start = chromosome_df['start'].iloc[0]
            end = start + 10000
            results = []
            while end < chromosome_df['end'].iloc[-1]:
                region_df = chromosome_df[(chromosome_df['start'] >= start) & (chromosome_df['end'] <= end)]
                if not region_df.empty:
                    mean_methylation = region_df['methylation_level'].mean()
                    results.append([chromosome, start, end, mean_methylation])
                start = end
                end += 10000
            region_df = chromosome_df[(chromosome_df['start'] >= start)]
            if not region_df.empty:
                mean_methylation = region_df['methylation_level'].mean()
                results.append([chromosome, start, region_df['end'].iloc[-1], mean_methylation])
        output_file_path = os.path.join(output_folder, filename)
        output_df = pd.DataFrame(results, columns=['chromosome', 'start', 'end', 'methylation_level'])
        output_df.to_csv(output_file_path, sep='\t', index=False, header=False)
print("！")