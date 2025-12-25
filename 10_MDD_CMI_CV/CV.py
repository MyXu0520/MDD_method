# import os
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# input_dir = '/student1/MDD/ESCA_MDD_meth'
# output_dir = '/student1/MDD/CV'
# os.makedirs(output_dir, exist_ok=True)
# all_data = pd.DataFrame()
# bed_files = [f for f in os.listdir(input_dir) if f.endswith('.bed')]
# for bed_file in bed_files:
#     file_path = os.path.join(input_dir, bed_file)
#     data = pd.read_csv(file_path, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Methylation_Level'])
#     sample_name = os.path.splitext(bed_file)[0]
#     data['Sample'] = sample_name
#     all_data = pd.concat([all_data, data], ignore_index=True)
# grouped = all_data.groupby(['Chromosome', 'Start', 'End'])['Methylation_Level'].agg(['mean', 'std']).reset_index()
# grouped['CV'] = grouped['std'] / grouped['mean']
# cv_threshold = np.percentile(grouped['CV'], 10)
# chromosomes = grouped['Chromosome'].unique()
# sns.set(style="whitegrid", context="notebook", palette="deep")
# for chromosome in chromosomes:
#     chromosome_data = grouped[grouped['Chromosome'] == chromosome]
#     chromosome_data = chromosome_data.sort_values(by=['Start'])
#     min_10_percent_indices = chromosome_data['CV'].nsmallest(int(len(chromosome_data) * 0.1)).index
#     min_10_percent_data = chromosome_data.loc[min_10_percent_indices]
#     plt.figure(figsize=(14, 8))
#     sns.lineplot(x='Start', y='CV', data=chromosome_data, marker='o', markersize=6, ci=None, label='All Regions')
#     sns.lineplot(x='Start', y='CV', data=min_10_percent_data, marker='o', markersize=6, color='red', label='Top 10% Lowest CV Regions')
#     plt.axhline(y=cv_threshold, color='gray', linestyle='--', label=f'10% CV Threshold = {cv_threshold:.2f}')
#     plt.title(f'Variation Coefficient of Methylation Levels by Position on Chromosome {chromosome}', fontsize=16)
#     plt.xlabel('Start Position', fontsize=14)
#     plt.ylabel('Variation Coefficient (CV)', fontsize=14)
#     plt.xticks(rotation=45, ha='right')
#     plt.grid(True, linestyle='--', alpha=0.7)
#     plt.yticks(np.arange(0, chromosome_data['CV'].max() + 0.1, 0.1))
#     plt.gca().spines['top'].set_visible(False)
#     plt.gca().spines['right'].set_visible(False)
#     plt.legend()
#     image_path = os.path.join(output_dir, f'cv_chromosome_{chromosome}.pdf')
#     plt.savefig(image_path, dpi=300, bbox_inches='tight')
#     plt.close()
# import os
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# input_dir = '/student1/MDD/ESCA_MDD_meth'
# output_dir = '/student1/MDD/CV'
# bed_output_dir = os.path.join(output_dir, 'top_10_percent_CV_regions')
# os.makedirs(output_dir, exist_ok=True)
# os.makedirs(bed_output_dir, exist_ok=True)
# all_data = pd.DataFrame()
# bed_files = [f for f in os.listdir(input_dir) if f.endswith('.bed')]
# for bed_file in bed_files:
#     file_path = os.path.join(input_dir, bed_file)
#     data = pd.read_csv(file_path, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Methylation_Level'])
#     sample_name = os.path.splitext(bed_file)[0]
#     data['Sample'] = sample_name
#     all_data = pd.concat([all_data, data], ignore_index=True)
# grouped = all_data.groupby(['Chromosome', 'Start', 'End'])['Methylation_Level'].agg(['mean', 'std']).reset_index()
# grouped['CV'] = grouped['std'] / grouped['mean']
# cv_threshold = np.percentile(grouped['CV'], 10)
# chromosomes = grouped['Chromosome'].unique()
# sns.set(style="whitegrid", context="notebook", palette="deep")
# for chromosome in chromosomes:
#     chromosome_data = grouped[grouped['Chromosome'] == chromosome]
#     chromosome_data = chromosome_data.sort_values(by=['Start'])
#     plt.figure(figsize=(14, 8))
#     sns.lineplot(x='Start', y='CV', data=chromosome_data, marker='o', markersize=6, ci=None, label='CV')
#     plt.axhline(y=cv_threshold, color='gray', linestyle='--', label=f'10% CV Threshold ({cv_threshold:.4f})')
#     low_cv_data = chromosome_data[chromosome_data['CV'] <= cv_threshold]
#     sns.scatterplot(x='Start', y='CV', data=low_cv_data, color='red', marker='o', s=50, label='Top 10% Low CV Regions')
#     plt.title(f'Variation Coefficient of Methylation Levels by Position on Chromosome {chromosome}', fontsize=16)
#     plt.xlabel('Start Position', fontsize=14)
#     plt.ylabel('Variation Coefficient (CV)', fontsize=14)
#     plt.xticks(rotation=45, ha='right')
#     plt.grid(True, linestyle='--', alpha=0.7)
#     plt.yticks(np.arange(0, max(chromosome_data['CV'].max(), cv_threshold) + 0.1, 0.1))
#     plt.gca().spines['top'].set_visible(False)
#     plt.gca().spines['right'].set_visible(False)
#     plt.legend()
#     image_path = os.path.join(output_dir, f'cv_chromosome_{chromosome}.pdf')
#     plt.savefig(image_path, dpi=300, bbox_inches='tight')
#     plt.close()
# top_10_percent_cv_regions = grouped[grouped['CV'] <= cv_threshold]
# for index, row in top_10_percent_cv_regions.iterrows():
#     chromosome = row['Chromosome']
#     start = row['Start']
#     end = row['End']
#     cv = row['CV']
#     bed_file_path = os.path.join(bed_output_dir, f'{chromosome}_{start}_{end}_CV_{cv:.4f}.bed')
#     with open(bed_file_path, 'w') as bed_file:
#         bed_file.write(f'{chromosome}\t{start}\t{end}\n')
# import os
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# input_dir = '/student1/MDD/ESCA_MDD_meth'
# output_dir = '/student1/MDD/CV'
# output_bed_dir = os.path.join(output_dir, 'top_10_percent_cv_regions')
# os.makedirs(output_dir, exist_ok=True)
# os.makedirs(output_bed_dir, exist_ok=True)
# all_data = pd.DataFrame()
# bed_files = [f for f in os.listdir(input_dir) if f.endswith('.bed')]
# for bed_file in bed_files:
#     file_path = os.path.join(input_dir, bed_file)
#     data = pd.read_csv(file_path, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Methylation_Level'])
#     sample_name = os.path.splitext(bed_file)[0]
#     data['Sample'] = sample_name
#     all_data = pd.concat([all_data, data], ignore_index=True)
# grouped = all_data.groupby(['Chromosome', 'Start', 'End'])['Methylation_Level'].agg(['mean', 'std']).reset_index()
# grouped['CV'] = grouped['std'] / grouped['mean']
# cv_threshold = np.percentile(grouped['CV'], 10)
# chromosomes = grouped['Chromosome'].unique()
# sns.set(style="whitegrid", context="notebook", palette="deep")
# for chromosome in chromosomes:
#     chromosome_data = grouped[grouped['Chromosome'] == chromosome]
#     chromosome_data = chromosome_data.sort_values(by=['Start'])
#     plt.figure(figsize=(14, 8))
#     sns.lineplot(x='Start', y='CV', data=chromosome_data, marker='o', markersize=6, ci=None, label='Overall CV')
#     plt.axhline(y=cv_threshold, color='gray', linestyle='--', label=f'10% CV Threshold ({cv_threshold:.4f})')
#     low_cv_data = chromosome_data[chromosome_data['CV'] <= cv_threshold]
#     sns.scatterplot(x='Start', y='CV', data=low_cv_data, color='red', marker='o', s=50, label='Top 10% Low CV Regions')
#     plt.title(f'Variation Coefficient of Methylation Levels by Position on Chromosome {chromosome}', fontsize=16)
#     plt.xlabel('Start Position', fontsize=14)
#     plt.ylabel('Variation Coefficient (CV)', fontsize=14)
#     plt.xticks(rotation=45, ha='right')
#     plt.grid(True, linestyle='--', alpha=0.7)
#     plt.yticks(np.arange(0, chromosome_data['CV'].max() + 0.1, 0.1))
#     plt.gca().spines['top'].set_visible(False)
#     plt.gca().spines['right'].set_visible(False)
#     plt.legend()
#     image_path = os.path.join(output_dir, f'cv_chromosome_{chromosome}.pdf')
#     plt.savefig(image_path, dpi=300, bbox_inches='tight')
#     plt.close()
# top_10_percent_cv_regions = grouped[grouped['CV'] <= cv_threshold]
# top_10_percent_cv_regions.to_csv(os.path.join(output_bed_dir, 'top_10_percent_cv_regions.bed'), sep='\t', index=False, header=False)
# for chromosome in top_10_percent_cv_regions['Chromosome'].unique():
#     chromosome_regions = top_10_percent_cv_regions[top_10_percent_cv_regions['Chromosome'] == chromosome]
#     output_file_path = os.path.join(output_bed_dir, f'top_10_percent_cv_regions_chr{chromosome}.bed')
#     chromosome_regions.to_csv(output_file_path, sep='\t', index=False, header=False)
# import os
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
# input_dir = '/student1/MDD/ESCA_MDD_meth'
# output_dir = '/student1/MDD/CV'
# output_bed_dir = os.path.join(output_dir, 'top_10_percent_CV_regions')
# os.makedirs(output_dir, exist_ok=True)
# os.makedirs(output_bed_dir, exist_ok=True)
# all_data = pd.DataFrame()
# bed_files = [f for f in os.listdir(input_dir) if f.endswith('.bed')]
# for bed_file in bed_files:
#     file_path = os.path.join(input_dir, bed_file)
#     data = pd.read_csv(file_path, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Methylation_Level'])
#     sample_name = os.path.splitext(bed_file)[0]
#     data['Sample'] = sample_name
#     all_data = pd.concat([all_data, data], ignore_index=True)
# grouped = all_data.groupby(['Chromosome', 'Start', 'End'])['Methylation_Level'].agg(['mean', 'std']).reset_index()
# grouped['CV'] = grouped['std'] / grouped['mean']
# cv_threshold = np.percentile(grouped['CV'], 10)
# chromosomes = grouped['Chromosome'].unique()
# sns.set(style="whitegrid", context="notebook", palette="deep")
# for chromosome in chromosomes:
#     chromosome_data = grouped[grouped['Chromosome'] == chromosome]
#     chromosome_data = chromosome_data.sort_values(by=['Start'])
#     plt.figure(figsize=(14, 8))
#     sns.lineplot(x='Start', y='CV', data=chromosome_data, marker='o', markersize=6, ci=None, label='All Regions')
#     plt.axhline(y=cv_threshold, color='r', linestyle='--', label=f'10% CV Threshold ({cv_threshold:.4f})')
#     top_10_percent_cv_regions = chromosome_data.nsmallest(int(len(chromosome_data) * 0.1), 'CV')
#     sns.lineplot(x='Start', y='CV', data=top_10_percent_cv_regions, marker='o', markersize=6, ci=None, label='Top 10% Lowest CV Regions', color='green')
#     plt.title(f'Variation Coefficient of Methylation Levels by Position on Chromosome {chromosome}', fontsize=16)
#     plt.xlabel('Start Position', fontsize=14)
#     plt.ylabel('Variation Coefficient (CV)', fontsize=14)
#     plt.xticks(rotation=45, ha='right')
#     plt.grid(True, linestyle='--', alpha=0.7)
#     plt.legend()
#     plt.gca().spines['top'].set_visible(False)
#     plt.gca().spines['right'].set_visible(False)
#     image_path = os.path.join(output_dir, f'cv_chromosome_{chromosome}.pdf')
#     plt.savefig(image_path, dpi=300, bbox_inches='tight')
#     plt.close()
#     top_10_percent_cv_regions_bed = top_10_percent_cv_regions[['Chromosome', 'Start', 'End']].copy()
#     top_10_percent_cv_regions_bed.to_csv(os.path.join(output_bed_dir, f'top_10_percent_cv_chromosome_{chromosome}.bed'), sep='\t', index=False, header=False)
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
input_dir = '/student1/MDD/ESCA_MDD_meth'
output_dir = '/student1/MDD/CV'
top_cv_dir = os.path.join(output_dir, 'top_10_percent_cv')
os.makedirs(output_dir, exist_ok=True)
os.makedirs(top_cv_dir, exist_ok=True)
all_data = pd.DataFrame()
bed_files = [f for f in os.listdir(input_dir) if f.endswith('.bed')]
for bed_file in bed_files:
    file_path = os.path.join(input_dir, bed_file)
    data = pd.read_csv(file_path, sep='\t', header=None, names=['Chromosome', 'Start', 'End', 'Methylation_Level'])
    sample_name = os.path.splitext(bed_file)[0]
    data['Sample'] = sample_name
    all_data = pd.concat([all_data, data], ignore_index=True)
grouped = all_data.groupby(['Chromosome', 'Start', 'End'])['Methylation_Level'].agg(['mean', 'std']).reset_index()
grouped['CV'] = grouped['std'] / grouped['mean']
cv_threshold = np.percentile(grouped['CV'], 10)
top_10_percent_cv = grouped[grouped['CV'] >= cv_threshold]
chromosomes = grouped['Chromosome'].unique()
sns.set(style="whitegrid", context="notebook", palette="deep")
for chromosome in chromosomes:
    chromosome_data = grouped[grouped['Chromosome'] == chromosome]
    chromosome_data = chromosome_data.sort_values(by=['Start'])
    plt.figure(figsize=(14, 8))
    sns.lineplot(x='Start', y='CV', data=chromosome_data, marker='o', markersize=6, ci=None, label='All Regions')
    plt.axhline(y=cv_threshold, color='r', linestyle='--', label=f'10% CV Threshold ({cv_threshold:.4f})')
    bottom_10_percent_cv = chromosome_data.nsmallest(int(len(chromosome_data) * 0.1), 'CV')
    sns.lineplot(x='Start', y='CV', data=bottom_10_percent_cv, marker='o', markersize=6, ci=None, label='Bottom 10% CV Regions', color='green')
    plt.title(f'Variation Coefficient of Methylation Levels by Position on Chromosome {chromosome}', fontsize=16)
    plt.xlabel('Start Position', fontsize=14)
    plt.ylabel('Variation Coefficient (CV)', fontsize=14)
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.yticks(np.arange(0, chromosome_data['CV'].max() + 0.1, 0.1))
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    image_path = os.path.join(output_dir, f'cv_chromosome_{chromosome}.pdf')
    plt.savefig(image_path, dpi=300, bbox_inches='tight')
    plt.close()
for _, row in top_10_percent_cv.iterrows():
    chromosome = row['Chromosome']
    start = row['Start']
    end = row['End']
    cv = row['CV']
    bed_line = f'{chromosome}\t{start}\t{end}\tCV={cv:.4f}'
    bed_file_path = os.path.join(top_cv_dir, f'{chromosome}_top_10_percent_cv.bed')
    with open(bed_file_path, 'a') as bed_file:
        bed_file.write(bed_line + '\n')
print("，。，，10%The region of the coefficient of variation has also been saved asBEDDoc.。")