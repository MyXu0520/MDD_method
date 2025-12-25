# import os
# import pandas as pd
# import matplotlib.pyplot as plt
# work_dir = '/student1/MDD'
# bed_files = [f for f in os.listdir(work_dir) if f.endswith('.bed')]
# for bed_file in bed_files:
#     df = pd.read_csv(os.path.join(work_dir, bed_file), sep='\t', header=None, names=['chromosome', 'start', 'end', 'methylation_level'])
#     methylation_levels = df['methylation_level']
#     plt.figure(figsize=(10, 6))
#     plt.title(f'Methylation Level Distribution for {os.path.splitext(bed_file)[0]}')
#     plt.xlabel('Methylation Level')
#     plt.ylabel('Frequency')
#     plt.grid(True)
#     output_file = os.path.join(work_dir, f'{os.path.splitext(bed_file)[0]}_methylation_histogram.pdf')
#     plt.savefig(output_file)
#     plt.close()
# print("Histograms have been generated and saved.")
import os
import pandas as pd
import matplotlib.pyplot as plt
work_dir = '/student1/MDD/ESCA_10kb_methylation'
color_dict = {
    "ESCC_Nonmalignant": "#4489C8",
    "GEJ_Nonmalignant": "#EE7C79",
    "ESCC": "#008F91",
    "EAC": "#FFCD44",
    "GEJ": "#FFCD44"
}
bed_files = [f for f in os.listdir(work_dir) if f.endswith('.bed')]
for bed_file in bed_files:
    file_name = os.path.splitext(bed_file)[0]
    for sample_type, color in color_dict.items():
        if sample_type in file_name:
            plot_color = color
            break
    else:
        plot_color = '#808080'
    df = pd.read_csv(os.path.join(work_dir, bed_file), sep='\t', header=None, names=['chromosome', 'start', 'end', 'methylation_level'])
    methylation_levels = df['methylation_level']
    plt.figure(figsize=(10, 6))
    plt.hist(methylation_levels, bins=30, edgecolor='black', alpha=0.7, color=plot_color)
    plt.title(f'Methylation Level Distribution for {file_name}')
    plt.xlabel('Methylation Level')
    plt.ylabel('Frequency')
    plt.grid(True)
    output_file = os.path.join(work_dir, f'{file_name}_methylation_histogram.pdf')
    plt.savefig(output_file)
    plt.close()
print("Colored histograms have been generated and saved.")