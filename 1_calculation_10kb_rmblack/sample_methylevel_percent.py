import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
MDD_dir = "/bioccNSW/MyXu/MDD/MDDs"
genome_size_file = "/bioccNSW/MyXu/simulation_read_level/sim/hg38.chrome.size"
output_csv = "/bioccNSW/MyXu/MDD/MDD_coverage_percentage.csv"
output_plot = "/bioccNSW/MyXu/MDD/MDD_coverage_plot.pdf"
def get_genome_size(genome_size_file):
    genome_df = pd.read_csv(genome_size_file, sep='\t', header=None, names=['chr', 'length'])
    return genome_df['length'].sum()
def calculate_MDD_length(bed_file):
    df = pd.read_csv(bed_file, sep='\t', header=None, names=['chr', 'start', 'end'])
    df['length'] = df['end'] - df['start']
    return df['length'].sum()
def main():
    genome_length = get_genome_size(genome_size_file)
    print(f": {genome_length}")
    results = []
    for file in os.listdir(MDD_dir):
        if file.endswith(".bed"):
            sample_path = os.path.join(MDD_dir, file)
            sample_name = file.replace(".bed", "")
            MDD_length = calculate_MDD_length(sample_path)
            percentage = (MDD_length / genome_length) * 100
            sample_type = 'Normal' if 'Nonmalignant' in file else 'Tumor'
            results.append({'Sample': sample_name, 'MDD Percentage': percentage, 'Type': sample_type})
            print(f"{sample_name}: {percentage:.4f}%")
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f": {output_csv}")
    plt.figure(figsize=(16, 9))
    sns.barplot(x='Sample', y='MDD Percentage', hue='Type', data=df, palette={'Normal': '#3498db', 'Tumor': '#e74c3c'})
    plt.xticks(rotation=45, ha='right')
    plt.title('MDD Coverage Percentage in Genome')
    plt.xlabel('Sample')
    plt.ylabel('MDD Percentage (%)')
    plt.tight_layout()
    plt.savefig(output_plot)
    print(f": {output_plot}")
if __name__ == '__main__':
    main()