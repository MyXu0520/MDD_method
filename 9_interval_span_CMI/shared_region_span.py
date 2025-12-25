import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
input_file = '/student1/MDD/result/1/ESCA_sharedMDDs.bed'
output_dir = '/student1/MDD/ESCA_MDD_meth'
os.makedirs(output_dir, exist_ok=True)
df = pd.read_csv(input_file, sep='\t', header=None, names=['chromosome', 'start', 'end'])
df['start'] = df['start'].astype(int)
df['end'] = df['end'].astype(int)
df_sorted = df.sort_values(by=['chromosome', 'start']).reset_index(drop=True)
df_sorted['span'] = df_sorted['end'] - df_sorted['start'] + 1
df_sorted['prev_end'] = df_sorted.groupby('chromosome')['end'].shift(1)
df_intervals = df_sorted.dropna(subset=['prev_end'])
df_intervals['interval'] = df_intervals['start'] - df_intervals['prev_end']
Q1 = df_intervals['interval'].quantile(0.25)
Q3 = df_intervals['interval'].quantile(0.75)
IQR = Q3 - Q1
df_intervals_filtered = df_intervals[~((df_intervals['interval'] < (Q1 - 1.5 * IQR)) | (df_intervals['interval'] > (Q3 + 1.5 * IQR)))]
plt.figure(figsize=(10, 6))
plt.hist(df_intervals_filtered['interval'], bins=50, edgecolor='black', alpha=0.7)
plt.xlabel('Interval Size (bp)')
plt.ylabel('Frequency')
plt.title('Interval Size Frequency Distribution Histogram (Outliers Removed)')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'interval_size_frequency_histogram.pdf'))
plt.close()
Q1_span = df_sorted['span'].quantile(0.25)
Q3_span = df_sorted['span'].quantile(0.75)
IQR_span = Q3_span - Q1_span
df_sorted_filtered = df_sorted[~((df_sorted['span'] < (Q1_span - 1.5 * IQR_span)) | (df_sorted['span'] > (Q3_span + 1.5 * IQR_span)))]
plt.figure(figsize=(10, 6))
plt.hist(df_sorted_filtered['span'], bins=50, edgecolor='black', alpha=0.7)
plt.xlabel('Span Size (bp)')
plt.ylabel('Frequency')
plt.title('Span Size Frequency Distribution Histogram (Outliers Removed)')
plt.grid(axis='y', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, 'span_size_frequency_histogram.pdf'))
plt.close()