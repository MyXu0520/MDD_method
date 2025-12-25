import pandas as pd
import os
input_dir = "/student1/MDD/result/1"
output_dir = "/student1/MDD/Shared_CMI"
os.makedirs(output_dir, exist_ok=True)
esca_bed = pd.read_csv(os.path.join(input_dir, "ESCA_sharedMDDs.bed"), sep='\t', header=None, names=["chr", "start", "end"])
esca_bed["start"] = esca_bed["start"].astype(int)
esca_bed["end"] = esca_bed["end"].astype(int)
esca_bed.sort_values(by=["chr", "start"], inplace=True)
esca_bed['next_start'] = esca_bed.groupby('chr')['start'].shift(-1)
esca_bed['interval'] = esca_bed['next_start'] - esca_bed['end']
interval_threshold_min = 30000  # 30kb
interval_threshold_max = 300000  # 03mb
valid_intervals = esca_bed[(esca_bed['interval'] > interval_threshold_min) & (esca_bed['interval'] < interval_threshold_max)]
valid_intervals = valid_intervals[['chr', 'end', 'next_start']].rename(columns={"end": "start", "next_start": "end"})
valid_intervals = valid_intervals.dropna()
output_file = os.path.join(output_dir, "intervals_30kb_to_03mb.bed")
valid_intervals.to_csv(output_file, sep='\t', index=False, header=False, float_format='%d')
print(f"Intervals BED file saved to {output_file}")