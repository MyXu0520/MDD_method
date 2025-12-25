# import pandas as pd
# import pybedtools
# from scipy.stats import hypergeom
# import matplotlib.pyplot as plt
# import venn
# def read_bed_file(file_path):
#     columns = ['chromosome', 'start', 'end']
#     df = pd.read_csv(file_path, sep='\t', header=None, names=columns)
#     return df
# def calculate_overlap(file1_path, file2_path):
#     df1 = read_bed_file(file1_path)
#     df2 = read_bed_file(file2_path)
#     bt1 = pybedtools.BedTool.from_dataframe(df1)
#     bt2 = pybedtools.BedTool.from_dataframe(df2)
#     intersection_df = intersection.to_dataframe(header=None, names=['chromosome', 'start', 'end'])
#     return len(intersection_df), df1.shape[0], df2.shape[0]
# def hypergeometric_test(overlap, total1, total2):
#     return p_value
# def plot_venn_diagram(set1, set2, overlap, output_file):
#     venn2 = venn.venn2([set(set1), set(set2)], ('Set1', 'Set2'))
#     plt.title("Venn Diagram of BED File Regions")
#     plt.savefig(output_file, format='png')
# def main():
#     file1_path = '/student1/MDD/meta/loci/ESCC_specificMDDs.bed'
#     file2_path = '/student1/MDD/meta/loci/ESCC_NT_DMR_bed.bed'
#     output_venn_file = 'venn_diagram.png'
#     overlap, total1, total2 = calculate_overlap(file1_path, file2_path)
#     df1 = read_bed_file(file1_path)
#     df2 = read_bed_file(file2_path)
#     set1 = {f"{row['chromosome']}:{row['start']}-{row['end']}" for _, row in df1.iterrows()}
#     set2 = {f"{row['chromosome']}:{row['start']}-{row['end']}" for _, row in df2.iterrows()}
#     p_value = hypergeometric_test(overlap, total1, total2)
#     print(f"Hypergeometric test p-value: {p_value}")
#     plot_venn_diagram(set1, set2, overlap, output_venn_file)
# if __name__ == "__main__":
#     main()
# bedtools intersect -a /student1/MDD/meta/loci/ESCA_sharedMDDs.bed -b /student1/MDD/meta/loci/EAC_ESCC_DMR_bed.bed -wa -wb > /student1/MDD/meta/loci/intersection1.bed
# bedtools intersect -a /student1/MDD/meta/loci/EAC_specificMDDs.bed -b /student1/MDD/meta/loci/EAC_SPEC_DMR.bed -wa -wb > /student1/MDD/meta/loci/intersection2.bed
# bedtools intersect -a /student1/MDD/meta/loci/ESCC_specificMDDs.bed -b /student1/MDD/meta/loci/ESCC_SPEC_DMR.bed -wa -wb > /student1/MDD/meta/loci/intersection3.bed
import pandas as pd
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import matplotlib_venn as venn
def count_bed_lines(file_path):
    with open(file_path, 'r') as file:
        return sum(1 for _ in file)
bed_file1 = '/student1/MDD/meta/loci/ESCA_sharedMDDs.bed'
bed_file2 = '/student1/MDD/meta/loci/EAC_ESCC_DMR_bed.bed'
intersection_file = '/student1/MDD/meta/loci/intersection1.bed'
len_file1 = count_bed_lines(bed_file1)
len_file2 = count_bed_lines(bed_file2)
len_intersection = count_bed_lines(intersection_file)
total_regions = len_file1 + len_file2 - len_intersection
k_successes = len_file1
n_draws = len_file2
p_value = 1 - hypergeom.cdf(len_intersection - 1, total_regions, k_successes, n_draws)
print(f"Total regions in file1: {len_file1}")
print(f"Total regions in file2: {len_file2}")
print(f"Total intersected regions: {len_intersection}")
print(f"Hypergeometric P-value: {p_value}")
# def read_bed_to_set(file_path):
#     with open(file_path, 'r') as file:
#         return set(tuple(line.strip().split()) for line in file)
# set1 = read_bed_to_set(bed_file1)
# set2 = read_bed_to_set(bed_file2)
# fig, ax = plt.subplots()
# venn_diagram = venn.venn2([set1, set2], ('ESCA_sharedMDDs.bed', 'EAC_ESCC_DMR_bed.bed'), ax=ax)
# plt.title("Venn Diagram of EAC_specificMDDs and EAC_SPEC_DMR.bed")
# plt.savefig('/student1/MDD/meta/loci/venn_diagram.png')
# plt.show()
def count_bed_lines(file_path):
    with open(file_path, 'r') as file:
        return sum(1 for _ in file)
bed_file1 = '/student1/MDD/meta/loci/EAC_specificMDDs.bed'
bed_file2 = '/student1/MDD/meta/loci/EAC_SPEC_DMR.bed'
intersection_file = '/student1/MDD/meta/loci/intersection2.bed'
len_file1 = count_bed_lines(bed_file1)
len_file2 = count_bed_lines(bed_file2)
len_intersection = count_bed_lines(intersection_file)
total_regions = len_file1 + len_file2 - len_intersection
k_successes = len_file1
n_draws = len_file2
p_value = 1 - hypergeom.cdf(len_intersection - 1, total_regions, k_successes, n_draws)
print(f"Total regions in file1: {len_file1}")
print(f"Total regions in file2: {len_file2}")
print(f"Total intersected regions: {len_intersection}")
print(f"Hypergeometric P-value: {p_value}")
def count_bed_lines(file_path):
    with open(file_path, 'r') as file:
        return sum(1 for _ in file)
bed_file1 = '/student1/MDD/meta/loci/ESCC_specificMDDs.bed'
bed_file2 = '/student1/MDD/meta/loci/ESCC_SPEC_DMR.bed'
intersection_file = '/student1/MDD/meta/loci/intersection3.bed'
len_file1 = count_bed_lines(bed_file1)
len_file2 = count_bed_lines(bed_file2)
len_intersection = count_bed_lines(intersection_file)
total_regions = len_file1 + len_file2 - len_intersection
k_successes = len_file1
n_draws = len_file2
p_value = 1 - hypergeom.cdf(len_intersection - 1, total_regions, k_successes, n_draws)
print(f"Total regions in file1: {len_file1}")
print(f"Total regions in file2: {len_file2}")
print(f"Total intersected regions: {len_intersection}")
print(f"Hypergeometric P-value: {p_value}")