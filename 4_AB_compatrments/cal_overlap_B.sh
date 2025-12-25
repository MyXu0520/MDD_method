#!/bin/bash
EAC_B_chr_all="/student1/MDD/result/1/ESCC_B_chr_all.bed"
ESCA_sharedMDDs="/student1/MDD/result/1/ESCA_sharedMDDs.bed"
output_file="/student1/MDD/result/1/EAC_B_chr_all_intersect.bed"
bedtools intersect -a $EAC_B_chr_all -b $ESCA_sharedMDDs -wa > $output_file
total_lines=$(wc -l < $EAC_B_chr_all)
intersect_lines=$(wc -l < $output_file)
overlap_rate=$(echo "scale=4; $intersect_lines / $total_lines" | bc)
echo "Overlap rate: $overlap_rate"
echo "Overlap rate: $overlap_rate" > /student1/MDD/result/1/ESCC_shared_overlap_rate.txt