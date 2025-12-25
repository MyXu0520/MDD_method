#!/bin/bash
input_folder="/fast/bio2/bw_files"
regions_file="/fast/bio2/ESCA_sharedMDDs.bed"
output_matrix="/fast/bio2/result/matrix.gz"
bin_size=50000
before_length=500000
after_length=500000
bw_files=$(ls ${input_folder}/*.bw)
computeMatrix reference-point \
  --regionsFileName $regions_file \
  --scoreFileName $bw_files \
  --outFileName $output_matrix \
  --referencePoint center \
  --beforeRegionStartLength $before_length \
  --afterRegionStartLength $after_length \
  --binSize $bin_size \
  --missingDataAsZero \
  --skipZeros
echo "ComputeMatrix Finished!"