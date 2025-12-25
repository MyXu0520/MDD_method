#!/bin/bash
input_dir="/fast/bio2/bed_deeptools"
output_dir="/fast/bio2/bw_files"
chrom_sizes="/fast/bio2/hg38.chrome.size"
log_file="${output_dir}/conversion.log"
mkdir -p "${output_dir}"
echo "===  $(date) ===" > "${log_file}"
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${log_file}"
}
export LC_ALL=C
process_file() {
    local bed_file="$1"
    local base_name=$(basename "${bed_file}" .bed)
    log ": ${base_name}"
    sorted_bed="${output_dir}/${base_name}_sorted.bed"
    if [[ -f "${sorted_bed}" ]]; then
        log "，"
    else
        log "BED"
        if ! sort -k1,1 -k2,2n -k3,3n "${bed_file}" > "${sorted_bed}" 2>> "${log_file}"; then
            log "：${bed_file} ，"
            return 1
        fi
    fi
    bedgraph="${output_dir}/${base_name}.bedgraph"
    if [[ -f "${bedgraph}" ]]; then
        log "bedGraph，"
    else
        log "bedGraph"
        if ! bedtools genomecov -i "${sorted_bed}" -g "${chrom_sizes}" -bg > "${bedgraph}" 2>> "${log_file}"; then
            log "：${base_name} bedGraph，"
            rm -f "${sorted_bed}"
            return 1
        fi
    fi
    log "bedGraph"
    if ! sort -k1,1 -k2,2n "${bedgraph}" -o "${bedgraph}.sorted" 2>> "${log_file}"; then
        log "：bedGraph，"
        rm -f "${sorted_bed}" "${bedgraph}"
        return 1
    fi
    mv "${bedgraph}.sorted" "${bedgraph}"
    bw_file="${output_dir}/${base_name}.bw"
    if [[ -f "${bw_file}" ]]; then
        log "BigWig，"
    else
        log "BigWig"
        if ! bedGraphToBigWig "${bedgraph}" "${chrom_sizes}" "${bw_file}" 2>> "${log_file}"; then
            log "：${base_name} BigWig，"
            rm -f "${sorted_bed}" "${bedgraph}"
            return 1
        fi
    fi
    log ""
    rm -f "${sorted_bed}" "${bedgraph}"
    log ": ${base_name}\n------------------------------------"
    return 0
}
find "${input_dir}" -type f -iname "*.bed" -print0 | while IFS= read -r -d $'\0' bed_file; do
    process_file "$bed_file"
done
log "！：${output_dir}"
input_folder="/fast/bio2/bw_files"
regions_file="/fast/bio2/ESCA_sharedMDDs.bed"
output_matrix="/fast/bio2/result/matrix.gz"
bin_size=1000
before_length=50000
after_length=50000
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
plotHeatmap -m /fast/bio2/result/matrix.gz \
  -o raw_heatmap.pdf \
  --colorMap RdBu \
  --zMin -3 --zMax 3 \
  --kmeans 3