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
find "${input_dir}" -maxdepth 1 -name "*.bed" -print0 | while IFS= read -r -d $'\0' bed_file; do
    base_name=$(basename "${bed_file}" .bed)
    log ": ${base_name}"
    sorted_bed="${output_dir}/${base_name}_sorted.bed"
    if [[ -f "${sorted_bed}" ]]; then
        log "，"
    else
        log "BED"
        sort -k1,1 -k2,2n -k3,3n "${bed_file}" > "${sorted_bed}" 2>> "${log_file}" || exit 1
    fi
    bedgraph="${output_dir}/${base_name}.bedgraph"
    if [[ -f "${bedgraph}" ]]; then
        log "bedGraph，"
    else
        log "bedGraph"
        bedtools genomecov -i "${sorted_bed}" -g "${chrom_sizes}" -bg > "${bedgraph}" 2>> "${log_file}" || exit 1
    fi
    log "bedGraph"
    sort -k1,1 -k2,2n "${bedgraph}" -o "${bedgraph}.sorted"
    mv "${bedgraph}.sorted" "${bedgraph}"
    bw_file="${output_dir}/${base_name}.bw"
    if [[ -f "${bw_file}" ]]; then
        log "BigWig，"
    else
        log "BigWig"
        bedGraphToBigWig "${bedgraph}" "${chrom_sizes}" "${bw_file}" 2>> "${log_file}" || exit 1
    fi
    log ""
    rm -f "${sorted_bed}" "${bedgraph}"
    log ": ${base_name}\n------------------------------------"
done
log "！：${output_dir}"
plotHeatmap -m /fast/bio2/result/matrix.gz \
  -o raw_heatmap.pdf \
  --colorMap RdBu \
  --zMin -3 --zMax 3 \
  --kmeans 3