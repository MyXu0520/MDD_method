#!/usr/bin/env bash

set -e  # Exit immediately on any error

# Configuration parameters
DATA_DIR="/data3/xumy_PMD/data"
OUT_DIR="/data3/xumy_PMD/MethyLasso_out"
mkdir -p $OUT_DIR

echo "[1] Automatically collecting sample lists..."

# Adenocarcinoma samples
ADENO_FILES=$(ls $DATA_DIR/AC*.bed $DATA_DIR/GEJ*.bed 2>/dev/null | tr '\n' ',' | sed 's/,$//')

# Squamous cell carcinoma samples
SQUAMOUS_FILES=$(ls $DATA_DIR/ESCC*.bed 2>/dev/null | tr '\n' ',' | sed 's/,$//')

echo "Adeno files: $ADENO_FILES"
echo "Squamous files: $SQUAMOUS_FILES"

echo ""
echo "[2] Running MethyLasso (Adenocarcinoma PMD detection)..."

Rscript MethyLasso.R \
    --n1 Adeno \
    --c1 $ADENO_FILES \
    --mC 4 --uC 5 \
    -o $OUT_DIR/Adeno \
    -s FALSE \
    -t 8

echo ""
echo "[3] Running MethyLasso (Squamous cell carcinoma PMD detection)..."

Rscript MethyLasso.R \
    --n1 Squamous \
    --c1 $SQUAMOUS_FILES \
    --mC 4 --uC 5 \
    -o $OUT_DIR/Squamous \
    -s FALSE \
    -t 8

echo ""
echo "[4] Extracting PMD BED files..."

# Output format: chr start end (PMD regions only)
awk 'BEGIN{OFS="\t"} $7=="PMD"{print $1,$2,$3}' $OUT_DIR/Adeno/Adeno_pmd.tsv     > $OUT_DIR/Adeno/Adeno.PMD.bed
awk 'BEGIN{OFS="\t"} $7=="PMD"{print $1,$2,$3}' $OUT_DIR/Squamous/Squamous_pmd.tsv > $OUT_DIR/Squamous/Squamous.PMD.bed

echo "[5] Identifying adenocarcinoma/squamous-specific and shared PMDs..."

# Adenocarcinoma-specific PMDs
bedtools subtract -a $OUT_DIR/Adeno/Adeno.PMD.bed -b $OUT_DIR/Squamous/Squamous.PMD.bed > $OUT_DIR/Adeno/Adeno_specific_PMD.bed

# Squamous cell carcinoma-specific PMDs
bedtools subtract -a $OUT_DIR/Squamous/Squamous.PMD.bed -b $OUT_DIR/Adeno/Adeno.PMD.bed > $OUT_DIR/Squamous/Squamous_specific_PMD.bed

# Shared PMDs between both cancer types
bedtools intersect -a $OUT_DIR/Adeno/Adeno.PMD.bed -b $OUT_DIR/Squamous/Squamous.PMD.bed > $OUT_DIR/shared_PMD.bed

echo "Analysis pipeline completed successfully!"