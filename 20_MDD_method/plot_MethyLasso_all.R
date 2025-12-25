#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(pheatmap)
  library(GenomicRanges)
  library(VennDiagram)
  library(rtracklayer)
  library(RColorBrewer)
})

#########################################
# PARAMETERS AND CONFIGURATION
#########################################
option_list <- list(
  make_option("--dir", type="character", help="MethyLASSO output directory"),
  make_option("--prefix", type="character", help="Prefix of MethyLASSO output"),
  make_option("--samples", type="character", help="Comma-separated sample names")
)
opt <- parse_args(OptionParser(option_list = option_list))
OUT_DIR <- opt$dir
PREFIX <- opt$prefix
SAMPLES <- unlist(strsplit(opt$samples,","))

#########################################
# 1. READ PMD FILES
#########################################
read_pmd <- function(sam){
  file <- file.path(OUT_DIR, paste0(sam, "_PMDs.bed"))
  dt <- fread(file, header=FALSE)
  colnames(dt) <- c("chr","start","end","score")
  dt$length <- dt$end - dt$start
  dt$sample <- sam
  return(dt)
}

PMDs <- rbindlist(lapply(SAMPLES, read_pmd))

#########################################
# 2. PMD OVERLAP ANALYSIS (VENN DIAGRAM)
#########################################
# Example analysis for two sample groups (specific/shared)
adeno_samples <- grep("GEJ|Adeno", SAMPLES, value=TRUE)
squa_samples <- grep("ESCC|Squamous", SAMPLES, value=TRUE)

# Merge PMDs within each group
adeno_gr <- reduce(GRanges(seqnames=PMDs[PMDs$sample %in% adeno_samples]$chr,
                           IRanges(start=PMDs[PMDs$sample %in% adeno_samples]$start,
                                   end=PMDs[PMDs$sample %in% adeno_samples]$end)))
squa_gr <- reduce(GRanges(seqnames=PMDs[PMDs$sample %in% squa_samples]$chr,
                           IRanges(start=PMDs[PMDs$sample %in% squa_samples]$start,
                                   end=PMDs[PMDs$sample %in% squa_samples]$end)))

# Generate Venn diagram
venn_file <- file.path(OUT_DIR, paste0(PREFIX,"_PMD_overlap_venn.pdf"))
pdf(venn_file, width=5, height=5)
draw.pairwise.venn(area1=length(adeno_gr), area2=length(squa_gr),
                   cross.area=length(intersect(adeno_gr, squa_gr)),
                   category=c("Adeno","Squamous"), fill=c("skyblue","salmon"))
dev.off()

#########################################
# 3. PMD HEATMAP (AVERAGE METHYLATION PER PMD)
#########################################
pmd_matrix <- dcast(PMDs, chr + start + end ~ sample, value.var="score", fun.aggregate=mean)
rownames(pmd_matrix) <- paste(pmd_matrix$chr, pmd_matrix$start, pmd_matrix$end, sep="_")
pmd_matrix <- as.matrix(pmd_matrix[,-c(1:3)])

heatmap_file <- file.path(OUT_DIR, paste0(PREFIX,"_PMD_heatmap_multi.pdf"))
pdf(heatmap_file, width=10, height=8)
pheatmap(pmd_matrix,
         cluster_rows=TRUE, cluster_cols=TRUE,
         color=colorRampPalette(brewer.pal(9,"Blues"))(100),
         show_rownames=FALSE, show_colnames=TRUE)
dev.off()

#########################################
# 4. PMD ANNOTATION (GENE OVERLAP ANALYSIS)
#########################################
# Requires gene annotation file (e.g., hg38 GTF)
gtf_file <- "/data3/xumy_PMD/annotation/hg38.gtf"
genes <- import(gtf_file)
genes <- genes[genes$type=="gene"]

PMD_gr <- GRanges(seqnames=PMDs$chr, ranges=IRanges(start=PMDs$start,end=PMDs$end))
overlap_genes <- findOverlaps(PMD_gr, genes)
PMD_gene_count <- as.data.table(table(queryHits(overlap_genes)))
fwrite(PMD_gene_count, file.path(OUT_DIR,paste0(PREFIX,"_PMD_gene_overlap_count.tsv")), sep="\t")

#########################################
# 5. ENHANCED DMR VISUALIZATION (VOLCANO PLOT + BOXPLOT)
#########################################
dmr_file <- file.path(OUT_DIR, paste0(PREFIX, "_DMRs.bed"))
dmr <- fread(dmr_file, header=FALSE)
colnames(dmr) <- c("chr","start","end","lasso_weight")
dmr$neglogP <- -log10(abs(dmr$lasso_weight)+1e-6)

# Volcano plot visualization
ggplot(dmr, aes(x=lasso_weight, y=neglogP)) +
  geom_point(alpha=0.5) + theme_bw() +
  labs(title="DMR Volcano (LASSO Weight)", x="LASSO Weight", y="-log10(abs(weight))") +
  ggsave(file.path(OUT_DIR,paste0(PREFIX,"_DMR_volcano_extended.pdf")), width=6, height=5)

cat("[INFO] Extended visualization module completed. Outputs saved to:", OUT_DIR, "\n")