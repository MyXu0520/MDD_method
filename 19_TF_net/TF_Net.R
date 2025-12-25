library(Matrix)
library(motifmatchr)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(dplyr)
library(rtracklayer)
run_motif_enrichment <- function(dmr, pwmList, bg.ratio = 10, upstream = 2000, downstream = 200) {
  dmr_df <- dmr[,1:3]
  colnames(dmr_df) <- c("chrom", "start", "end")
  dmr_df$start <- as.numeric(as.character(dmr_df$start))
  dmr_df$end <- as.numeric(as.character(dmr_df$end))
  valid_chroms <- c(paste0("chr", 1:22))
  dmr_df <- dmr_df[dmr_df$chrom %in% valid_chroms, ]
  dmr.gr <- GRanges(
    seqnames = dmr_df$chrom,
    ranges = IRanges(start = dmr_df$start, end = dmr_df$end)
  )
  #seqlevels(dmr.gr, pruning.mode = "coarse") <- valid_chroms
  standard_chroms <- setdiff(standardChromosomes(BSgenome.Hsapiens.UCSC.hg38), c("chrM","chrX","chrY"))
  used_chroms <- intersect(standard_chroms, seqlevels(dmr.gr))
  dmr.gr <- keepSeqlevels(dmr.gr, used_chroms, pruning.mode = "coarse")
  seqlengths(dmr.gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[used_chroms]
  dmr.gr <- trim(dmr.gr)
  any(end(dmr.gr) > seqlengths(dmr.gr)[as.character(seqnames(dmr.gr))])
  promoter_file <- "D:\\CTC\\\\promoter.bed"
  cgi_file <- "D:\\CTC\\\\CGI.bed"
  enhancer_file <- "D:\\CTC\\\\Enhancer.bed"
  promoters.gr <- import.bed(promoter_file)
  cpg.gr <- import.bed(cgi_file)
  enhancer.gr <- import.bed(enhancer_file)
  seqlevels(promoters.gr, pruning.mode = "coarse") <- valid_chroms
  seqlevels(cpg.gr, pruning.mode = "coarse") <- valid_chroms
  seqlevels(enhancer.gr, pruning.mode = "coarse") <- valid_chroms
  seqnames(cpg.gr)
  promoters.reduced <- GenomicRanges::reduce(promoters.gr)
  cpg.reduced <- GenomicRanges::reduce(cpg.gr)
  enhancer.reduced <- GenomicRanges::reduce(enhancer.gr)
  bg.candidate <- c(promoters.reduced, cpg.reduced, enhancer.reduced)
  set.seed(123)
  dmr.lengths <- width(dmr.gr)
  bg.n <- length(dmr.lengths) * bg.ratio
  bg.gr <- GRanges()
  while (length(bg.gr) < bg.n) {
    sampled <- sample(bg.candidate, size = bg.n, replace = TRUE)
    start(sampled) <- start(sampled) + sample(-100:100, size = bg.n, replace = TRUE)
    width(sampled) <- sample(dmr.lengths, size = bg.n, replace = TRUE)
    bg.gr <- c(bg.gr, sampled)
  }
  bg.gr <- bg.gr[width(bg.gr) > 10]
  bg.gr <- bg.gr[seq_len(bg.n)]
  used_bg_chroms <- intersect(standard_chroms, seqlevels(bg.gr))
  bg.gr <- keepSeqlevels(bg.gr, used_bg_chroms, pruning.mode = "coarse")
  seqlengths(bg.gr) <- seqlengths(BSgenome.Hsapiens.UCSC.hg38)[used_bg_chroms]
  bg.gr <- trim(bg.gr)
  any(end(bg.gr) > seqlengths(bg.gr)[as.character(seqnames(bg.gr))])
  motif_ix_dmr <- matchMotifs(pwmList, dmr.gr, genome = BSgenome.Hsapiens.UCSC.hg38)
  motif_ix_bg <- matchMotifs(pwmList, bg.gr, genome = BSgenome.Hsapiens.UCSC.hg38)
  dmr_hits <- colSums(motifMatches(motif_ix_dmr))
  bg_hits <- colSums(motifMatches(motif_ix_bg))
  total_dmr <- length(dmr.gr)
  total_bg <- length(bg.gr)
  res <- data.frame(
    motif_id = names(dmr_hits),
    dmr_hits = dmr_hits,
    bg_hits = bg_hits,
    stringsAsFactors = FALSE
  )
  res <- res %>%
    mutate(
      fold_change = (dmr_hits / total_dmr) / (bg_hits / total_bg),
      p_value = mapply(function(h, b) {
        mat <- matrix(c(h, total_dmr - h, b, total_bg - b), nrow = 2)
        if (any(mat < 0) || any(!is.finite(mat))) return(1)
        fisher.test(mat, alternative = "greater")$p.value
      }, dmr_hits, bg_hits),
      FDR = p.adjust(p_value, method = "BH"),
      log10p = -log10(p_value + 1e-300),
      log10FDR = -log10(FDR + 1e-300),
      TF_name = sapply(names(pwmList), function(id) name(pwmList[[id]]))
    ) %>%
    arrange(desc(log10FDR))
  return(res)
}
library(dplyr)
library(data.table)
library(DSS)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
data <- read.csv("C:\\Users\\dell\\Desktop\\_MDD\\ESCA_sharedMDDs_covered_genes.csv", header = TRUE)
data <-data[,2:4]
colnames(data)<-c("chr","start","end")
data<-unique(data)
load("D:\\CTC\\_\\PWMatrixLists_5594.rds")
PWMatrixLists
dmr<-data
res_data <- run_motif_enrichment(dmr = dmr, pwmList = PWMatrixLists)
write.table(res_data,"C:\\Users\\dell\\Desktop\\_MDD\\motif.txt",quote = F)
res_filtered <- res_data[res_data$FDR < 0.01 & res_data$fold_change > 1.2 & res_data$dmr_hits >= 10, ]
motif_annotation_file <- "D:\\CTC\\_\\Human_TF_MotifList_v_1.01.csv"
motif_annotation <- read.csv(motif_annotation_file)
motif_annotation <-motif_annotation[,c(2,7)]
colnames(motif_annotation) <- c("TF_name", "motif_id")
TF_annotated <- merge(res_filtered, motif_annotation, by = "motif_id", all.x = TRUE)
TF_annotated$TF_name <- ifelse(TF_annotated$TF_name.x == "Unknown" & !is.na(TF_annotated$TF_name.y),
                               TF_annotated$TF_name.y,
                               TF_annotated$TF_name.x)
TF_annotated <- TF_annotated[, !(names(TF_annotated) %in% c("TF_name.x", "TF_name.y"))]
library(ggplot2)
library(dplyr)
TF<-unique(TF_annotated$TF_name)
tf<-c("E2F1","STAT1", "ETS1", "ETV4", "ELK1")
df<-TF_annotated[which(TF_annotated$TF_name%in%tf),]
df_max_fc <- df %>%
  group_by(TF_name) %>%
  summarise(max_fc = max(fold_change, na.rm = TRUE)) %>%
  top_n(5, max_fc)
df_max_fc_long <- reshape2::melt(df_max_fc, id.vars = "TF_name")
ggplot(df_max_fc_long, aes(x = reorder(TF_name, value), y = value)) +
  geom_bar(stat = "identity", fill = "#CDAA7D") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Top 5 TFs with Highest Fold Change", x = "TF Name", y = "Fold Change")
#ggsave("top_20_tfs.png", width = 10, height = 7)
library(dorothea)
data(dorothea_hs, package = "dorothea")
tf_targets <- dorothea_hs
tf_targets<-tf_targets[which(tf_targets$confidence=="A"|tf_targets$confidence=="B"),]
TF_names<-unique(df$TF_name)
targets <- subset(tf_targets, tf %in% TF_names)
table(targets$tf)
MDD_targets <- unique(targets$target)
library(clusterProfiler)
library(org.Hs.eg.db)
gene.df <- bitr(MDD_targets, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
ego <- enrichGO(
  gene         = gene.df$ENTREZID,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
ekegg <- enrichKEGG(
  gene         = gene.df$ENTREZID,
  organism     = 'hsa',
  pvalueCutoff = 0.05
)
####hypo
#go  c(4,12,66,67,97,111,129,237,348,717)
#kegg  c(2,4,5,6,16,38,47,79,152,154)
####hyper
#go  c(6,7,85,124,493,510,885,927,1349,1822)
#kegg  c(22,35,75,89,92,102,116,122,123,149)
ego1<-as.data.frame(ego)
rownames(ego1)<-1:nrow(ego1)
write.table(ego1,"C:\\Users\\dell\\Desktop\\_MDD\\ego.txt",quote = F,row.name=F)
ego1<-ego1[c(1,3,5,9,26,50,52,54,57,297),]
ekegg1<-as.data.frame(ekegg)
rownames(ekegg1)<-1:nrow(ekegg1)
write.table(ekegg1,"C:\\Users\\dell\\Desktop\\_MDD\\kegg.txt",quote = F,row.name=F)
ekegg1<-ekegg1[c(2,7,11,12,14,15,19,24,30,39),]
library(ggplot2)
#df_subset<-ego1
df_subset<-ekegg1
df_subset$GeneRatio <- sapply(df_subset$GeneRatio, function(x) {
  eval(parse(text = x))
})
df_subset$Description <- factor(df_subset$Description, levels = df_subset$Description[order(df_subset$GeneRatio)])
ggplot(df_subset, aes(x = GeneRatio, y = Description)) +
  geom_point(aes(size = Count, color = -log10(p.adjust))) +
  scale_color_gradient(low = "#0033FF", high = "#FF3300") +
  labs(
    x = "Gene Ratio",
    y = "KEGG",
    color = "-log10(p.adjust)",
    size = "Gene Count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.position = "right"
  )
library(dorothea)
data(dorothea_hs, package = "dorothea")
tf_targets <- dorothea_hs
tf_targets<-tf_targets[which(tf_targets$confidence=="A"|tf_targets$confidence=="B"),]
TF_names<-unique(df$TF_name)
targets <- subset(tf_targets, tf %in% TF_names)
table(targets$tf)
targets<-targets[,c(1,3)]
tf<-unique(targets$tf)  #TF
target<-unique(targets$target)
setwd("C:\\Users\\dell\\Desktop\\_MDD\\TF\\")
write.table(targets,"TF.txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(tf,"TF.txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(target,".txt",row.names =F,sep = "\t",quote = F,col.names =F)
oncogene<-read.table("C:\\Users\\dell\\Desktop\\_MDD\\TF\\Oncogene().txt")
onco_tf<-intersect(oncogene$V1,tf)
onco_targets<-intersect(oncogene$V1,target)
write.table(onco_tf,"TF.txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(onco_targets,".txt",row.names =F,sep = "\t",quote = F,col.names =F)
library(biomaRt)
df <- read.csv("C:/Users/dell/Desktop/_MDD/ESCA_sharedMDDs_covered_genes.csv")
head(df)
df$gene_ID_clean <- gsub("\\..*$", "", df$gene_ID)
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  filters = "ensembl_gene_id",
  values = unique(df$gene_ID_clean),
  mart = ensembl
)
df_annotated <- merge(df, gene_info, by.x = "gene_ID_clean", by.y = "ensembl_gene_id", all.x = TRUE)
head(df_annotated)
MDD_gene<-intersect(df_annotated$hgnc_symbol,target)
write.table(MDD_gene,"mdd.txt",row.names =F,sep = "\t",quote = F,col.names =F)
aa<-intersect(MDD_gene,oncogene$V1)
table(targets$tf)
#E2F1
E2F1_targets<-targets[which(targets$tf=="E2F1"),]
E2F1_targets<-E2F1_targets[which(E2F1_targets$target%in%df_annotated$hgnc_symbol),]
E2F1_tf<-unique(E2F1_targets$tf)  #TF
E2F1_target<-unique(E2F1_targets$target)
E2F1_onco<-intersect(oncogene$V1,E2F1_target)
setwd("C:\\Users\\dell\\Desktop\\_MDD\\TF\\TFMDD\\E2F1\\")
write.table(E2F1_targets,"TF.txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(E2F1_tf,"TF.txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(E2F1_target,".txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(E2F1_onco,".txt",row.names =F,sep = "\t",quote = F,col.names =F)
#ETS1
ETS1_targets<-targets[which(targets$tf=="ETS1"),]
ETS1_targets<-ETS1_targets[which(ETS1_targets$target%in%df_annotated$hgnc_symbol),]
ETS1_tf<-unique(ETS1_targets$tf)  #TF
ETS1_target<-unique(ETS1_targets$target)
ETS1_onco<-intersect(oncogene$V1,ETS1_target)
setwd("C:\\Users\\dell\\Desktop\\_MDD\\TF\\TFMDD\\ETS1\\")
write.table(ETS1_targets,"TF.txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(ETS1_tf,"TF.txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(ETS1_target,".txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(ETS1_onco,".txt",row.names =F,sep = "\t",quote = F,col.names =F)
#STAT1
STAT1_targets<-targets[which(targets$tf=="STAT1"),]
STAT1_targets<-STAT1_targets[which(STAT1_targets$target%in%df_annotated$hgnc_symbol),]
STAT1_tf<-unique(STAT1_targets$tf)  #TF
STAT1_target<-unique(STAT1_targets$target)
STAT1_onco<-intersect(oncogene$V1,STAT1_target)
setwd("C:\\Users\\dell\\Desktop\\_MDD\\TF\\TFMDD\\STAT1\\")
write.table(STAT1_targets,"TF.txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(STAT1_tf,"TF.txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(STAT1_target,".txt",row.names =F,sep = "\t",quote = F,col.names =F)
write.table(STAT1_onco,".txt",row.names =F,sep = "\t",quote = F,col.names =F)