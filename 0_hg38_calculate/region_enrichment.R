#rm(list = ls())
options(stringAsFactors = F)
library(GenomicFeatures)
dat <- data.table::fread(file = "D:/MDD_result/1/hg38_rmblack.bed", header = F, sep = "\t", fill = T, skip = "", data.table=F)
dat[1:4, ]
bin <- GRanges(seqnames = dat[,1],
               ranges = IRanges(start = dat[,2]+1,
                                end = dat[,3]))
bin
saveRDS(bin,file = "D:/MDD_result/1/windowsNoBlack.rds")
# #test region
rm(list = ls())
options(stringAsFactors = F)
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("D:/EdgeDownload/gencode.v47.annotation (2).gtf.gz")
# saveDb(txdb, file=file)
tss <- promoters(txdb, upstream=200, downstream=200, columns=c("TXNAME", "GENEID"))
tss
bin_tss <- bin[bin %over% tss]
# saveRDS(bin_tss, file = "tss_bin.rds")
enrichment <- function(interest_regions,
                       Total_number,
                       interest_Terms,
                       anno_interest_Terms,
                       filter = F,
                       log2fc = 2,
                       fdr = 1e-3) {
  # All DEGs
  n <- length(interest_regions)
  # All genes
  N <- Total_number
  # k and M
  res_l <- sapply(interest_Terms, function(x){
    # M: All genes in TermA
    M <- length(x)
    # k: DEGs in TermA
    k <- sum(interest_regions %over% x)
    # enrichment P: phyper(k-1,M,N-M,n,lower.tail = F)
    p <- phyper(k-1,M,N-M,n,lower.tail = F)
    return(c(p,N,M,n,k))
  })
  res_l <- as.data.frame(t(res_l))
  colnames(res_l) <- c("P_hyper","All_genes_N","All_genes_in_term_M","All_DEGs_n","DEGs_in_Term_k")
  # P adjust
  FDR = p.adjust(res_l$P,method = "bonferroni")
  # enrichment fold change: (k/n)/(N/M)=kN/nM
  term_in_DEG <- res_l$DEGs_in_Term_k / res_l$All_DEGs_n
  term_in_BG <- res_l$All_genes_in_term_M / res_l$All_genes_N
  fc <- term_in_DEG / term_in_BG
  fc_res <- data.frame(Log2FC = log2(fc),
                       term_in_DEG = term_in_DEG,
                       term_in_BG = term_in_BG)
  # annotation
  anno <- cbind(fc_res,FDR,res_l,anno_interest_Terms)
  if (filter) {
    anno <- anno[anno$Log2FC>=log2fc & anno$FDR<fdr,]
  }
  return(anno)
}
anno <- enrichment(interest_regions, Total_number, interest_Terms, anno_interest_Terms, filter=F)
anno