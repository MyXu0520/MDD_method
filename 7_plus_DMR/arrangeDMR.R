load("D:/MDD_result/7/ESCC_dmrs_dml.RData")
dmrs <- dmrs[order(dmrs$chr,dmrs$start),]
ESCC_NT_DMR_bed <- dmrs[,c(1:3)]
write.table(ESCC_NT_DMR_bed,file = "ESCC_NT_DMR_bed.bed",sep = "\t",col.names = F,row.names = F,quote = F)
load("D:/MDD_result/7/EAC_dmrs_dml.RData")
dmrs <- dmrs[order(dmrs$chr,dmrs$start),]
EAC_NT_DMR_bed <- dmrs[,c(1:3)]
write.table(EAC_NT_DMR_bed,file = "EAC_NT_DMR_bed.bed",sep = "\t",col.names = F,row.names = F,quote = F)
load("D:/MDD_result/7/EAC_ESCC_dmrs_dml.RData")
dmrs <- dmrs[order(dmrs$chr,dmrs$start),]
EAC_ESCC_DMR_bed <- dmrs[,c(1:3)]
write.table(EAC_ESCC_DMR_bed,file = "EAC_ESCC_DMR_bed.bed",sep = "\t",col.names = F,row.names = F,quote = F)