shared_cover <- read.table("D:/MDD_result/17/shared_cover.bed",sep = "\t")
shared_cover <- shared_cover[,c(1,2,3,6)]
shared_MDD <- read.table("D:/MDD_result/17/ESCA_sharedMDDs.bed")
new_name_MDD <- paste0(shared_MDD$V1,"_",shared_MDD$V2,"_",shared_MDD$V3)
# new_name_cover <- paste0(shared_cover$V1,"_",shared_cover$V2,"_",shared_cover$V3)
rownames(shared_MDD) <- new_name_MDD
# rownames(shared_cover) <- new_name_cover
# distinct(df, column_name, .keep_all= TRUE)
shared_cover_filter <- shared_cover[!duplicated(shared_cover[,c(1,2,3)]),]
table(shared_cover_filter$V6)
# color <- c("#AA2D36","#2F56A3","#64329E")
# c("#334894","#81b1d0","#e4f2e9","#f6d386","#e3653f","#9c1e26")
c("#6CA3D4","#BF95C1","#61C1BF","#0074B3","#A5D395","#EDD283","#E1807E")
library(ggplot2)
library(tidyverse)
data <- data.frame(Shared_MDD = c(3177,3911,2714,3674,337,2951,2988), EAC_SPEC = c(1759,2521,11398,2147,80,1534,1583),ESCC_SPEC = c(259,426,190,340,9,238,233),
                   group = paste0("group",1:7)) %>%   pivot_longer(cols = !group, names_to = "X", values_to = "count")
data$group <- as.factor(data$group)
pdf("D:/MDD_result/17/Region_elements1.pdf")
p1 <- ggplot(data)+
  geom_bar(aes(X, count, fill = group), color = "#f3f4f4", alpha = 1.0,
           position = "fill", stat = "identity", size = 1)+
  scale_fill_manual(values = c("
                    labels=c("promoter","intron","intergenetic","exon","enhancer","5'UTR","3'UTR"))+
  # scale_color_manual(values = "gray")+
  # annotate("text", x = 2, y = 0.85, label="*", size = 5)+
  ggtitle("The location of MDD on the genome", subtitle = "(by type of MDD)")+
  scale_x_discrete(labels = c(expression(atop(bold("Shared_MDD"))),
                              expression(atop(bold("EAC_SPEC"))),
                              expression(atop(bold("ESCC_SPEC")))))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"),
        panel.background = element_rect(fill = "
  guides(fill=guide_legend(reverse=TRUE))+
  labs(fill="Region")
p1
dev.off()
data <- data.frame(core_MDD = c(265,371,237,330,12,224,254), core_CMI = c(195,195,189,195,11,195,193),
                   group = paste0("group",1:7)) %>%   pivot_longer(cols = !group, names_to = "X", values_to = "count")
data$group <- as.factor(data$group)
pdf("D:/MDD_result/17/Region_elements2.pdf")
p1 <- ggplot(data)+
  geom_bar(aes(X, count, fill = group), color = "#f3f4f4", alpha = 1.0,
           position = "fill", stat = "identity", size = 1)+
  scale_fill_manual(values = c("
                    labels=c("promoter","intron","intergenetic","exon","enhancer","5'UTR","3'UTR"))+
  # scale_color_manual(values = "gray")+
  # annotate("text", x = 2, y = 0.85, label="*", size = 5)+
  ggtitle("The location of MDD on the genome", subtitle = "(by type of Region)")+
  scale_x_discrete(labels = c(expression(atop(bold("core_MDD"))),
                              expression(atop(bold("core_CMI")))))+
  xlab("")+
  ylab("")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, face = "italic"),
        panel.background = element_rect(fill = "
  guides(fill=guide_legend(reverse=TRUE))+
  labs(fill="Region")
p1
 dev.off()
 data <- data.frame(ESCA_DMR = c(4156,5535,1547,4194,25,3036,1668), EAC_DMR = c(2360,1962,664,2338,17,1487,790),ESCC_SDMR = c(2105,3813,1664,1972,11,1160,762),
                    group = paste0("group",1:7)) %>%   pivot_longer(cols = !group, names_to = "X", values_to = "count")
 data$group <- as.factor(data$group)
 pdf("D:/MDD_result/17/Region_elements3.pdf")
 p1 <- ggplot(data)+
   geom_bar(aes(X, count, fill = group), color = "#f3f4f4", alpha = 1.0,
            position = "fill", stat = "identity", size = 1)+
   scale_fill_manual(values = c("
                     labels=c("promoter","intron","intergenetic","exon","enhancer","5'UTR","3'UTR"))+
   # scale_color_manual(values = "gray")+
   # annotate("text", x = 2, y = 0.85, label="*", size = 5)+
   ggtitle("The location of DMR on the genome", subtitle = "(by type of Region)")+
   scale_x_discrete(labels = c(expression(atop(bold("ESCA_DMR"))),
                               expression(atop(bold("EAC_DMR"))),
                               expression(atop(bold("ESCC_DMR")))))+
   xlab("")+
   ylab("")+
   theme_bw()+
   theme(panel.grid = element_blank(),
         plot.title = element_text(hjust = 0.5, face = "bold"),
         plot.subtitle = element_text(hjust = 0.5, face = "italic"),
         panel.background = element_rect(fill = "
   guides(fill=guide_legend(reverse=TRUE))+
   labs(fill="Region")
 p1
 dev.off()