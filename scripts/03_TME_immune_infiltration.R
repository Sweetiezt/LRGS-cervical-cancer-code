rm(list=ls())
dir='F:\\liulu\\TCGA\\'
od=getwd()
setwd(dir)
library(limma)
library(stringr)
library(rms)
library(survival)
library(openxlsx)
library(VIM)
rt3<-read.table("TCGAuse.tsv",sep="\t",quote = "")
library(dplyr)
rt<-as.data.frame(rt3)
out=rbind(ID=colnames(rt),rt)
write.table(out,file="precitumor.txt",sep="\t",quote=F,col.names=F)
library(e1071)
source("CIBERSORT.R")
results=CIBERSORT("ref.txt", "precitumor.txt", perm=1000, QN=TRUE)
write.csv(results,"CIBERSORT_Output.csv")

library(ggplot2)
library(openxlsx)
library(ggpubr)
library(ggExtra)
library(survival)
library(survivalROC)
library(survminer)
library(IOBR)
expr_coad<-rt
dim(expr_coad)
tme_deconvolution_methods

# MCPcounter
im_mcpcounter <- deconvo_tme(eset = expr_coad,
                             method = "mcpcounter"
)
# im_mcpcounter<-as.data.frame(im_mcpcounter)
# rownames(im_mcpcounter)<-im_mcpcounter$ID
# im_mcpcounter<-im_mcpcounter[,-1]
write.table(im_mcpcounter,"immcpcounter.tsv",sep="\t",quote = F)
#AAA<-read.table("mcpcounter.tsv",sep="\t")
# a<-as.data.frame(t(im_mcpcounter))
# b<-cor(a,resultsmcp)

mm<-merge(n,im_mcpcounter,by="ID")
library(ggplot2)
library(ggpubr)
r<-mm
TME.cells <- colnames(r)[24:33]
TME.cells
TME.data<-r
#Cell componment boxplot
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(TME.data) == TME.cells[i])
  sub <- data.frame(CellType = TME.cells[i],
                    Type= TME.data$risk2,
                    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}
#箱式图boxplot
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Immune cell composition",
  main = ""
) +
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

#tiff(file="免疫细胞.tiff",width = 20,height =15,units ="cm",compression="lzw",bg="white",res=400)
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Type",
  xlab = "",
  ylab = "Immune cell composition",
  main = "",palette =c("#FC4E07","#00AFBB" )
) +
  stat_compare_means(aes(group = Type),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
dev.off()
w1<-compare_means( Composition~Type , data = plot.info,   group.by = "CellType")
write.table(w1,"mianyixibaomcp-risk2.tsv",sep="\t")


# EPIC
im_epic <- deconvo_tme(eset = expr_coad,
                       method = "epic",
                       arrays = F
)
write.table(im_epic,"imepic.tsv",sep="\t",quote = F)

mepic<-merge(n,im_epic,by="ID")
library(ggplot2)
library(ggpubr)
r<-mepic
TME.cells <- colnames(r)[24:31]
TME.cells
TME.data<-r
#Cell componment boxplot
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(TME.data) == TME.cells[i])
  sub <- data.frame(CellType = TME.cells[i],
                    Type= TME.data$risk,
                    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}
#箱式图boxplot
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Immune cell composition",
  main = ""
) +
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Type",
  xlab = "",
  ylab = "Immune cell composition",
  main = "",palette =c("#FC4E07","#00AFBB" )
) +
  stat_compare_means(aes(group = Type),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
dev.off()
w1<-compare_means( Composition~Type , data = plot.info,   group.by = "CellType")
write.table(w1,"mianyixibaoepic-risk.tsv",sep="\t")

# xCell
im_xcell <- deconvo_tme(eset = expr_coad,
                        method = "xcell",
                        arrays = F
)
write.table(im_xcell,"imxcell.tsv",sep="\t",quote = F)

mxcell<-merge(n,im_xcell,by="ID")
library(ggplot2)
library(ggpubr)
r<-mxcell
TME.cells <- colnames(r)[24:90]
TME.cells
TME.data<-r
#Cell componment boxplot
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(TME.data) == TME.cells[i])
  sub <- data.frame(CellType = TME.cells[i],
                    Type= TME.data$risk2,
                    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}
#箱式图boxplot
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Immune cell composition",
  main = ""
) +
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

#tiff(file="免疫细胞.tiff",width = 20,height =15,units ="cm",compression="lzw",bg="white",res=400)
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Type",
  xlab = "",
  ylab = "Immune cell composition",
  main = "",palette =c("#FC4E07","#00AFBB" )
) +
  stat_compare_means(aes(group = Type),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
dev.off()
w1<-compare_means( Composition~Type , data = plot.info,   group.by = "CellType")
a<-aggregate(x=plot.info$Composition, by=list(plot.info$Type,plot.info$CellType),mean)

write.table(w1,"mianyixibaoxcell-risk2.tsv",sep="\t")
write.table(a,"mianyixibaoxcell-a.tsv",sep="\t")

sigxcell<-w1[w1$p.signif!="ns",]$CellType
sigxcell<-append("ID",sigxcell)

xcellsig<-im_xcell[,sigxcell]

mxcell2<-merge(n,xcellsig,by="ID")
colnames(im_estimate)

p1<-ggplot(mxcell2,aes(x=risk2,y=Astrocytes_xCell))
p1+geom_boxplot(aes(fill=risk2))+xlab('risk')+scale_fill_brewer('risk',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = max(mxcell2$Astrocytes_xCell))+
  theme_classic()


p1<-ggplot(mxcell2,aes(x=risk2,y=mxcell2$HSC_xCell))
p1+geom_boxplot(aes(fill=risk2))+xlab('risk')+scale_fill_brewer('risk',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = max(mxcell2$HSC_xCell))+#ylim(0,10^(-12)) +
  theme_classic()


p1<-ggplot(mxcell2,aes(x=risk2,y=mxcell2$Macrophages_M2_xCell))
p1+geom_boxplot(aes(fill=risk2))+xlab('risk')+scale_fill_brewer('risk',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = max(mxcell2$Macrophages_M2_xCell))+#ylim(0,10^(-12)) +
  theme_classic()


p1<-ggplot(mxcell2,aes(x=risk2,y=mxcell2$`Memory_B-cells_xCell`))
p1+geom_boxplot(aes(fill=risk2))+xlab('risk')+scale_fill_brewer('risk',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = max(mxcell2$`Memory_B-cells_xCell`))+#ylim(0,10^(-23)) +
  theme_classic()

p1<-ggplot(mxcell2,aes(x=risk2,y=MEP_xCell))
p1+geom_boxplot(aes(fill=risk2))+xlab('risk')+scale_fill_brewer('risk',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = max(mxcell2$MEP_xCell))+#ylim(0,10^(-23)) +
  theme_classic()

p1<-ggplot(mxcell2,aes(x=risk2,y=MPP_xCell))
p1+geom_boxplot(aes(fill=risk2))+xlab('risk')+scale_fill_brewer('risk',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = max(mxcell2$MPP_xCell))+#ylim(0,10^(-23)) +
  theme_classic()

p1<-ggplot(mxcell2,aes(x=risk2,y=Smooth_muscle_xCell))
p1+geom_boxplot(aes(fill=risk2))+xlab('risk')+scale_fill_brewer('risk',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = max(mxcell2$Smooth_muscle_xCell))+#ylim(0,10^(-23)) +
  theme_classic()

p1<-ggplot(mxcell2,aes(x=risk2,y=Pericytes_xCell))
p1+geom_boxplot(aes(fill=risk2))+xlab('risk')+scale_fill_brewer('risk',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = max(mxcell2$Pericytes_xCell))+#ylim(0,10^(-10)) +
  theme_classic()




# CIBERSORT
im_cibersort <- deconvo_tme(eset = expr_coad,
                            method = "cibersort",
                            arrays = F,
                            perm = 1000
)
write.table(im_cibersort,"im_cibersort.tsv",sep="\t",quote = F)


mciber<-merge(n,im_cibersort,by="ID")
library(ggplot2)
library(ggpubr)
r<-mciber
TME.cells <- colnames(r)[24:45]
TME.cells
TME.data<-r
#Cell componment boxplot
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(TME.data) == TME.cells[i])
  sub <- data.frame(CellType = TME.cells[i],
                    Type= TME.data$risk2,
                    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}
#箱式图boxplot
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Immune cell composition",
  main = ""
) +
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

#tiff(file="免疫细胞.tiff",width = 20,height =15,units ="cm",compression="lzw",bg="white",res=400)
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Type",
  xlab = "",
  ylab = "Immune cell composition",
  main = "",palette =c("#FC4E07","#00AFBB" )
) +
  stat_compare_means(aes(group = Type),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
dev.off()
w1<-compare_means( Composition~Type , data = plot.info,   group.by = "CellType")
write.table(w1,"mianyixibaoepic-risk2.tsv",sep="\t")


# IPS
im_ips <- deconvo_tme(eset = expr_coad,
                      method = "ips",
                      plot = F
)
write.table(im_ips,"im_ips.tsv",sep="\t",quote = F)

mips<-merge(n,im_ips,by="ID")
library(ggplot2)
library(ggpubr)
r<-mips
TME.cells <- colnames(r)[24:29]
TME.cells
TME.data<-r
#Cell componment boxplot
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(TME.data) == TME.cells[i])
  sub <- data.frame(CellType = TME.cells[i],
                    Type= TME.data$risk2,
                    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}
#箱式图boxplot
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Immune cell composition",
  main = ""
) +
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

#tiff(file="免疫细胞.tiff",width = 20,height =15,units ="cm",compression="lzw",bg="white",res=400)
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Type",
  xlab = "",
  ylab = "Immune cell composition",
  main = "",palette =c("#FC4E07","#00AFBB" )
) +
  stat_compare_means(aes(group = Type),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
dev.off()
w1<-compare_means( Composition~Type , data = plot.info,   group.by = "CellType")
write.table(w1,"mianyixibaoepic-risk2.tsv",sep="\t")


# quanTIseq
im_quantiseq <- deconvo_tme(eset = expr_coad,
                            method = "quantiseq",
                            scale_mrna = T
)
write.table(im_quantiseq,"im_quantiseq.tsv",sep="\t",quote = F)


mquantiseq<-merge(n,im_quantiseq,by="ID")
library(ggplot2)
library(ggpubr)
r<-mquantiseq
TME.cells <- colnames(r)[24:34]
TME.cells
TME.data<-r
#Cell componment boxplot
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(TME.data) == TME.cells[i])
  sub <- data.frame(CellType = TME.cells[i],
                    Type= TME.data$risk2,
                    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}
#箱式图boxplot
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Immune cell composition",
  main = ""
) +
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

#tiff(file="免疫细胞.tiff",width = 20,height =15,units ="cm",compression="lzw",bg="white",res=400)
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Type",
  xlab = "",
  ylab = "Immune cell composition",
  main = "",palette =c("#FC4E07","#00AFBB" )
) +
  stat_compare_means(aes(group = Type),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
dev.off()
w1<-compare_means( Composition~Type , data = plot.info,   group.by = "CellType")
write.table(w1,"mianyixibaoquantiseq-risk2.tsv",sep="\t")


# ESTIMATE
im_estimate <- deconvo_tme(eset = expr_coad,
                           method = "estimate"
)

im_estimate$ID<-gsub("-",".",im_estimate$ID)
write.table(im_estimate,"im_estimate.tsv",sep="\t",quote = F)


mestimate<-merge(n,im_estimate,by="ID")
library(ggplot2)
library(ggpubr)
r<-mestimate
TME.cells <- colnames(r)[24:27]
TME.cells
TME.data<-r
#Cell componment boxplot
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(TME.data) == TME.cells[i])
  sub <- data.frame(CellType = TME.cells[i],
                    Type= TME.data$risk2,
                    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}
#箱式图boxplot
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Immune cell composition",
  main = ""
) +
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

#tiff(file="免疫细胞.tiff",width = 20,height =15,units ="cm",compression="lzw",bg="white",res=400)
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Type",
  xlab = "",
  ylab = "Immune cell composition",
  main = "",palette =c("#FC4E07","#00AFBB" )
) +
  stat_compare_means(aes(group = Type),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
dev.off()
w1<-compare_means( Composition~Type , data = plot.info,   group.by = "CellType")
write.table(w1,"mianyixibaoestimate-risk2.tsv",sep="\t")


# TIMER
im_timer <- deconvo_tme(eset = expr_coad
                        ,method = "timer"
                        ,group_list = rep("coad",dim(expr_coad)[2])
)
write.table(im_timer,"imtimer.tsv",sep="\t",quote = F)

mtimer<-merge(n,im_timer,by="ID")
library(ggplot2)
library(ggpubr)
r<-mtimer
TME.cells <- colnames(r)[24:29]
TME.cells
TME.data<-r
#Cell componment boxplot
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(TME.data) == TME.cells[i])
  sub <- data.frame(CellType = TME.cells[i],
                    Type= TME.data$risk2,
                    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}
#箱式图boxplot
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Immune cell composition",
  main = ""
) +
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

#tiff(file="免疫细胞.tiff",width = 20,height =15,units ="cm",compression="lzw",bg="white",res=400)
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Type",
  xlab = "",
  ylab = "Immune cell composition",
  main = "",palette =c("#FC4E07","#00AFBB" )
) +
  stat_compare_means(aes(group = Type),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
dev.off()
w1<-compare_means( Composition~Type , data = plot.info,   group.by = "CellType")
write.table(w1,"mianyixibaotimer-risk2.tsv",sep="\t")


# 需要提供reference，暂不演示！
#im_svr <- deconvo_tme(eset = expr_coad
#                      ,method = "svr"
#                      ,arrays = F
#                      ,reference = 
#                      )
#im_lsei <- deconvo_tme(eset = expr_coad
#                       ,method = "lsei"
#                       ,arrays = F
#                       )
dim(im_cibersort)
im_cibersort[1:4,1:4]
dim(im_xcell)
im_xcell[1:4,1:4]

library(tidyr)
# 取前12个样本做演示
res<-cell_bar_plot(input = im_cibersort, title = "CIBERSORT Cell Fraction")

tme_combine <- im_mcpcounter %>% 
  inner_join(im_epic, by="ID") %>% 
  inner_join(im_xcell, by="ID") %>% 
  inner_join(im_cibersort, by="ID") %>% 
  inner_join(im_ips, by= "ID") %>% 
  inner_join(im_quantiseq, by="ID") %>% 
  inner_join(im_estimate, by= "ID") %>% 
  inner_join(im_timer, by= "ID")
tme_combine[1:4,1:4]
dim(tme_combine)
cellmarker<-list
im_ssgsea <- calculate_sig_score(eset = expr_coad
                                 , signature = cellmarker # 这个28种细胞的文件需要自己准备
                                 , method = "ssgsea" # 选这个就好了
)
im_ssgsea[1:4,1:4]
write.table(im_ssgsea,"imssgsea.tsv",sep="\t",quote = F)

mssgsea<-merge(n,im_ssgsea,by="ID")
library(ggplot2)
library(ggpubr)
r<-mssgsea
TME.cells <- colnames(r)[24:51]
TME.cells
TME.data<-r
#Cell componment boxplot
plot.info <- NULL
for (i in 1:length(TME.cells)) {
  idx.sub <- which(colnames(TME.data) == TME.cells[i])
  sub <- data.frame(CellType = TME.cells[i],
                    Type= TME.data$risk2,
                    Composition = TME.data[, idx.sub]
  )
  plot.info <- rbind(plot.info, sub)
}
#箱式图boxplot
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "CellType",
  xlab = "",
  ylab = "Immune cell composition",
  main = ""
) +
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))

#tiff(file="免疫细胞.tiff",width = 20,height =15,units ="cm",compression="lzw",bg="white",res=400)
ggboxplot(
  plot.info,
  x = "CellType",
  y = "Composition",
  color = "black",
  fill = "Type",
  xlab = "",
  ylab = "Immune cell composition",
  main = "",palette =c("#FC4E07","#00AFBB" )
) +
  stat_compare_means(aes(group = Type),
                     method = "wilcox.test",
                     label = "p.signif",
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "ns")))+
  
  theme(axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    vjust = 1
  ))
dev.off()
w1<-compare_means( Composition~Type , data = plot.info,   group.by = "CellType")
write.table(w1,"mianyixibaossgsea-risk2.tsv",sep="\t")


# 总的
signature_collection
# 代谢相关
signature_metabolism
# 微环境相关
signature_tme
# 肿瘤相关
signature_tme
tme_combine <- im_mcpcounter %>% 
  inner_join(im_epic, by="ID") %>% 
  inner_join(im_xcell, by="ID") %>% 
  inner_join(im_cibersort, by="ID") %>% 
  
  inner_join(im_quantiseq, by="ID") %>% 

  inner_join(im_timer, by= "ID")

tme_combine <- tme_combine %>% 
  inner_join(im_ssgsea, by = "ID")
write.table(tme_combine,"tmecombine.tsv",sep="\t",quote = F)

