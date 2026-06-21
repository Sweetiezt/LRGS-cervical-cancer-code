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
library(knitr)
library(rmarkdown)
library(IMvigor210CoreBiologies)
data(cds)
expMatrix <- counts(cds)
eff_length2 <- fData(cds)[,c("entrez_id","length","symbol")]
rownames(eff_length2) <- eff_length2$entrez_id
head(eff_length2)
feature_ids <- rownames(expMatrix)
expMatrix <- expMatrix[feature_ids %in% rownames(eff_length2),]
mm <- match(rownames(expMatrix),rownames(eff_length2))
eff_length2 <- eff_length2[mm,]

x <- expMatrix/eff_length2$length
eset <- t(t(x)/colSums(x))*1e6
summary(duplicated(rownames(eset)))

eset <- IOBR::anno_eset(eset = eset,annotation = eff_length2,symbol = "symbol",probe = "entrez_id",method = "mean")
tumor_type <- "blca"
if(max(eset)>100) eset <- log2(eset+1)

m2 <- eset
#rownames(m2) <- geneNames[rownames(m2)] # 注释
m2[1:4,1:4]

data(human_gene_signatures)
ind_genes <- human_gene_signatures
goi <- names(ind_genes)
goi # 20个基因集
for (sig in goi) {
  pData(cds)[, sig] <- NA # 在pData新增一列列名为sig，内容NA
  genes <- ind_genes[[sig]] # 在基因集的基因
  genes <- genes[genes %in% rownames(m2)] # 与原矩阵基因求交集
  tmp <- m2[genes, , drop=FALSE] # 交集基因的表达量
  pData(cds)[, sig] <- gsScore(tmp) # gsScore计算分数，并写入该对应列中
}

irf <- "Best Confirmed Overall Response"
feat <- "IC Level" # IC：免疫染色评分
table(pData(cds)[, irf])

tmpDat <- pData(cds)[ !is.na(pData(cds)[, feat]) & pData(cds)[, irf] != "NE", ]
table(tmpDat[, irf])
tmpDat[, irf] <- droplevels(tmpDat[, irf]) # 将level=NE删除掉
table(tmpDat[, irf])
ic <- table(tmpDat[, feat], tmpDat[, "binaryResponse"])
ic
pval <- signif(fisher.test(ic)$p.value, 2) # fisher.test检验，p值取两位有效数字
print(paste("Fisher P for IC by binary response:", pval))

ph<-tmpDat[,c(21,22,1:5,25:45)]
colnames(ph)[1]<-"futime"
colnames(ph)[2]<-"fustat"

m3<-as.data.frame(t(m2))

library(survival)
library(limma)
library(stringr)
library(openxlsx)
library(survival)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(pheatmap)
library(ggcorrplot)
library(survminer)
library(dplyr)
rt=read.table("riskcut2.txt",header=T,sep="\t",check.names=F,row.names=1)
rt1<-rt[,1:7]
rt2<-ph
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt1)
# multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
outTab=gsub("`","",outTab)

a = predict(multiCox,type="lp",rt)
b = apply(rt[,names(multiCox$coefficients)], 1,function(k){ sum(multiCox$coefficients * k)})
dat = data.frame(a,b,k = a-b)
riskScore=apply(rt1[,names(multiCox$coefficients)], 1,function(k){ sum(multiCox$coefficients * k)})
coxGene=rownames(multiCoxSum$coefficients)
gene<-row.names(outTab)


dat<-m3[,gene]
dat$ID<-rownames(dat)
ph$ID<-rownames(ph)
n<-merge(ph,dat,by="ID")
row.names(n)<-n[,1]
n<-n[,-1]
n1<-as.data.frame(t(n))
colnames(n)[8]<-"subtype"

#save.image("myzl.Rdata")

#计算出每个患者的风险评分，并按cut-off分为高低危
risk_scores=apply(n[,names(multiCox$coefficients)], 1,function(k){ sum(multiCox$coefficients * k)})
#risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
re<-cbind(id=rownames(cbind(n[,outCol],risk_scores)),cbind(n[,outCol],risk_scores))
res.cut <-surv_cutpoint(re, time = "futime", event = "fustat",variables = c("risk_scores"))
summary(res.cut)
res.cat <- surv_categorize(res.cut)

fit <- survfit(Surv(futime, fustat) ~risk_scores, data = res.cat)
ggsurvplot(fit, 
           data = res.cat, 
           risk.table = TRUE, 
           conf.int = F,
           pval = T)

cutoff<-res.cut[["cutpoint"]]
risk_scores<-apply(n[,names(multiCox$coefficients)], 1,function(k){ sum(multiCox$coefficients * k)})
risk_levels<-as.vector(ifelse(risk_scores>cutoff$cutpoint,"high","low"))
table(risk_levels)

fit <- survfit(Surv(futime, fustat) ~risk_levels, data = n)
ggsurvplot(fit, 
           data = n, 
           risk.table = TRUE, 
           conf.int = F,
           pval = T)

outCol=colnames(n)

re1<-cbind(id=rownames(cbind(n[,outCol],risk_scores,risk_levels)),cbind(n[,outCol],risk_scores,risk_levels))
re2<-re1[,12:36]
write.table(re2,"ccrcc.tsv",sep="\t",quote = F)


colnames(re1)[9]<-"response"

fit <- survfit(Surv(futime, fustat) ~re1$binaryResponse, data = re1)
ggsurvplot(fit,
           data = re1,
           risk.table = TRUE,
           conf.int = TRUE,
           pval = T)



rm<-re1
library(rms)
library(survival)
library(openxlsx)
library(VIM)
library(pROC)
rm$res<-rm$response
rm$res<-as.character(rm$res)
rm$res[rm$res!="CR"]<-"non_CR"
fit1 <- glm(res=="CR"~ risk_levels ,data = rm,family = "binomial")
summary(fit1)
1/exp(coef(fit1))
pre <- predict(fit1,type='response')
plot.roc(rm$res,pre,
         main="ROC Curve", percent=TRUE,
         print.auc=TRUE,
         ci=TRUE, of="thresholds",
         thresholds="best",
         print.thres="best")
#--------画图-------
tmpDat=re1
feat<-"risk_levels"
for (rr in c("CR")) {
  #rr=c("CR")
  tmpDat$group <- ifelse(tmpDat[, "response"] == rr, "CR", "NCR")
  tmpDat$group <- factor(tmpDat$group)
  ic <- table(tmpDat[, feat], tmpDat$group)
  ic # 列就是R,NR，行IC0，IC1，IC2+
  print(rr)
  print("all")
  print(signif(fisher.test(ic)$p.value, 2))
}

ic <- table(tmpDat[, feat], tmpDat[, "response"])
ic

nSamples <- rowSums(ic)
nSamples
#IC0  IC1 IC2+ 
#  83  112  102

ic1 <- prop.table(t(ic),margin=2) # 按列求频率
oldMar <- par()$mar
par(mar=c(.5, 1, 2, 0.5))
# 四个数字分别表示，下、左、上、右四个方向的内外边距，数值愈大距离越远
# 设置颜色colors
data(color_palettes)
# font sizes
labCex <- 0.9
namesCex <- 0.9
legendCex <- 0.9
titleCex <- 1
axisCex <- 0.9
titleF <- 1

# 画图
a <- barplot(ic1,
             ylab="fraction of patients", # y轴标签
             cex.names=namesCex, # 设置条形标签（barlabels）的大小
             cex.axis=axisCex, # x轴刻度字体大虚小
             cex.lab=labCex, # 坐标轴标签文字大小
             legend.text=rownames(ic1), # 图例标签
             col=color_palettes$response_palette, # 颜色
             width=0.16, # 柱子宽度
             xlim=c(0,0.5), # x轴范围
             args.legend=list(bty="n", # 布局设置，n：图例不加边框
                              cex=legendCex, # 图例字体大小
                              x="topright") # 图例位置
) # x轴不显示

a
# 0.112 0.304 0.496 # 对应三根柱子的横坐标
# x轴下添加标签
text(x = a, 
     y = par("usr")[2] -0.03, # x,y 文字坐标
     labels = levels(tmpDat[, feat]), # 标签内容："IC0"  "IC1"  "IC2+"
     srt = -45, # 倾斜角度
     xpd = TRUE, # 允许文本越出绘图区域
     adj=0, # 偏倚率
     cex=namesCex) # 字体大小
# >nSamples
#IC0  IC1 IC2+ 
#  83  112  102
# 在柱子上方添加文本
mtext(nSamples,
      side=3, # 在上方添加字体：83  112  102
      at=a,
      line=0,
      cex=1)
mtext(paste("Immune response","\n", "PD-L1"),
      side=3, # 取1234 ，意味着放置文本的边下左上右
      at=a[1], # 第二个柱子上
      line=1, # 外移文本
      font=titleF, # 字体样式
      cex=titleCex) # 字体大小
par(mar=oldMar) # 恢复默认边距

dev.off()


compare_means(risk_scores ~response, data = rm,
              method = "t.test")

rmhigh<-rm[with(rm,(rm$risk_levels=="high")),]
p4<-ggplot(rmhigh,aes(x=res,y=risk_scores))
compare_means(risk_scores ~ res, data = rmhigh,
              method = "wilcox.test")
p4+geom_boxplot(aes(fill=res))+xlab('risk')+scale_fill_brewer('risk',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = 15)+
  theme_classic()

fit <- survfit(Surv(futime, fustat) ~res, data = rmhigh)
ggsurvplot(fit,
           data = rmhigh,
           risk.table = TRUE,
           conf.int = TRUE,
           pval = T)

rmlow<-rm[with(rm,(rm$risk_levels=="low")),]
p4<-ggplot(rmlow,aes(x=res,y=risk_scores))
compare_means(risk_scores ~ res, data = rmlow,
              method = "wilcox.test")
p4+geom_boxplot(aes(fill=res))+xlab('response')+scale_fill_brewer('response',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = 12)+
  theme_classic()

fit <- survfit(Surv(futime, fustat) ~rmlow$res, data = rmlow)
ggsurvplot(fit,
           data = rmlow,
           risk.table = TRUE,
           conf.int = TRUE,
           pval = T)



rmhigh<-rm[with(rm,(rm$risk_levels=="high")),]
p4<-ggplot(rmhigh,aes(x=res,y=risk_scores))
compare_means(risk_scores ~ res, data = rmhigh,
              method = "wilcox.test")
p4+geom_boxplot(aes(fill=res))+xlab('risk')+scale_fill_brewer('risk',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = 15)+
  theme_classic()

fit <- survfit(Surv(futime, fustat) ~res, data = rmhigh)
ggsurvplot(fit,
           data = rmhigh,
           risk.table = TRUE,
           conf.int = TRUE,
           pval = T)

rmcr<-rm[with(rm,(rm$res=="non_CR")),]
p4<-ggplot(rmcr,aes(x=risk_levels,y=risk_scores))
compare_means(risk_scores ~ risk_levels, data = rmcr,
              method = "wilcox.test")
p4+geom_boxplot(aes(fill=risk_levels))+xlab('response')+scale_fill_brewer('response',palette='Set1')+ 
  stat_compare_means(method = "wilcox.test",label = "p.format",label.y = 12)+
  theme_classic()

fit <- survfit(Surv(futime, fustat) ~risk_levels, data = rmcr)
ggsurvplot(fit,
           data = rmcr,
           risk.table = TRUE,
           conf.int = TRUE,
           pval = T)


table(rm$res,rm$risk_levels)

ggplot(rm,aes(x=response,y=risk_levels))
