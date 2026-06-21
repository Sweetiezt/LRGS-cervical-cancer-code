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
rt2<-read.table("TCGATPM.txt",header = T)
library(dplyr)
exp<-rt2
library(limma)
exp= avereps(exp[,-c(1)],             # 多个探针对应一个基因，取均值
             ID = exp$gene_id)
exp=as.data.frame(exp)
data<-as.data.frame(t(exp))
group_list <-ifelse(as.numeric(str_sub(rownames(data),14,15)) <10,'tumor','normal')
cc<-cbind(data,group_list)
cc<-as.data.frame(cc)
exp<-cc[cc$group_list=="tumor",]
exp<-subset(exp,select=-c(group_list))
exp<-as.data.frame(t(exp))
exp<-log2(exp+1)
exp = normalizeBetweenArrays(exp) # 分位数标准化 method默认为quantile
exp=as.data.frame(exp)

p1 <- boxplot(exp,outline=FALSE,las=2,col = 'red',xaxt = 'n',ann = F)
p1
dev.off()
colSums(rt2[,-1])
TCGAsur<-read.table(file = 'sur.tsv',
                sep = '\t',
                header = T,
                stringsAsFactors = F,quote="")
TCGAsur<-TCGAsur[,c(1:3,9)]

library(readxl)
#lactate <- read_excel("lactate.xlsx")
lactate<-gene
gene3<-rownames(exp)
length(unique(gene3))
genec<-unique(gene)
common<-intersect(genec,gene3)
d<-setdiff(gene,common)
d
intersect("H2AZ1",gene3)
intersect("H1-5",gene3)
intersect("H4C3",gene3)
intersect("DYNC2I2",gene3)
intersect("POLR1G",gene3)
intersect("H1-2",gene3)
intersect("GFUS",gene3)
intersect("H1-3",gene3)

dd<-c(intersect("H2AZ1",gene3),
      intersect("H1-5",gene3),
      intersect("H4C3",gene3),
      intersect("DYNC2I2",gene3),
      intersect("POLR1G",gene3),
      intersect("H1-2",gene3),
      intersect("GFUS",gene3),
      intersect("H1-3",gene3))
common<-append(common,dd)
a3<-intersect(common,gene3)
rm<-exp[common,]
data<-na.omit(rm)
data<-t(data)
rownames(data)<-c(as.character(substr(rownames(data),1,16)))

TCGA<-as.data.frame(data)
TCGA$ID<-rownames(TCGA)
TCGA1<-merge(TCGAsur,TCGA,by="ID")  
row.names(TCGA1)<-TCGA1$ID
colnames(TCGA1)[3]<-"futime"
colnames(TCGA1)[2]<-"fustat"
r<-TCGA1
library(survival)
outTab=data.frame()
for(i in colnames(r[,5:ncol(r)])){
  cox <- coxph(Surv(futime, fustat) ~ r[,i], data = r)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab=rbind(outTab,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                     N=coxSummary$n
               )
  )
}


outTab = outTab[is.na(outTab$pvalue)==FALSE,]
outTab=outTab[order(as.numeric(as.vector(outTab$pvalue))),]
write.table(outTab,file="Lactate-related-uniCoxResultall.tsv",sep="\t",quote=F)

sigTab=outTab[as.numeric(as.vector(outTab$pvalue))<0.05,] 
write.table(sigTab,file="Lactate-related-uniCoxResult.tsv",sep="\t",quote=F)
sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
uniSigExp=r[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
write.table(uniSigExp,file="uniSigExpgenetcga.txt",sep="\t",row.names=F,quote=F)


outTab$pvalue<-as.numeric(outTab$pvalue)
outTab$HR<-as.numeric(outTab$HR)
yMax=max(-log10(outTab$pvalue))   
xMax=max(abs(outTab$HR))
library(ggplot2)

outTab$change <- ifelse(outTab$pvalue < 0.05,
                         ifelse(outTab$HR > 1,'Risk','Protective'),
                         'None')

outTab<-outTab1

ggplot(data= outTab, aes(y = -log10(pvalue), x = HR, color = change)) +
  geom_point(alpha=0.8, size = 2) +
  theme_bw(base_size = 15) +
  theme(plot.title=element_text(hjust=0.5),   #  标题居中
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) + # 网格线设置为空白
  geom_hline(yintercept= 0 ,linetype= 2 ) +
  scale_color_manual(name = "", 
                     values = c("#F44336", "#5C6BC0", "#BDC3C7"),
                     limits = c("Risk", "Protective", "None")) +
  ylim(0,yMax) + 
  xlim(0,2) +
  labs(title = 'Volcano', y = '-Log10(P.Value)', x = 'Hazard Ratio')
dev.off()



library("glmnet")
library("survival")
rt=read.table("uniSigExpgenetcga.txt",header=T,sep="\t",row.names=1,check.names=F)
rt$futime[rt$futime<=0]=1
rt$futime=rt$futime/1
set.seed(202507)
x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))

fit <- glmnet(x, y, family = "cox", maxit = 100000)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 100000,alpha = 1)
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
aa<-data.frame(lassoGene,actCoef)
write.table(aa,file="lassocoef.tsv",sep="\t",quote=F)
lassoGene=c("futime","fustat",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExptcga.txt",sep="\t",row.names=F,quote=F)


# 加载必要的包
library(randomForestSRC)
# 合并表达数据和生存数据
data_rfsrc <- read.table("uniSigExpgenetcga.txt",header=T,sep="\t",row.names=1,check.names=F)
exp_data<-data_rfsrc[,-c(1:2)]
# 构建随机森林模型
set.seed(123)
rf_fit <- rfsrc(Surv(futime, fustat) ~ ., 
                data = data_rfsrc,
                ntree = 1000,  # 树的数量
                mtry = floor(sqrt(ncol(exp_data))),  # 每次分裂考虑的变量数
                importance = TRUE)  # 计算变量重要性

# 获取基因重要性评分
gene_importance <- as.data.frame(rf_fit$importance)
gene_importance_df <- data.frame(Gene = rownames(gene_importance), 
                                 Importance = gene_importance$`rf_fit$importance`)
gene_importance_df <- gene_importance_df[order(gene_importance_df$Importance, decreasing = TRUE), ]

write.table(gene_importance_df,file="randomforest.tsv",sep="\t",quote=F)

# 可视化基因重要性
library(ggplot2)
top_genes <- gene_importance_df[1:10, ]  # 取前20个重要基因
ggplot(top_genes, aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  labs(x = "Gene", y = "Importance", title = "Top 10 Important Genes") +
  theme_minimal()

library(survival)
library(survivalROC)
library(survminer)
rt=read.table("lassoSigExptcga.txt",header=T,sep="\t",check.names=F,row.names=1)
rt=data_rfsrc[,ab]
rt<-rt[,1:7]
multiCox=coxph(Surv(futime, fustat) ~ ., data = rt)
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
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)

test_ph <- cox.zph(multiCox)
test_ph
plot(test_ph)
dev.off()

outTab<-as.data.frame(outTab)
outTab$HR<-as.numeric(outTab$HR)
outTab$pvalue<-as.numeric(outTab$pvalue)
library(tidyverse)
outTab$group <- ifelse( outTab$HR> 1 & outTab$pvalue < 0.05 ,"Risk" ,
                      ifelse(outTab$HR< 1 &  outTab$pvalue<0.05 ,"Protective","noSig")
)

outTab2<-outTab

out2 <- outTab %>% 
  mutate(hjust2 = ifelse(coef>0,1,0)) %>% 
  mutate(nudge_y = ifelse(coef>0,-0.1,0.1)) %>% 
  arrange(coef) %>% 
  rownames_to_column("ID") 

limt = max(abs(out2$coef))
ggplot(out2, aes(reorder(ID,coef), coef,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("Protective","Risk"), #设置颜色
                    values = c("navy","red3"))+
  geom_text(data = out2, aes(label =ID, #添加通路标签
                              y = nudge_y),
            nudge_x =0,nudge_y =0,hjust =out2$hjust,
            size = 6)+ #设置字体大小
  labs(x = "", #设置标题 和 坐标
       y=c("cox coefficient\n"),
       title = "")+
  scale_y_continuous(limits=c(-limt,limt))+
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank(), #主题微调
        panel.border = element_rect(size = 0.6),
        plot.title = element_text(hjust = 0.5,size = 18),
        axis.text.y = element_blank(),
        axis.title = element_text(hjust = 0.5,size = 18),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        #legend.position = "none" #去掉legend
  )


library(AnnotationHub)	#library导入需要使用的数据包
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(GOplot)
gene5<-colnames(rt)[3:7]

eg <- bitr(gene5, 
           fromType="SYMBOL", 
           toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
           OrgDb="org.Hs.eg.db")

head(eg)
go <- enrichGO(eg$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.05,
               keyType = 'ENTREZID',
               readable = T)

dim(go)
ego2<-data.frame(go)
write.table(ego2,"gogene5.tsv",sep = "\t")
library(stringr)

dotplot(go,showCategory=20)
dev.off()

# 加载必需包
library(ggplot2)
library(dplyr)
library(stringr)  # 用于简化GO条目名称（避免过长）

go_plot_data <-ego2 %>%
  # 1. 筛选显著条目（通常p.adjust < 0.05，可调整）
  filter(p.adjust < 0.05) %>%
  # 2. 按GO类别（BP/CC/MF）分组，每组取前10个最显著条目（避免图表拥挤）
  group_by(ONTOLOGY) %>%
  slice_min(p.adjust, n = 10,with_ties = FALSE) %>%  # 按p.adjust升序取前10（最显著）
  ungroup() %>%
  # 3. 简化GO条目名称（过长会导致标签重叠，保留前50个字符）
  #mutate(Description = str_trunc(Description, width = 50, side = "right")  # 截断长名称
  #) %>%
  # 4. 排序：按GO类别分组，再按p.adjust升序（确保显著条目在上方）
  arrange(ONTOLOGY, p.adjust) %>%
  # 5. 转换为因子（确保绘图时顺序正确，不自动按字母排序）
  mutate(
    Description = factor(Description, levels = unique(Description)),
    ONTOLOGY = factor(ONTOLOGY, levels = c("BP", "CC", "MF"))  # 固定分面顺序
  )

# 检查预处理后的数据
cat("预处理后GO数据预览：\n")
print(go_plot_data[, c("ONTOLOGY", "Description", "p.adjust", "Count")])


go_facet_plot <- ggplot(go_plot_data, 
                        aes(x = -log10(p.adjust),  # x轴：-log10(P值)，越大越显著
                            y = Description,       # y轴：GO条目描述
                            size = Count,       # 点大小：基因数越多，点越大
                            color = ONTOLOGY)) +    # 点颜色：按GO类别区分
  # 核心：散点（透明度过渡，避免重叠）
  geom_point(alpha = 0.8) +
  # 分面：按ONTOLOGY（BP/CC/MF）分面，每行一个类别
  facet_wrap(~ONTOLOGY, ncol = 1, scales = "free_y") +  # scales="free_y"：每组y轴独立（避免空条目）
  # 颜色配置（经典GO配色：BP-蓝色，CC-红色，MF-绿色）
  scale_color_manual(values = c("BP" = "#3498db", "CC" = "#e74c3c", "MF" = "#2ecc71")) +
  # 点大小配置（基因数越多，点越大，范围可调整）
  scale_size_continuous(range = c(2, 6)) +
  # 标签配置（清晰简洁）
  labs(
    x = expression(-log[10](Adjusted~P-value)),  # x轴标签：-log10(校正P值)
    y = "GO Term Description",                  # y轴标签：GO条目描述
    size = "Gene Count",                        # 图例：点大小对应基因数
    color = "GO Category",                      # 图例：点颜色对应GO类别
    title = "GO Enrichment Analysis Results (Faceted by Category)",  # 标题
    subtitle = "Significant Terms: Adjusted P-value < 0.05, Top 10 per Category"  # 副标题（说明筛选条件）
  ) +
  # 主题优化（避免标签重叠，提升可读性）
  theme_bw() +
  theme(
    # 分面标题：加粗、增大字体
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "gray90", color = NA),  # 分面背景色
    # 轴标签：调整大小，y轴标签左对齐（避免长名称截断）
    axis.text.y = element_text(size = 8, hjust = 1),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    # 图例：调整位置（右侧）和大小
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    # 标题：居中、加粗
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50"),
    # 调整边距（避免内容超出）
    plot.margin = margin(10, 20, 10, 20)
  ) +
  # 可选：添加P值阈值线（如p.adjust=0.05，对应-log10(0.05)=1.3）
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "gray50", size = 0.8)

# 显示分面图
print(go_facet_plot)

library(R.utils)
R.utils::setOption("clusterProfiler.download.method","auto")
kegg <- enrichKEGG(eg$ENTREZID, organism = "hsa", keyType = 'kegg', 
                   pvalueCutoff = 0.05,
                   pAdjustMethod = 'BH', 
                   minGSSize = 3,
                   maxGSSize = 500,
                   qvalueCutoff = 0.02,
                   use_internal_data = FALSE)
head(kegg)
ekegg<-as.data.frame(kegg)
write.table(ekegg,"kegggene5.tsv",sep="\t")
#可视化，和上面的一样

#tiff(file="DOTKEGG.tiff",width = 15,height =10,units ="cm",compression="lzw",bg="white",res=400)
dotplot(kegg,font.size=8, showCategory=20)	# 画气泡图
dev.off()

browseKEGG(ekk,'hsa04657')	# 显示通路图
p<-dotplot(kegg, showCategory=20) #气泡图
barplot(kegg,showCategory=20,drop=T) #柱状图

p

#riskScore=predict(multiCox,type="lp",newdata=rt)
#a = predict(multiCox,type="risk",rt)
b = apply(rt[,names(multiCox$coefficients)], 1,function(k){ sum(multiCox$coefficients * k)})
dat = data.frame(a,b,k = a-b)
#计算出每个患者的风险评分，并按cut-off分为高低危
riskScore=apply(rt[,names(multiCox$coefficients)], 1,function(k){ sum(multiCox$coefficients * k)})
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("futime","fustat",coxGene)
re<-cbind(id=rownames(cbind(rt[,outCol],riskScore)),cbind(rt[,outCol],riskScore))
risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))
fit <- survfit(Surv(futime, fustat) ~risk, data = re)
ggsurvplot(fit,
           data = re,
           risk.table = TRUE,
           conf.int = TRUE,
           pval = T)
write.table(cbind(id=rownames(cbind(rt[,outCol],riskScore,risk)),cbind(rt[,outCol],riskScore,risk)),
            file="riskcut.txt",
            sep="\t",
            quote=F,
            row.names=F)




library(survival)
library("survminer")
rt=read.table("riskcut2.txt",header=T,sep="\t",row.names = 1)
res.cut <-surv_cutpoint(rt, time = "futime", event = "fustat",variables = c("riskScore"))
summary(res.cut)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(futime, fustat) ~riskScore, data = res.cat)
ggsurvplot(fit,
           data = res.cat,
           risk.table = TRUE,
           conf.int = TRUE,
           pval = T)
risk2=res.cat$riskScore
table(risk2==rt$risk2)


diff=survdiff(Surv(futime, fustat) ~risk2,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(futime, fustat) ~ risk2, data = rt)
#tiff(file="survival.tiff",width = 20,height = 20,units ="cm",compression="lzw",bg="white",res=400)
ggsurvplot(fit, 
           data=rt,
           conf.int=TRUE,
           pval=paste0("p=",pValue),
           pval.size=4,
           risk.table=T,
           legend.labs=c("High risk", "Low risk"),
           legend.title="Risk",
           xlab="Time(Days)",
           break.time.by = 12*30,
           ggtheme = theme_light(),
           risk.table.y.text.col = T,
           risk.table.height = 0.18,
           risk.table.y.text = F,
           ncensor.plot = T,
           ncensor.plot.height = 0.18,
           conf.int.style = "ribbon")
dev.off()
summary(fit)

write.table(cbind(id=rownames(cbind(rt,risk2)),cbind(rt[,],risk2)),
            file="riskcut2.txt",
            sep="\t",
            quote=F,
            row.names=F)

rt=read.table("riskcut2.txt",header=T,sep="\t",row.names = 1)


par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =12*30, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#3年ROC
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =12*30*3, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#5年ROC
#pdf(file="ROC-5.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                predict.time =12*30*5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

roc1=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                 predict.time =12*30, method="KM")
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

roc2=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                 predict.time =12*30*3, method="KM")   #在此更改时间，单位为年
lines(roc2$FP,roc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="blue",lwd=2)

roc3=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$riskScore, 
                 predict.time =12*30*5, method="KM")   #在此更改时间，单位为年
lines(roc3$FP,roc3$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="green",lwd=2)

legend("bottomright", 
       c("1-year AUC:0.801","3-year AUC:0.765","5-year AUC:0.774"),
       lwd=2,lty = 1,
       col=c("red","blue","green"))
dev.off()

library(survcomp)

rt=read.table("riskcut2.txt",header=T,sep="\t",check.names=F,row.names=1)
cindex <- concordance.index(x=rt$riskScore,
                            surv.time = rt$futime, 
                            surv.event = rt$fustat,
                            method = "noether")
cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

######临床信息############
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

TCGAsur<-read.table(file = 'sur.tsv',
                    sep = '\t',
                    header = T,
                    stringsAsFactors = F,quote="")


rt=read.table("riskcut2.txt",header=T,sep="\t",check.names=F,row.names=1)
rt$ID<-rownames(rt)
sur1<-TCGAsur[,c(1:3,6,9,13,14,16,17,19:21,24,26,32,33)]
write.table(sur1,"TCGAclinical1.tsv",sep="\t",quote=F,row.names = F)###EXCEL

library(readxl)
sur<- read_excel("TCGAclinical2.xls")
sur$stage<-as.character(sur$stage)
sur$M<-as.character(sur$M)
sur$N<-as.character(sur$N)
sur$T<-as.character(sur$T)
sur$cancer_status<-as.character(sur$cancer_status)
sur$radiation<-as.character(sur$radiation)
sur$ethnicity<-as.character(sur$ethnicity)
sur$grade<-as.character(sur$grade)
sur$histology<-as.character(sur$histology)
sur2<-sur[,-c(2:3)]

rt$risk2[rt$risk2=="low"]<-"0"
rt$risk2[rt$risk2=="high"]<-"1"
rt1<-merge(rt,sur2,by="ID")

rt1$stage[rt1$stage=="2"]<-"1"

rt1$stage[rt1$stage=="3"]<-"4"
rt1$stage[rt1$stage=="4"]<-"2"

rt1$age <- as.numeric(rt1$age)
rt1$histology <- factor(rt1$histology)
rt1$stage<- factor(rt1$stage)
rt1$risk <- factor(rt1$risk)
rt1$risk2 <- factor(rt1$risk2)

outTab1=data.frame()
for(i in colnames(rt1[,4:ncol(rt1)])){
  cox <- coxph(Surv(futime, fustat) ~ rt1[,i], data = rt1)
  coxSummary = summary(cox)
  coxP=coxSummary$coefficients[,"Pr(>|z|)"]
  outTab1=rbind(outTab1,
               cbind(id=i,
                     z=coxSummary$coefficients[,"z"],
                     HR=coxSummary$conf.int[,"exp(coef)"],
                     HR.95L=coxSummary$conf.int[,"lower .95"],
                     HR.95H=coxSummary$conf.int[,"upper .95"],
                     pvalue=coxSummary$coefficients[,"Pr(>|z|)"],
                     N=coxSummary$n
               )
  )
}

outTab1 = data.frame()
for(i in colnames(rt1[,4:ncol(rt1)])){
  # 单因素Cox回归
  cox <- coxph(Surv(futime, fustat) ~ rt1[,i], data = rt1)
  coxSummary = summary(cox)
  
  # Schoenfeld 残差检验（PH假设）
  ph_test <- cox.zph(cox)
  ph_p <- ph_test$table[1, "p"] # 提取该基因的PH检验p值
  
  # 拼接结果
  outTab1 = rbind(outTab1,
                  cbind(id = i,
                        z = coxSummary$coefficients[,"z"],
                        HR = coxSummary$conf.int[,"exp(coef)"],
                        HR.95L = coxSummary$conf.int[,"lower .95"],
                        HR.95H = coxSummary$conf.int[,"upper .95"],
                        pvalue = coxSummary$coefficients[,"Pr(>|z|)"],
                        ph_pvalue = ph_p,  # 新增：PH检验p值
                        N = coxSummary$n
                  )
  )
}

# 查看最终表格
outTab1

outTab1 = outTab1[is.na(outTab1$pvalue)==FALSE,]
outTab1=outTab1[order(as.numeric(as.vector(outTab1$pvalue))),]

write.table(outTab1,file="uniclinical_CoxResult_revised.tsv",sep="\t",row.names=T,quote=F)
sigTab=outTab1[as.numeric(as.vector(outTab1$pvalue))<0.05,] #P值的筛选阈值可设为0.10
write.table(sigTab,file="uniCoxResult.Sig1_revised.tsv",sep="\t",row.names=F,quote=F)
sigGenes=c("futime","fustat")
sigGenes=c(sigGenes,as.vector(sigTab[,1]))
sigGenes<-unique(sigGenes)

uniSigExp=rt1[,sigGenes]
uniSigExp=cbind(id=row.names(uniSigExp),uniSigExp)
uniSigExp=cbind(id=rt1$ID,uniSigExp)
write.table(uniSigExp,file="TCGAuniSigExpgene11.tsv",sep="\t",row.names=T,quote=F)


aa<-rt1
library(tableone)  
library(survival)
##画森林图的包
library(forestplot)
library(stringr)
library(dplyr)
library(plyr)
y<- Surv(time=aa$futime,event=aa$fustat==1)#1为复发

#2.批量单因素回归模型建立：Uni_cox_model
Uni_cox_model<- function(x){
  FML <- as.formula(paste0 ("y~",x))
  cox<- coxph(FML,data=aa)
  cox1<-summary(cox)
  ph_test <- cox.zph(cox)
  ph_p <- ph_test$table[1, "p"] # 提取该基因的PH检验p值
  
  HR <- round(cox1$coefficients[,2],2)
  PValue <- round(cox1$coefficients[,5],3)
  CI5 <-round(cox1$conf.int[,3],2)
  CI95 <-round(cox1$conf.int[,4],2)
  ph_p<-round(ph_p,3)
  N<-cox1$n
  
  Uni_cox_model<- data.frame('Characteristics' = x,
                             'HR' = HR,
                             'CI5' = CI5,
                             'CI95' = CI95,
                             'p' = PValue,
                             "N"=N,
                             ph_pvalue = ph_p
                          )
  return(Uni_cox_model)}  

#3.将想要进行的单因素回归变量输入模型
#3-(1)查看变量的名字和序号
names(aa)
#3-(2)输入变量序号
variable.names<- colnames(aa)[c(4:21)] #例：这里选择了3-10号变量
variable.names
#4.输出结果
Uni_cox <- lapply(variable.names, Uni_cox_model)
Uni_cox<- ldply(Uni_cox,data.frame)

#5.优化表格，这里举例HR+95% CI+P 风格
Uni_cox$CI<-paste(Uni_cox$HR,"[", Uni_cox$CI5,'-',Uni_cox$CI95,"]")
#Uni_cox<-Uni_cox[,-3:-4]
#查看单因素cox表格
View(Uni_cox)
write.table(Uni_cox,"TCGAclinicalunicox-revised.tsv",sep="\t",quote=F)
result1 <-Uni_cox
a<-result1$Characteristics
unique(a)

b2<-c(  "RIBC2","HSD17B10","MIF","MSMO1","SPINT2","Riskscore","high", "Age","III/IV",  
 "With_tumor",  "Yes" , "Hispanic or latino" , "ACC","ASCC","BMI","LNE","Present" , "G2","G3","G4","LNP")   
result1$Characteristics<-b2

# result$Characteristics<-str_remove(result$Characteristics,"stage")
ins <- function(x) {c(x, rep(NA, ncol(result1)-1))}
result1<-result1[,c(1,2,3,4,8,5,6,7)]
##2-2：插入空行，形成一个新表
for(i in 5:7) {result1[, i] = as.character(result1[, i])}
result2<-rbind(c("Characteristics", NA, NA, NA, "HR(95%CI)","pvalue","N","ph_value"),
               result1[1:6, ],
              ins("Risk_levels"),
              ins("low"),
              result1[7, ],
              result1[8, ],
               ins("Stage"),
               ins("I/II"),  
               result1[9, ], 
               ins("Cancer_status"),
               ins("Tumor free"),
               result1[10, ],
               ins("Radiation"),
               ins("No"),
               result1[11, ],
               ins("Ethnicity"),
               ins("Not hispanic or latino"),
               result1[12, ],
               ins("Histology"),
               ins("SCC"),
               result1[13:16, ],
              ins("LVSI"),
              ins("Absent"),
               result1[17, ],
               ins("Grade"),
               ins("G1"),
               result1[18:21, ],
               c(NA, NA, NA, NA, NA,NA,NA,NA)#
)
for(i in 2:4) {result2[, i] = as.numeric(result2[, i])}

result2$CI95[result2$CI95=="Inf"]<-0
fig2 <- forestplot(result2[,c(1,5,6,7,8)], 
                   mean=result2[,2],   
                   lower=result2[,3],  
                   upper=result2[,4], 
                   zero=1,        
                   boxsize=0.1, 
                   col=fpColors(box='red',summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
                   
                   graph.pos=6)     


#tiff(file="foruni.tiff",width = 30,height = 25,units ="cm",compression="lzw",bg="white",res=400)
fig2
dev.off()


colnames(rt1)[10]<-"Risk_levels"
#1-1
myVars <-colnames(rt1)[4:21]
#1-2
catVars <-  c("Risk_levels", "stage" ,"cancer_status", "radiation", "ethnicity","histology" , "LVI", "grade")
#1-3
table1<- print(CreateTableOne(vars=myVars,
                              data = rt1,
                              factorVars = catVars),
               showAllLevels=TRUE)

#2. 在基线表table1里插入空行，使它的行数和变量跟result一致
N<-rbind(c(NA,NA),
         table1[2:7, ],
         c(NA,NA),
         table1[8:10,],
         c(NA,NA),
         table1[11:12,],
         c(NA,NA), 
         table1[13:14,],
         c(NA,NA), 
         table1[15:16,],
         c(NA,NA), 
         table1[17:18,],
         c(NA,NA), 
         table1[19:23,],
         c(NA,NA), 
         table1[24:25,],
         c(NA,NA), 
         table1[26:30,],
         c(NA,NA)
         
         )

#N<-N[,-1]
N<-data.frame(N)
#3.把N表和result表合在一起
result3<-cbind(result2,N)
#调顺序。变为:变量-N-HR......顺序
result3<-result3[,c(1,10,2:4,5,6,7,8)]

#4.优化第一行。第一行行名中加入"Number(%)"
for(i in 2:7) {result3[, i] = as.character(result3[, i])}
result3<-rbind(c("Characteristics","Number (%)",NA,NA,NA,"HR (95%CI)","P.value","N","ph_value"),
               result3[2:nrow(result3),])
for(i in 3:5) {result3[, i] = as.numeric(result3[, i])}


#画图fig-3，注：因为多了一列，所以要注意改代码数字
fig3<-forestplot(result3[,c(1,2,6,7,8,9)], 
                 mean=result3[,3],   
                 lower=result3[,4],
                 upper=result3[,5],  
                 zero=1,          
                 boxsize=0.4,  
                 graph.pos=3)

fig3






fig3_1<-forestplot(result3[,c(1,2,6,7,8,9)], 
                   mean=result3[,3],   
                   lower=result3[,4],  
                   upper=result3[,5], 
                   zero=1,            
                   boxsize=0.6,      
                   graph.pos= "right" ,
                   hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                   "2" = gpar(lty=2),
                                   "40"= gpar(lwd=2,lty=1,columns=c(1:4)) ),
                   graphwidth = unit(.25,"npc"),
                   xlab="",
                   xticks=c(0.4,1,3,5,7,10) ,
                   is.summary=c(T,T,T,T,T,T,T,T,F,F,T,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,F,T,T,T,F,F,T,F,F,F,F,T),
                   txt_gp=fpTxtGp(
                     label=gpar(cex=1),
                     ticks=gpar(cex=1), 
                     xlab=gpar(cex=1.5), 
                     title=gpar(cex=2)),
                   lwd.zero=1,
                   lwd.ci=1.5,
                   lwd.xaxis=2, 
                   lty.ci=1.5,
                   ci.vertices =T,
                   ci.vertices.height=0.2, 
                   clip=c(0.1,8),
                   ineheight=unit(8, 'mm'), 
                   line.margin=unit(8, 'mm'),
                   colgap=unit(6, 'mm'),
                   fn.ci_norm="fpDrawDiamondCI", 
                   col=fpColors(box ='#021eaa', 
                                lines ='#021eaa', 
                                zero = "black"))

fig3_1
dev.off()
result4<-result3[c(1:6),]

result4$Overall[1]<-"Median(SE)"
fig4_1<-forestplot(result4[,c(1,2,6,7,8,9)], 
                   mean=result4[,3],   
                   lower=result4[,4],  
                   upper=result4[,5], 
                   zero=1,            
                   boxsize=0.6,      
                   graph.pos= "right" ,
                   hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                   "2" = gpar(lty=2),
                                   "7"= gpar(lwd=2,lty=1,columns=c(1:4)) ),
                   graphwidth = unit(.25,"npc"),
                   xlab="",
                   xticks=c(0.4,1,3,5,7,10,20) ,
                   is.summary=c(T,T,T,T,T,T,T),
                   txt_gp=fpTxtGp(
                     label=gpar(cex=1),
                     ticks=gpar(cex=1), 
                     xlab=gpar(cex=1.5), 
                     title=gpar(cex=2)),
                   lwd.zero=1,
                   lwd.ci=1.5,
                   lwd.xaxis=2, 
                   lty.ci=1.5,
                   ci.vertices =T,
                   ci.vertices.height=0.2, 
                   clip=c(0.1,8),
                   ineheight=unit(8, 'mm'), 
                   line.margin=unit(8, 'mm'),
                   colgap=unit(6, 'mm'),
                   fn.ci_norm="fpDrawDiamondCI", 
                   col=fpColors(box ='#021eaa', 
                                lines ='#021eaa', 
                                zero = "black"))

fig4_1

result5<-result3[-c(2:6),]

fig5_1<-forestplot(result5[,c(1,2,6,7,8,9)], 
                   mean=result5[,3],   
                   lower=result5[,4],  
                   upper=result5[,5], 
                   zero=1,            
                   boxsize=0.6,      
                   graph.pos= "right" ,
                   hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                                   "2" = gpar(lty=2),
                                   "35"= gpar(lwd=2,lty=1,columns=c(1:4)) ),
                   graphwidth = unit(.25,"npc"),
                   xlab="",
                   xticks=c(0.4,1,3,5,7,10,20) ,
                   is.summary=c(T,T,T,F,F,T,T,F,F,T,F,F,T,F,F,T,F,F,T,F,F,F,T,T,T,F,F,T,F,F,F,F,T),
                   txt_gp=fpTxtGp(
                     label=gpar(cex=1),
                     ticks=gpar(cex=1), 
                     xlab=gpar(cex=1.5), 
                     title=gpar(cex=2)),
                   lwd.zero=1,
                   lwd.ci=1.5,
                   lwd.xaxis=2, 
                   lty.ci=1.5,
                   ci.vertices =T,
                   ci.vertices.height=0.2, 
                   clip=c(0.1,8),
                   ineheight=unit(8, 'mm'), 
                   line.margin=unit(8, 'mm'),
                   colgap=unit(6, 'mm'),
                   fn.ci_norm="fpDrawDiamondCI", 
                   col=fpColors(box ='#021eaa', 
                                lines ='#021eaa', 
                                zero = "black"))

fig5_1





ddist <- datadist(rt1)
options(datadist='ddist')
aa<-rt1
library(tableone)  
library(survival)
##画森林图的包
library(forestplot)
library(stringr)
library(dplyr)
library(plyr)
y<- Surv(time=aa$futime,event=aa$fustat==1)#1为复发
#2.批量单因素回归模型建立：Uni_cox_model
Uni_cox_model<- function(x){
  FML <- as.formula(paste0 ("y~",x))
  cox<- coxph(FML,data=aa)
  cox1<-summary(cox)
  HR <- round(cox1$coefficients[,2],2)
  PValue <- round(cox1$coefficients[,5],3)
  CI5 <-round(cox1$conf.int[,3],2)
  CI95 <-round(cox1$conf.int[,4],2)
  Uni_cox_model<- data.frame('Characteristics' = x,
                             'HR' = HR,
                             'CI5' = CI5,
                             'CI95' = CI95,
                             'p' = PValue)
  return(Uni_cox_model)}  
names(aa)
variable.names<- colnames(aa)[c(4:21)] #例：这里选择了3-10号变量
variable.names
Uni_cox <- lapply(variable.names, Uni_cox_model)
Uni_cox<- ldply(Uni_cox,data.frame)

#5.优化表格，这里举例HR+95% CI+P 风格
Uni_cox$CI<-paste(Uni_cox$CI5,'-',Uni_cox$CI95)
#Uni_cox<-Uni_cox[,-3:-4]

View(Uni_cox)
#1.提取单因素p<0.05变量
unique(Uni_cox$Characteristics[Uni_cox$p<0.05][c(7,8,10)])

# #2.多因素模型建立
mul_cox_model<- as.formula(paste0("y~",
                                   paste0(unique(Uni_cox$Characteristics[Uni_cox$p<0.05][c(7,8,10)]),
                                          collapse = "+")))
mul_cox<-coxph(mul_cox_model,data=aa)
cox4<-summary(mul_cox) 
cox4


#一-1.cox多因素回归分析
mul_cox<-coxph(mul_cox_model,data=aa)
#一-2 multi1：提取：变量+HR+95%CI+95%CI
mul_cox1 <- summary(mul_cox)
colnames(mul_cox1$conf.int)
multi1<-as.data.frame(round(mul_cox1$conf.int[, c(1, 3, 4)], 2))
#一-3、multi2：提取：HR(95%CI)和P
multi2<-ShowRegTable(mul_cox, 
                     exp=TRUE, 
                     digits=2, 
                     pDigits =3,
                     printToggle = TRUE, 
                     quote=FALSE, 
                     ciFun=confint)

ph_test <- cox.zph(mul_cox) # 多因素模型PH检验
ph_test 
plot(ph_test)

#一-4.将两次提取结果合并成表；取名result
result <-cbind(multi1,multi2);result
#一-5.行名转为表格第一列，并给予命名"Characteristics"
result<-tibble::rownames_to_column(result, var = "Characteristics");result
b3<-c(  "high","III/IV","BMI")   

result$Characteristics<-b3

fig1<- forestplot(result[,c(1,5,6)], #告诉函数，合成的表格result的第1，5，6列还是显示数字
                  mean=result[,2],   #告诉函数，表格第2列为HR，它要变成森林图的小方块
                  lower=result[,3],  #告诉函数表格第3列为5%CI，
                  upper=result[,4],  #表格第5列为95%CI，它俩要化作线段，穿过方块
                  zero=1,            #告诉函数，零线或参考线为HR=1即x轴的垂直线
                  boxsize=0.3,       #设置小黑块的大小
                  graph.pos=2)       #森林图应插在图形第2列


fig1
#2. 给参考变量插入空行
#2-1.这步代码不用改
ins <- function(x) {c(x, rep(NA, ncol(result)-1))}
##2-2：插入空行，形成一个新表
for(i in 5:6) {result[, i] = as.character(result[, i])}
result<-rbind(c("Characteristics", NA, NA, NA, "HR(95%CI)","pvalue"),
              ins("Risk_levels"),
              ins("low"),  
              result[1, ], 
              ins("Stage"),
              ins("I/II"),  
              result[2, ],
              result[3:nrow(result), ],
              c(NA, NA, NA, NA, NA,NA)#
)
for(i in 2:4) {result[, i] = as.numeric(result[, i])}



fig2 <- forestplot(result[,c(1,5,6)], 
                   mean=result[,2],   
                   lower=result[,3],  
                   upper=result[,4], 
                   zero=1,        
                   boxsize=0.5, 
                   col=fpColors(box='red',summary="#8B008B",lines = 'black',zero = '#7AC5CD'),
                   
                   graph.pos=4)     

fig2

#step1，筛选P<0.1者
f <- as.formula(Surv(futime,fustat) ~Risk_levels+stage )
cox <- coxph(f,data=rt1)
summary(cox)

ph_test <- cox.zph(cox) # 多因素模型PH检验
ph_test 
plot(ph_test)

ddist <- datadist(rt1)
options(datadist='ddist')

rt1$Risk_levels<-factor(rt1$Risk_levels,labels =c( "low","high"))
rt1$stage<-factor(rt1$stage,labels = c("I/II","III/IV"))

cox <- cph(Surv(futime,fustat) ~Risk_levels+stage,surv=T,x=T, y=T,data=rt1) 
surv <- Survival(cox)
sur_1_year<-function(x)surv(1*365,lp=x)
sur_3_year<-function(x)surv(1*365*3,lp=x)
sur_5_year<-function(x)surv(1*365*5,lp=x)
nom_sur <- nomogram(cox,fun=list(sur_1_year,sur_3_year,sur_5_year),lp= F,funlabel=c('1-Year survival','3-Year survival','5-Year survival'),maxscale=100,fun.at= c('1.0','0.9','0.8','0.7','0.6','0.5','0.4','0.3','0.2','0.1','0'))
#pdf("nomogram.pdf")
plot(nom_sur,xfrac=0.4)
dev.off()

library(survivalROC)
library(survival)

#计算Nomogram模型的risk score，并分为高低危组
cox_m <- coxph(Surv(futime,fustat) ~ Risk_levels+stage, data = rt1)
ph_test <- cox.zph(cox_m) # 多因素模型PH检验
ph_test 
plot(ph_test)


risknom<-predict(cox_m,type="lp",newdata=rt1)


re<-cbind(id=rownames(cbind(rt,risknom)),cbind(rt,risknom))
#risk=as.vector(ifelse(riskScore>median(riskScore),"high","low"))

res.cut <-surv_cutpoint(re, time = "futime", event = "fustat",variables = c("risknom"))
summary(res.cut)
res.cat <- surv_categorize(res.cut)
fit <- survfit(Surv(futime, fustat) ~risknom, data = res.cat)
ggsurvplot(fit, 
           data = res.cat, 
           risk.table = TRUE, 
           conf.int = TRUE,
           pval = T)

risknomlevel<-res.cat$risknom


write.table(cbind(id=rownames(cbind(rt1,risknom,risknomlevel)),cbind(rt1,risknom,risknomlevel)),"risknomlevel.txt",sep="\t",quote=F,row.names=F)

#绘制ROC曲线
rt=read.table("risknomlevel.txt",header=T,sep="\t",check.names=F,row.names=1)
#1年ROC
#pdf(file="ROC-1.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risknom, 
                predict.time =1*365, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#3年ROC
#pdf(file="ROC-3.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risknom, 
                predict.time =365*3, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#5年ROC
#pdf(file="ROC-5.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risknom, 
                predict.time =365*5, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

#整合1，3，5年ROC
#tiff(file="ROC.tiff",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)

roc1=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risknom, 
                 predict.time =365, method="KM")
plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="1-Specificity", ylab="Sensitivity",
     main=paste("ROC curve"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)

roc2=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risknom, 
                 predict.time =365*3, method="KM")   #在此更改时间，单位为年
lines(roc2$FP,roc2$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="blue",lwd=2)
#text(locator(1), paste("1 year",round(roc2$AUC,3),sep=":"),col="blue")

roc3=survivalROC(Stime=rt$futime, status=rt$fustat, marker = rt$risknom, 
                 predict.time =365*5, method="KM")   #在此更改时间，单位为年
lines(roc3$FP,roc3$TP,type="l",xlim=c(0,1),ylim=c(0,1),col="green",lwd=2)
#text(locator(1), paste("2 year",round(roc2$AUC,3),sep=":"),col="green")

legend("bottomright", 
       c("1-year AUC:0.844","3-year AUC:0.734","5-year AUC:0.715"),
       lwd=2,
       col=c("red","blue","green"))
dev.off()

library(survcomp)
rt=read.table("risknomlevel.txt",header=T,sep="\t",check.names=F,row.names=1)
rt<-na.omit(rt)
cindex <- concordance.index(x=rt$risknom,
                            surv.time = rt$futime, 
                            surv.event = rt$fustat,
                            method = "noether")
cindex$c.index
cindex$se
cindex$lower
cindex$upper
cindex$p.value

#1-year
cox1 <- cph(Surv(futime,fustat) ~ stage+Risk_levels,surv=T,x=T, y=T,time.inc = 1*365,data=rt) 
cal <- calibrate(cox1, cmethod="KM", method="boot", u=1*365, m= 10, B=1000)
#pdf("calibration1.pdf",12,8)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal,lwd=3,lty=2,errbar.col="black",xlim = c(0,1.0),ylim = c(0,1.0),xlab ="Nomogram-Predicted Probability of 1-Year Survival",ylab="Actual 1-Year Survival",col="blue")
#lines(cal,c('mean.predicted','KM'),type = 'l',lwd = 3,col ="black" ,pch = 16)
box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
dev.off()

#3-year
cox1 <- cph(Surv(futime,fustat) ~ stage+Risk_levels,surv=T,x=T, y=T,time.inc = 3*365,data=rt) 
cal <- calibrate(cox1, cmethod="KM", method="boot", u=3*165, m= 10, B=1000)
#pdf("calibration3.pdf",12,8)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal,lwd=3,lty=2,errbar.col="black",xlim = c(0,1.0),ylim = c(0,1.0),xlab ="Nomogram-Predicted Probability of 3-Year Survival",ylab="Actual 3-Year Survival",col="blue")
lines(cal,c('mean.predicted','KM'),type = 'a',lwd = 3,col ="black" ,pch = 16)
box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
dev.off()

#5-year
cox1 <- cph(Surv(futime,fustat) ~stage+Risk_levels,surv=T,x=T, y=T,time.inc = 5*365,data=rt) 
cal <- calibrate(cox1, cmethod="KM", method="boot", u=5*365, m= 10, B=1000)
#pdf("calibration5.pdf",12,8)
par(mar = c(10,5,3,2),cex = 1.0)
plot(cal,lwd=3,lty=2,errbar.col="black",xlim = c(0,1.0),ylim = c(0,1.0),xlab ="Nomogram-Predicted Probability of 5-Year Survival",ylab="Actual 5-Year Survival",col="blue")
lines(cal,c('mean.predicted','KM'),type = 'a',lwd = 3,col ="black" ,pch = 16)
box(lwd = 1)
abline(0,1,lty = 3,lwd = 3,col = "black")
dev.off()

library(timeROC)
library(ggplot2)
library(dplyr)


ROC.nom<-timeROC(T=rt$futime,delta=rt$fustat,
                 weighting="marginal", iid = TRUE,
                   marker=rt$risknom,cause=1,
                   times=quantile(rt$futime,probs=seq(0.2,0.8,0.1)),ROC=T)
ROC.nom

plotAUCcurve(ROC.nom, conf.int = F, col="darkcyan")

ROC.stage<-timeROC(T=rt$futime,delta=rt$fustat,
                 weighting="marginal", iid = TRUE,
                 marker=as.factor(rt$stage),cause=1,
                 times=quantile(rt$futime,probs=seq(0.2,0.8,0.1)),ROC=T)

ROC.risk<-timeROC(T=rt$futime,delta=rt$fustat,
                 weighting="marginal", iid = TRUE,
                 marker=as.factor(rt$Risk_levels),cause=1,
                 times=quantile(rt$futime,probs=seq(0.2,0.8,0.1)),ROC=T)

plotAUCcurve(ROC.nom, conf.int = F, col="green")

plotAUCcurve(ROC.stage, conf.int = F, col="red",add=TRUE)

plotAUCcurve(ROC.risk, conf.int = F, col="blue",add=TRUE)


legend("bottomright",c('stage',"risk_levels",'nomogram'),col=c('red','blue','green'),lty=1,lwd=2)
dev.off()
