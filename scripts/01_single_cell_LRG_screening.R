rm(list=ls());gc()
options(stringsAsFactors = F) 
setwd("F:\\liulu\\singlecell\\")
library(devtools)
library(Seurat)
library(ggplot2)
library(clustree)
library(cowplot)
library(dplyr)
library(hdf5r) 
library(patchwork)
library(Seurat)
assays = dir("GSE279998/")
assays
assay <- assays[2:17]
assay
dir <- paste0("GSE279998/",assay)

samples_name = c( "GSM858384", "GSM858385", "GSM858386" ,"GSM858387", "GSM858388" ,"GSM858389", "GSM858390",
                  "GSM858391" ,"GSM858392", "GSM858393" ,"GSM858394" ,"GSM858395", "GSM858396" ,"GSM858397",
                  "GSM858398", "GSM858399" )

scRNAlist <- list()
for (i in 1:length(dir)){
  counts =  Read10X(data.dir = dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts = counts,project = samples_name[i],min.cells = 3, min.features = 200)
  } 

dir.create("Integrate")
setwd("./Integrate")
names(scRNAlist)<-samples_name
scRNA=merge(x=scRNAlist[[1]],
            y=scRNAlist[ -1 ],add.cell.ids = samples_name)
scRNA<-JoinLayers(scRNA)   #5.0 need this
library(mgsub)
sampleType<-mgsub(scRNA@meta.data$orig.ident,c( "GSM858384", "GSM858385", "GSM858386" ,"GSM858387", "GSM858388" ,"GSM858389", "GSM858390",
                                                "GSM858391" ,"GSM858392", "GSM858393" ,"GSM858394" ,"GSM858395", "GSM858396" ,"GSM858397",
                                                "GSM858398", "GSM858399"), 
                        c( "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "SCC", "ADS","ADC","ADC","ADC"), recycle = T)
scRNA<-AddMetaData(scRNA,sampleType,col.name = "SampleType")
table(scRNA$orig.ident,scRNA$SampleType)

scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
HB.genes <- c("HBA1","HBA2","HBB","HBD","HBE1","HBG1","HBG2","HBM","HBQ1","HBZ")
HB_m <- match(HB.genes, rownames(scRNA@assays$RNA)) 
HB.genes <- rownames(scRNA@assays$RNA)[HB_m] 
HB.genes <- HB.genes[!is.na(HB.genes)] 
scRNA[["percent.HB"]]<-PercentageFeatureSet(scRNA, features=HB.genes) 
col.num <- length(levels(scRNA@active.ident))
violin <- VlnPlot(scRNA,
                  features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                  cols =rainbow(col.num), 
                  pt.size = 0, #不需要显示点，可以设置pt.size = 0
                  ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin
dev.off()
plot1 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(scRNA, feature1 = "nCount_RNA", feature2 = "percent.HB")
plot3
pearplot <- CombinePlots(plots = list(plot1, plot2, plot3), nrow=1, legend="none") 
pearplot
dev.off()
theme.set2 = theme(axis.title.x=element_blank())
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.HB")
group = "orig.ident"
# 质控前小提琴图
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA, group.by=group, pt.size = 0,
                       features = plot.featrures[i]) + theme.set2 + NoLegend()}
violin <- wrap_plots(plots = plots, nrow=2) 
#tiff("orig.ident_vlnplot.tiff",width=40,height = 38,units = "cm",pointsize = 12,res = 400)
violin
dev.off()

minGene=200
maxGene=5000
pctMT=10
minCounts=200
scRNA <- subset(scRNA, subset = nFeature_RNA > minGene &percent.mt < pctMT&nCount_RNA>minCounts)
col.num <- length(levels(scRNA@active.ident))
violin <-VlnPlot(scRNA,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"), 
                 cols =rainbow(col.num), 
                 pt.size = 0, 
                 ncol = 4) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 
violin
dev.off()
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)
saveRDS(scRNA1, file="scRNArsh1.rds")
#--------jiangweijulei------
library(Seurat)
library(tidyverse)
library(patchwork)
dir.create("cluster")
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000) 
top10 <- head(VariableFeatures(scRNA), 10) 
plot1 <- VariableFeaturePlot(scRNA) 
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5) 
plot <- CombinePlots(plots = list(plot1, plot2),legend="bottom") 
plot
dev.off()
##如果内存足够最好对所有基因进行中心化
scale.genes <-  rownames(scRNA)
scRNA <- ScaleData(scRNA, features = scale.genes)

library(harmony)
scRNA <- RunPCA(scRNA, features = VariableFeatures(scRNA),reduction.name = "pca") 
scRNA <- RunHarmony(scRNA, group.by.vars = c("orig.ident"),reduction.save = "harmony")
plot1 <- DimPlot(scRNA, reduction = "pca", group.by="orig.ident") 
plot2 <- ElbowPlot(scRNA, ndims=50, reduction="pca") 
plot3 <- DimPlot(scRNA, reduction = "harmony", group.by="orig.ident") 
plotc <- plot1+plot3
plotc
dev.off()


pc.num=1:20
scRNA1 <-scRNA %>% 
  RunUMAP(reduction = "pca", dims = pc.num) %>% 
  FindNeighbors(reduction = "pca", dims = pc.num) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

scRNA1 <- scRNA1 %>% 
  RunTSNE(reduction = "pca", dims = pc.num)

table(scRNA@meta.data$seurat_clusters)
DimPlot(scRNA1, reduction = "tsne", group.by = "orig.ident",   pt.size=0.5, label = F)

#----------细胞聚类--------
pc.num=1:20
scRNA <-scRNA %>% 
  RunUMAP(reduction = "harmony", dims = pc.num) %>% 
  FindNeighbors(reduction = "harmony", dims = pc.num) %>% 
  FindClusters(resolution = 0.8) %>% 
  identity()

scRNA <- scRNA %>% 
  RunTSNE(reduction = "harmony", dims = pc.num)

table(scRNA@meta.data$seurat_clusters)
# 创建图片对象
p3 <- DimPlot(scRNA, reduction = "tsne", group.by = "SampleType", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
p4 <- DimPlot(scRNA, reduction = "tsne", group.by = "orig.ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)

dev.off()
#-----------细胞类型鉴定-----------
dir.create("cell_identify")
diff.wilcox = FindAllMarkers(scRNA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "cell_identify/diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "cell_identify/top10_diff_genes_wilcox.csv", row.names = F)

# ##top10基因绘制热图
top10_genes <- read.csv("cell_identify/top10_diff_genes_wilcox.csv")
top10_genes = CaseMatch(search = as.vector(top10_genes$gene), match = rownames(scRNA))
plot1 = DoHeatmap(scRNA, features = top10_genes, group.by = "seurat_clusters", group.bar = T, size = 4)
ggsave("cell_identify/top10_markers.pdf", plot=plot1, width=20, height=30)
ggsave("cell_identify/top10_markers.png", plot=plot1, width=20, height=30)

###########SingleR鉴定细胞类型##########
library(SingleR)
library(celldex)
refdata <- get(load("E:/Desktop/singlecell/hyp/singleRref/HumanPrimaryCellAtlas_hpca.se_human.RData"))
refdata
testdata <- GetAssayData(scRNA, slot="counts")
clusters <- scRNA@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels, stringsAsFactors = F)

celltype$celltype[celltype$ClusterID==10]<-"Myeloid_cells"
celltype$celltype[celltype$ClusterID==15]<-"Myeloid_cells"
celltype$celltype[celltype$ClusterID==17]<-"Myeloid_cells"
celltype$celltype[celltype$ClusterID==21]<-"Myeloid_cells"
celltype$celltype[celltype$ClusterID==24]<-"Fibroblasts"
celltype$celltype[celltype$ClusterID==25]<-"Epithelial_cells"
celltype$celltype[celltype$ClusterID==27]<-"B_cell"
celltype$celltype[celltype$ClusterID==28]<-"NK_cell"
celltype$celltype[celltype$ClusterID==29]<-"Endothelial_cells"
celltype$celltype[celltype$ClusterID==32]<-"Fibroblasts"

celltype$celltype[celltype$celltype=="T_cell"]<-"T_cells"
celltype$celltype[celltype$celltype=="NK_cell"]<-"NK_cells"
celltype$celltype[celltype$celltype=="B_cell"]<-"B_cells"
write.csv(celltype,"cell_identify/celltype_singleR1.csv",row.names = F)
scRNA@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  scRNA@meta.data[which(scRNA@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}

table(scRNA$seurat_clusters,scRNA$celltype)
table(scRNA$celltype)
library(ggplot2) 

bcell.feat <- c("MS4A1", "CD79A")
nk.feat <- c("CCL5", "GNLY", "PRF1", "NKG7")
myeloid.feat <- c("FCER1G", "CSF1R", "CD68", "CD14")
DC.feat<-c("HLA-DQB1","HLA-DPB1","BIRC3")
neucell.feat <- c("FCGR3B")
t.feat <- c("CD3E", "CD3D")
fibro.feat <- c("THY1", "DCN", "LUM")
smc.feat <- c("COL6A2", "ACTA2", "CNN1", "MYH11")
endo.feat <- c( "EGFL7","ACKR1", "PECAM1", "VWF")
epithe.feature<-c("EPCAM","CDH1","KRT8")
feats <- c(endo.feat, neucell.feat, t.feat,DC.feat,nk.feat,fibro.feat, myeloid.feat, smc.feat, bcell.feat,epithe.feature)
seurat_object<-scRNA

seurat_object <- SetIdent(seurat_object, value="celltype")
DotPlot(seurat_object, assay = "RNA", features=feats, cols = "RdYlBu", dot.scale = 8) + #scale_colour_gradient(low = LIGHTGRAY, high = "#009933") + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1, size=12), axis.title = element_text(face="bold"))

DotPlot(scRNA, assay = "RNA", features=feats, cols = "RdYlBu", dot.scale = 8) + #scale_colour_gradient(low = LIGHTGRAY, high = "#009933") + 
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1, size=12), axis.title = element_text(face="bold"))
 

library(readxl)
lactate <- read_excel("lactate.xlsx")
gene<-lactate$Lactylate
markers <- list()
markers$lactate=lactate$Lactylate

library(UCell)
library(irGSEA)
scRNA1 <- irGSEA.score(object = scRNA, assay = "RNA", maxGSSize = 5000,
                       slot = "data", seeds = 123, ncores = 1, geneset = markers, 
                       custom = T, msigdb = F, 
                       method = c("AUCell", "UCell", "singscore"),  kcdf = 'Gaussian')


scRNA2 <- irGSEA.score(object = scRNA, assay = "RNA", 
                        slot = "data", seeds = 123, ncores = 1,
                        min.cells = 3, min.feature = 0,
                        custom = F, geneset = NULL, msigdb = T, 
                        species = "Homo sapiens", category = "C2",  
                        subcategory = NULL, geneid = "symbol",
                        method = c("AUCell", "UCell", "singscore", 
                                   "ssgsea"),
                        aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                        kcdf = 'Gaussian')


scRNAc3 <- irGSEA.score(object = scRNA, assay = "RNA", 
                        slot = "data", seeds = 123, ncores = 1,
                        min.cells = 3, min.feature = 0,
                        custom = F, geneset = NULL, msigdb = T, 
                        species = "Homo sapiens", category = "C5",  
                        subcategory = NULL, geneid = "symbol",
                        method = c("AUCell", "UCell", "singscore", 
                                   "ssgsea"),
                        aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                        kcdf = 'Gaussian')

scRNAH <- irGSEA.score(object = scRNA, assay = "RNA", 
                        slot = "data", seeds = 123, ncores = 1,
                        min.cells = 3, min.feature = 0,
                        custom = F, geneset = NULL, msigdb = T, 
                        species = "Homo sapiens", category = "H",  
                        subcategory = NULL, geneid = "symbol",
                        method = c("AUCell", "UCell", "singscore", 
                                   "ssgsea"),
                        aucell.MaxRank = NULL, ucell.MaxRank = NULL, 
                        kcdf = 'Gaussian')

result.degtypecell <- irGSEA.integrate(object = scRNA1, 
                                       group.by = "type", col.name = NULL,
                                       method = c("AUCell","UCell"))
library(Seurat)

Idents(scRNA) <- scRNA$celltype
ridgeplot <- irGSEA.ridgeplot(object = scRNA,
                              method = "UCell",
                              show.geneset = "lactate")
ridgeplot


Idents(scRNA) <- scRNA$SampleType

ridgeplot <- irGSEA.ridgeplot(object = scRNA,
                              method = "UCell",
                              show.geneset = "lactate")
ridgeplot


scatterplot <- irGSEA.density.scatterplot(object = scRNA1,
                                          method = "AUCell",
                                          show.geneset = "lactate",
                                          reduction = "tsne")
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)  # 用于统计检验和多重比较
library(rstatix) # 辅助统计检验

scRNA <- AddMetaData(
  object = scRNA,
  metadata = scRNA@assays$AUCell$scale.data,  # 你的 Aucell 打分矩阵（行=细胞ID，列=打分）
  col.name = "lactate"  # 自定义列名（后续代码用这个列名）
)

auc_data <- scRNA@meta.data %>%
  select(celltype, lactate) %>% 
  filter(!is.na(celltype) & !is.na(lactate))  

p1 <- ggplot(auc_data, aes(x = celltype, y = lactate, fill = celltype)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +  # 箱线图
  #geom_jitter(size = 0.5, alpha = 0.3, color = "black") +  # 散点（避免点重叠）
  labs(x = "celltype", y = "lactate", fill = "celltype") +  # 替换基因集名称
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # x轴标签旋转
    legend.position = "none",  # 移除重复图例
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

print(p1)
cell_stats <- auc_data %>%
  group_by(celltype) %>%  # 按样本类型（ADC/ADS/SCC）分组
  summarise(
    样本数 = n(),  # 每组的样本数量
    平均值 = mean(lactate, na.rm = TRUE),  # 平均值（lactate为Aucell打分列）
    标准差 = sd(lactate, na.rm = TRUE),    # 标准差（反映数据离散程度）
    中位数 = median(lactate, na.rm = TRUE),# 中位数（配合箱线图参考）
    最小值 = min(lactate, na.rm = TRUE),   # 最小值
    最大值 = max(lactate, na.rm = TRUE),   # 最大值
    .groups = "drop"  # 取消分组，返回普通数据框
  ) %>%
  # 可选：保留3位小数，让结果更简洁
  mutate(
    across(c(平均值, 标准差, 中位数, 最小值, 最大值), ~round(., 3))
  )

cell_stats




auc_data1 <- scRNA@meta.data %>%
  select(SampleType, lactate) %>% 
  filter(!is.na(SampleType) & !is.na(lactate))  

# 基础箱线图（显示中位数、四分位距、异常值）
p <- ggplot(auc_data1, aes(x = SampleType, y = lactate, fill = SampleType)) +
  geom_boxplot(alpha = 0.7, width = 0.6) +  # 箱线图
  #geom_jitter(size = 0.5, alpha = 0.3, color = "black") +  # 散点（避免点重叠）
  labs(x = "sample_type", y = "lactate", fill = "celltype") +  # 替换基因集名称
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  # x轴标签旋转
    legend.position = "none",  # 移除重复图例
    plot.title = element_text(hjust = 0.5, size = 12, face = "bold")
  )

print(p1)

kruskal_result <- auc_data1 %>%
  kruskal_test(lactate~ SampleType) 

kruskal_result <- auc_data1 %>% kruskal_test(lactate ~ SampleType)
cat("\nKruskal-Wallis 整体差异检验结果：\n")
print(kruskal_result)

wilcox_pairwise_all <- auc_data1 %>%
  pairwise_wilcox_test(
    lactate ~ SampleType,
    p.adjust.method = "fdr",
    paired = FALSE,
    alternative = "two.sided"
  ) 


my_comparisons <- list(c("ADC","ADS"), c("ADC", "SCC"),c("ADS", "SCC"))

p3 <- p + stat_compare_means(comparisons=my_comparisons,
                             label.y = c(0.4, 0.45, 0.5),
                             method="wilcox.test"
)
p3

p4 <- p + stat_compare_means(comparisons=my_comparisons,
                             label.y = c(0.4, 0.45, 0.5),
                             method="wilcox.test",
                             label="p.signif"
)
p3+p4
group_stats <- auc_data1 %>%
  group_by(SampleType) %>%  # 按样本类型（ADC/ADS/SCC）分组
  summarise(
    样本数 = n(),  # 每组的样本数量
    平均值 = mean(lactate, na.rm = TRUE),  # 平均值（lactate为Aucell打分列）
    标准差 = sd(lactate, na.rm = TRUE),    # 标准差（反映数据离散程度）
    中位数 = median(lactate, na.rm = TRUE),# 中位数（配合箱线图参考）
    最小值 = min(lactate, na.rm = TRUE),   # 最小值
    最大值 = max(lactate, na.rm = TRUE),   # 最大值
    .groups = "drop"  # 取消分组，返回普通数据框
  ) %>%
  # 可选：保留3位小数，让结果更简洁
  mutate(
    across(c(平均值, 标准差, 中位数, 最小值, 最大值), ~round(., 3))
  )
group_stats




a<-median(scRNA1@assays$AUCell$scale.data)
gene1_expr <-scRNA1@assays$AUCell$scale.data
# 创建共表达状态分类
expression_status <- ifelse(gene1_expr >=a, "high_level", "low_level")
# 添加到Seurat对象
status_colname <- paste0("lactate", "_status")
scRNA1[[status_colname]] <- expression_status

# 正确设置命名向量
status_names <- c("high_level", "low_level")
status_colors <- c( "#FF5733", "#338DFF")
# status_colors <- c("#F0F0F0", "#F0F0F0", "#F0F0F0", "#9933FF")
colors <- setNames(status_colors, status_names)
print(table(expression_status,scRNA1$SampleType))
# 绘制特征图
p <- DimPlot(scRNA1, reduction ="tsne", group.by = status_colname,
             cols = colors, pt.size = 1) 

# 显示图形
print(p)


scRNA <- readRDS("D:/liulu/singlecell/Integrate/scRNArsh1.rds")
table(scRNA$celltype,scRNA$lactate_status)
scRNA1$celltype.group <- paste(scRNA1$celltype, scRNA1$lactate_status, sep = "_")

Idents(scRNA1) <- "celltype.group"
table(scRNA1$celltype.group)

cellfordeg<-levels(scRNA1$celltype)
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(scRNA1, ident.1 = paste0(cellfordeg[i],"_high_level"), ident.2 = paste0(cellfordeg[i],"_low_level"), verbose = T, test.use = 'wilcox',min.pct = 0.1)
  write.csv(CELLDEG,paste0("cell_identify/",cellfordeg[i],".csv"))
}

list.files()
Idents(scRNA) <- "lactate_status"
mydeg <- FindMarkers(scRNA,ident.1 = 'high_level',ident.2 = 'low_level', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)
mydeg$symbol<-rownames(mydeg)
write.csv(mydeg,"cell_identify/alllactate.csv")

mydeg_pos <- mydeg %>%
  dplyr::filter(p_val_adj <0.01)%>%dplyr::filter(avg_log2FC>1)%>%
  dplyr::arrange(desc(avg_log2FC))

a<-median(scRNA@assays$AUCell$scale.data)
gene1_expr <-scRNA@assays$AUCell$scale.data

########相关性计算####
exprSet <- scRNA@assays$RNA$scale.data
exprSet<-as.data.frame(t(exprSet))#转置

library(dplyr)
# 2. 计算某个基因和其它基因的相关性（以S100A8为例）-----#####
#x <- as.numeric(exprSet[,"RPS7"])
length(x)

y<-gene1_expr
y
length(y)
colnames<-colnames(exprSet)
cor_data_df<- data.frame(colnames)
for(i in 1:length(colnames)){
  test<-cor.test(as.numeric(exprSet[,i]),y,type="spearman") #可更换pearson
  cor_data_df[i,2]<- test$estimate
  cor_data_df[i,3]<- test$p.value
  
}
names(cor_data_df)<-c("symbol","correlation","pvalue")
cor_data_df %>% head()
cor_data_df$cor<-abs(cor_data_df$correlation)
write.csv(cor_data_df,"cell_identify/allcorlactate.csv")
a<-cor.test(as.numeric(x),y,type="spearman")

# 3. 筛选有意义的正相关和负相关的基因-----####
library(dplyr)
library(tidyr)
cor_data_sig_pos <- cor_data_df %>%
  dplyr::filter(pvalue <0.01)%>%dplyr::filter(correlation >0.5)%>%
  dplyr::arrange(desc(correlation))

gene1<-cor_data_sig_pos$symbol
mydeg_pos$symbol<-rownames(mydeg_pos)
gene2<-mydeg_pos$symbol
a1<-intersect(mydeg_pos$symbol,cor_data_sig_pos$symbol)
gene<-append(gene1,gene2)
write.csv(gene,"cell_identify/lactate-related.csv")
write.table(gene1,"cell_identify/deg_lactate-related.tsv")
write.csv(gene2,"cell_identify/cor_lactate-related.csv")
intersect(gene1,gene2)

length(unique(gene1))
length(unique(gene2))
length(unique(gene))
table(scRNA$SampleType)


scRNA$celltype.group <- paste(scRNA$celltype, scRNA$lactate_status, sep = "_")
Idents(scRNA) <- "celltype.group"
table(scRNA$celltype.group)
expr <- AverageExpression(scRNA, assays = "RNA", slot = "data")[[1]]
expr <- expr[rowSums(expr)>0,]  #过滤细胞表达量全为零的基因
expr <- as.matrix(expr)
head(expr)

library(msigdbr)
msigdbr_species() #列出有的物种

H_df_all <- msigdbr(species = "Homo sapiens",
                    category = "H")  
H_df <- dplyr::select(H_df_all, gs_name, gene_symbol, gs_exact_source, gs_subcat)
H_list <- split(H_df$gene_symbol, H_df$gs_name) 


library(GSVA)
gsvaParh <- gsvaParam(as.matrix(expr), H_list,maxDiff = TRUE)
gsvaParh 
expr_geneseth <- gsva(gsvaParh)
colnames(expr_geneseth)


library(pheatmap)
pheatmap(expr_geneseth, show_colnames = T, 
         scale = "row",angle_col = "45",
         cluster_row = T,cluster_col = T,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50))



library(limma)

type = as.data.frame(colnames(expr_geneseth))
original_names <- type$`colnames(expr_geneseth)`
group_names <- sapply(strsplit(original_names, "-"), function(x) {
  if (length(x) >= 2) {
    trimws(x[length(x)-1])  # 去除可能的空格
  } else {
    warning(paste("列名", x, "拆分后长度不足2，无法提取组名"))
    NA
  }
})
type$group<-group_names
colnames(type)[1]<-"type"
group_list<-type
  
high = expr_geneseth[,grep('high',group_list$group)]
low  =  expr_geneseth[,grep('low',group_list$group)]
gsvahall<-cbind(high,low)

group_list = c(rep('high',ncol(high)),
               rep('low',ncol(low)))


design <- model.matrix(~ 0 + factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(gsvahall)

# 构建差异比较矩阵
contrast.matrix <- makeContrasts(high-low, levels = design)

# 差异分析，case vs. con
fit <- lmFit(gsvahall, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
allDiff=topTable(fit2,adjust='fdr',number=100000)
diff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
head(diff)

library(tidyverse)
diff$group <- ifelse( diff$logFC > 0 & diff$P.Value < 0.05 ,"up" ,
                      ifelse(diff$logFC < 0 & diff$P.Value < 0.05 ,"down","noSig")
)

diff2 <- diff %>% 
  mutate(hjust2 = ifelse(t>0,1,0)) %>% 
  mutate(nudge_y = ifelse(t>0,-0.1,0.1)) %>% 
  #filter(group != "noSig") %>% #注释掉该行即可
  arrange(t) %>% 
  rownames_to_column("ID")

diff2$ID <- factor(diff2$ID, levels = diff2$ID)
limt = max(abs(diff2$t))

ggplot(diff2, aes(ID, t,fill=group)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  scale_fill_manual(breaks=c("down","up"), #设置颜色
                    values = c("#008020","#08519C"))+
  geom_text(data = diff2, aes(label = diff2$ID, #添加通路标签
                              y = diff2$nudge_y),
            nudge_x =0,nudge_y =0,hjust =diff2$hjust,
            size = 3)+ #设置字体大小
  labs(x = "Pathways", #设置标题 和 坐标
       y=c("t value of GSVA score\n"),
       title = "GSVA")+
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
        legend.position = "none" #去掉legend
  )
write.table(diff2,"GSVA.tsv",sep = "\t")

save.image("liulu.Rdata")

library(AnnotationHub)	#library导入需要使用的数据包
library(org.Hs.eg.db)   #人类注释数据库
library(clusterProfiler)
library(dplyr)
library(ggplot2)
library(GOplot)
genes<-unique(gene)
eg <- bitr(genes, 
           fromType="SYMBOL", 
           toType=c("ENTREZID","ENSEMBL",'SYMBOL'),
           OrgDb="org.Hs.eg.db")


head(eg)
go <- enrichGO(eg$ENTREZID, 
               OrgDb = org.Hs.eg.db, 
               ont='ALL',
               pAdjustMethod = 'BH',
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,
               keyType = 'ENTREZID',
               readable = T)

dim(go)
ego2<-data.frame(go)
write.table(ego2,"go.tsv",sep = "\t")
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
  slice_min(p.adjust, n = 10) %>%  # 按p.adjust升序取前10（最显著）
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
write.table(go_plot_data,"goplot.tsv",sep="\t")



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
write.table(ekegg,"kegg.tsv",sep="\t")
#可视化，和上面的一样

#tiff(file="DOTKEGG.tiff",width = 15,height =10,units ="cm",compression="lzw",bg="white",res=400)
dotplot(kegg,font.size=8, showCategory=20)	# 画气泡图
dev.off()



