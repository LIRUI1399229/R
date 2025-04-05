
library(GEOmirror)
library(tidyverse)
library(stringr)
library(GEOquery)
library(dplyr)
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(WGCNA)
#参数设置
GSE="GSE128708"    
C="Control"          
P="COPD"              
Ccol="blue"        
Pcol="red"        
gset = getGEO('GSE128708', destdir=".", AnnotGPL = F, getGPL = F)
exprSet1 <- read.table(file = 'GSE128708_series_matrix.txt.gz',sep = '\t', header = T,quote = '""',fill = T,comment.char = "!") 
rownames(exprSet1) = exprSet1[,1] 
exprSet1 <- exprSet1[,-1]        
exp=exprSet1
class(gset)
###提取子集
gset[[1]]
#通过pData函数获取分组信息#
pdata <- pData(gset[[1]])
table(pdata$title)
write.table(pdata,file=paste0("pdata222.txt"),sep="\t",quote=F,col.names = T)
pdata2=read.table("sample.txt",header=F,row.names = 1)
pdata=pdata2
#设置参考水平,
group_list <- ifelse(str_detect(pdata$V2, "Control"), "Control",
                     "COPD")
#因子型
group_list = factor(group_list,
                    levels = c("Control","COPD"))
group_list
table(group_list)
##2.2 通过exprs函数获取表达矩阵并校正
exp <- exprs(gset[[1]])
sample=read.table("sample.txt",header=F,row.names = 1)
exp=exp[,rownames(sample)]


a=exp
library(hgu133plus2.db)
ids=toTable(hgu133plus2SYMBOL)
exprSet=a
table(rownames(exprSet) %in% ids$probe_id)
dim(exprSet)
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]
dim(exprSet)
ids=ids[match(rownames(exprSet),ids$probe_id),]
head(ids)
exprSet[1:5,1:5]

kelly <- function(exprSet,ids){
  tmp = by(exprSet,ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}
new_exprSet <- kelly(exprSet,ids)
exprSet=new_exprSet
exp=exprSet
rt=exp
######################################################################################################################1.数据准备
#分组
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
rt=rt[,rownames(sample)]
afcon=sum(sample[,1]==C)
#判断原始数据是否去了log
max(rt)
if(max(rt)>50) rt=log2(rt+1)     #rt最大值大于30则取log

#使用normalizeBetweenArrays进行矫正，矫正后赋值为rt1
rt1=normalizeBetweenArrays(as.matrix(rt))
range(rt1)
#未标准化,mar调整画布范围下左上右，自行调整哈
cols=rainbow(ncol(rt)) ###针对24个样本，设置颜色，整体呈现彩虹色
pdf(file = "1.raw.pdf",width=10,height = 8)
par(cex = 0.7,mar=c(8,8,8,8))
if(ncol(rt)>40) par(cex = 0.5,mar=c(8,8,8,8))   ###设置字体大小
boxplot(rt,las=2,col =cols ) ###绘图
dev.off()

#标准化
cols=rainbow(ncol(rt1)) ###针对24个样本，设置颜色，整体呈现彩虹色
pdf(file = "1.nor.pdf",width=10,height = 8)
par(cex = 0.5,mar=c(8,8,8,8))
if(ncol(rt1)>40) par(cex = 0.5,mar=c(8,8,8,8))   ###设置字体大小
boxplot(rt1,las=2,col =cols ) ###绘图
dev.off()

#保存标准化后结果
rt2=rbind(ID=colnames(rt1),rt1)
write.table(rt2,file=paste0("1.","norexp_",GSE,".txt"),sep="\t",quote=F,col.names = F)

#保留原始结果
rt3=rbind(ID=colnames(rt),rt)
write.table(rt3,file=paste0("1.","rawexp_",GSE,".txt"),sep="\t",quote=F,col.names = F)

#######################################################################################################################2.差异分析
exp=rt1
rt=exp
#####差异分析#####
library(tidyverse)
library(GEOquery)
#差异分析
library(limma)
design=model.matrix(~0+group_list)
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exp)
design
contrast.matrix<-makeContrasts(COPD-Control,levels = design)
contrast.matrix
#线性拟合模型构建#
fit <- lmFit(exp,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2) 
deg<-topTable(fit2, coef=1, n=Inf)
write.table(deg, file = "deg_all.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

##标记上下调基因
#logFC=1
#P.Value = 0.05
k1 = (deg$adj.P.Val < 0.05)&(deg$logFC < -0.585)
k2 = (deg$adj.P.Val < 0.05)&(deg$logFC > 0.585)
deg$change = ifelse(k1,"down",ifelse(k2,"up","not"))
table(deg$change)#计数
##热图##
cg = rownames(deg)[deg$change !="not"]
diff=exp[cg,]
library(pheatmap)
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=colnames(diff) 
pheatmap(diff,
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = F,
         show_colnames =T,
         cluster_cols = F,
         color = colorRampPalette(c("navy", "white", "red"))(50),
         fontsize = 10,
         fontsize_row=5,
         fontsize_col=5)
dev.off()
##火山图##
library(ggplot2)
head(deg)
ggplot(deg,aes(x=logFC,y=-log10(P.Value)))+
  geom_point(aes(color=change))+
  scale_colour_manual(values = c("dodgerblue","gray","firebrick"))+
  theme_bw()+
  labs(title="differentially expressed mRNAs",y="-log10(adj.P.Val)",x="log(Fold Change)")
dev.off()


#########################WGCNA#########################
###准备
afdir <- paste0(getwd(),"/5.WGCNA")           #路径必须中文
dir.create(afdir)

traitData=sample
traitData[,2]=traitData[,1]
traitData[,1]=ifelse(traitData[,1]==C,1,0)
traitData[,2]=ifelse(traitData[,2]==P,1,0)
#修改性状名称
colnames(traitData)=c(C,P)

###############正式分析
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
#Read in the female liver data set
datExpr0 = data.frame(t(rt))
colnames(datExpr0) <- rownames(rt)
rownames(datExpr0) <- colnames(rt)
datExpr1<-datExpr0
#筛选筛选方差前25%的基因
m.vars=apply(datExpr0,2,var)
expro.upper=datExpr0[,which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25))[4])]
datExpr1<-data.matrix(expro.upper)
datExpr0<-data.frame(datExpr1)
##check missing value
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

##filter过滤
meanFPKM=0.5  ####the threshold can be changed---过滤标准
n=nrow(datExpr0)
datExpr0[n+1,]=apply(datExpr0[c(1:nrow(datExpr0)),],2,mean)
datExpr0=datExpr0[1:n,datExpr0[n+1,] > meanFPKM]  # for meanFpkm in row n+1 and it must be above what you set--select meanFpkm>opt$meanFpkm(by rp)
filtered_fpkm=t(datExpr0)
filtered_fpkm=data.frame(rownames(filtered_fpkm),filtered_fpkm)
names(filtered_fpkm)[1]="sample"
head(filtered_fpkm)
write.table(filtered_fpkm, file=paste0(afdir,"/FPKM_filter.xls"),row.names=F, col.names=T,quote=FALSE,sep="\t")

sampleTree = hclust(dist(datExpr0), method = "average")
pdf(file =paste0(afdir,"/1_sampleClustering.pdf"), width = 12, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

### Plot a line to show the cut
abline(h = 400, col = "red")##是否选择剪切
dev.off()

### Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 400, minSize = 30)
table(clust)


### clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr0 = datExpr0[keepSamples, ]

#
nGenes = ncol(datExpr0)
nGenes
nSamples = nrow(datExpr0)
nSamples
dim(datExpr0)
datExpr0<-datExpr0 %>% filter(rownames(datExpr0)%in% rownames(sample))
nSamples = nrow(datExpr0)
nSamples
dim(datExpr0)
#
####载入性状数据
#Loading clinical trait data
for (df in colnames(traitData)) {
  traitData[,df]=traitData[,df]/max(traitData[,df])
  print(sd(traitData[,df]))
}
max(traitData)
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
rownames(datTraits) 
collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.

#sizeGrWindow(12,12)
pdf(file=paste0(afdir,"/2_Sample dendrogram and trait heatmap.pdf"),width=12,height=11)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()

#############################network constr########################################
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads()
# Choose a set of soft-thresholding powers
powers = c(c(1:10),seq(from=12,to=20,by=2))

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results:
#sizeGrWindow(9, 5)
pdf(file=paste0(afdir,"/3_Scale independence.pdf"),width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.9,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


######chose the softPower
softPower =sft$powerEstimate
#显示软阈值
print(softPower)

adjacency = adjacency(datExpr0, power = softPower)

##### Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)

#sizeGrWindow(12,9)
pdf(file=paste0(afdir,"/4_Gene clustering on TOM-based dissimilarity.pdf"),width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 100
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file=paste0(afdir,"/5_Dynamic Tree Cut.pdf"),width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file=paste0(afdir,"/6_Clustering of module eigengenes.pdf"),width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25######模块剪切高度
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()


# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#sizeGrWindow(12, 9)
pdf(file=paste0(afdir,"/7_merged dynamic.pdf"), width = 9, height = 6.5)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# Save module colors and labels for use in subsequent parts
#save(MEs, TOM, dissTOM,  moduleLabels, moduleColors, geneTree, sft, file = "networkConstruction-stepByStep.RData")




##############################relate modules to external clinical triats######################################
# Define numbers of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(10,6)
pdf(file=paste0(afdir,"/8_Module-trait relationships.pdf"),width=7,height=7.5)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


######## Define variable weight containing all column of datTraits

###MM and GS


# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#names of those trait
traitNames=names(datTraits)

geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")


####plot MM vs GS for each trait vs each module


##########example:royalblue and CK
#module="royalblue"
#column = match(module, modNames)
#moduleGenes = moduleColors==module

#trait="CK"
#traitColumn=match(trait,traitNames)

#sizeGrWindow(7, 7)

######

for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){
      
      #sizeGrWindow(7, 7)
      pdf(file=paste(afdir,"/9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}

#####
names(datExpr0)
probes = names(datExpr0)


#################export GS and MM############### 

geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.table(geneInfo, file = paste0(afdir,"/10_GS_and_MM.xls"),sep="\t",row.names=F)



####################################################Visualizing the gene network#######################################################


nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)


# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA



#随机选择基因数展示TOM热图
nSelect = 500
# For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select, select]
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select]

# Open a graphical window
#sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7
diag(plotDiss) = NA

pdf(file=paste0(afdir,"/13_Network heatmap plot_selected genes.pdf"),width=9, height=9)
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes", col=gplots::colorpanel(250,'red',"orange",'lemonchiffon'))
dev.off()



####################################################Visualizing the gene network of eigengenes####################################################


#sizeGrWindow(5,7.5)
pdf(file=paste0(afdir,"/14_Eigengene dendrogram and Eigengene adjacency heatmap.pdf"), width=5, height=7.5)
par(cex = 0.9)
plotEigengeneNetworks(MEs, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle= 90)
dev.off()

####venn图####
write.csv(deg,file = "diff_up2.csv",row.names = T)#下载两个基因集合
####venn####
library(randomcoloR)
library(venn) 
#文件名
A="DIFF"
B="WGCNA"
C="MitoCarta3.0"
#构建一个列表
geneList=list()
rt=read.table(paste0(A,".txt"),header=F,sep="\t",check.names=F)      
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)      
uniqGene=unique(geneNames)                 
geneList[[A]]=uniqGene                   
uniqLength=length(uniqGene)
print(paste("1",uniqLength,sep=" "))
rt=read.table(paste0(B,".txt"),header=F,sep="\t",check.names=F)    
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)       
uniqGene=unique(geneNames)                
geneList[[B]]=uniqGene
uniqLength=length(uniqGene)
print(paste("3",uniqLength,sep=" "))
rt=read.table(paste0(C,".txt"),header=F,sep="\t",check.names=F)    
geneNames=as.vector(rt[,1])                
geneNames=gsub("^ | $","",geneNames)       
uniqGene=unique(geneNames)                
geneList[[C]]=uniqGene
uniqLength=length(uniqGene)
print(paste("3",uniqLength,sep=" "))
mycol <- distinctColorPalette(3)
pdf(file="disease.pdf",width=5,height=5)                                                
venn(geneList,col=mycol[1:length(geneList)],zcolor=mycol[1:length(geneList)],box=F)
dev.off()
intersectGenes=Reduce(intersect,geneList)          
write.table(file="disease.txt",intersectGenes,sep="\t",quote=F,col.names=F,row.names=F)
###富集分析###
##GO##
library("org.Hs.eg.db")  
library("clusterProfiler")
library("enrichplot")
library("ggplot2")
library("ggnewscale")
library("enrichplot")
library("DOSE")
library(stringr)

pvalueFilter=0.05         
qvalueFilter=1  
showNum=10

rt=read.table("disease.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)

write.csv(rt,file = "rt.csv",row.names = T)
rt<-read.csv("rt.csv",header=T,row.names = 1)

colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}


kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO,file="GO.xls",sep="\t",quote=F,row.names = F)


if(nrow(GO)<30){
  showNum=nrow(GO)
}

pdf(file="GO_barplot.pdf",width = 9,height =9)
bar=barplot(kk, drop = TRUE, showCategory =showNum,split="ONTOLOGY",color = colorSel) + facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
print(bar)
dev.off()


pdf(file="GO_bubble.pdf",width = 12,height =9)
bub=dotplot(kk,showCategory = showNum, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
print(bub)
dev.off()
####KEGG####
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("pathview")
library("ggnewscale")
library("DOSE")
library(stringr)

pvalueFilter=0.05        
qvalueFilter=1        
showNum=20

rt=read.table("deg_all.txt",sep="\t",check.names=F,header=T)  
rownames(rt)=rt[,1]
rt=rt[read.table("DIFF.txt",sep="\t",check.names=F,header=F)[,1] ,c(1,2)]

genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","logFC","entrezID") 
#查看转换失败的ID
rt[rt[,"entrezID"]=="NA",]
rt=rt[rt[,"entrezID"]!="NA",]                        
gene=rt$entrezID
gene=unique(gene)
aflogfc=rt$logFC
names(aflogfc)=rt$symbol
colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG.xls",sep="\t",quote=F,row.names = F)

if(nrow(KEGG)<showNum){
  showNum=nrow(KEGG)
}

pdf(file="KEGG_barplot.pdf",width = 9,height = 7)
barplot(kk, drop = TRUE, showCategory = showNum, color = colorSel) +scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()

pdf(file="KEGG_bubble.pdf",width = 9,height = 7)
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))
dev.off()


Enrichment_KEGG <- enrichKEGG(gene = rt$entrezID,
                              organism = "hsa", 
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 0.05,
                              pAdjustMethod = "BH",
                              minGSSize = 10,
                              maxGSSize = 5000)
#将ENTREZ重转为symbol：
Enrichment_KEGG <- setReadable(Enrichment_KEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
#计算Rich Factor（富集因子）：
Enrichment_KEGG2 <- mutate(Enrichment_KEGG,
                           RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
#计算Fold Enrichment（富集倍数）：
Enrichment_KEGG2 <- mutate(Enrichment_KEGG2, 
                           FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
head(Enrichment_KEGG2@result)
str(Enrichment_KEGG2@result)
devtools::install_github("YuLab-SMU/GOSemSim")
devtools::install_github('GuangchuangYu/clusterProfiler')
packageVersion("clusterProfiler")

Enrichment_KEGG <- enrichKEGG(gene = rt$entrezID,
                              organism = "hsa", 
                              pvalueCutoff = 0.05,
                              qvalueCutoff = 1,
                              pAdjustMethod = "BH",
                              minGSSize = 10,
                              maxGSSize = 500)
#将ENTREZ重转为symbol：
Enrichment_KEGG <- setReadable(Enrichment_KEGG,OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
#计算Rich Factor（富集因子）：
Enrichment_KEGG2 <- mutate(Enrichment_KEGG,
                           RichFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
#计算Fold Enrichment（富集倍数）：
Enrichment_KEGG2 <- mutate(Enrichment_KEGG2, 
                           FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
head(Enrichment_KEGG2@result)
str(Enrichment_KEGG2@result)
as.data.frame(sort(table(Enrichment_KEGG2@result$category)))
as.data.frame(sort(table(Enrichment_KEGG2@result$subcategory)))
head(Enrichment_KEGG2@result[,1:3])
library(ggplot2)
p1 <- dotplot(Enrichment_KEGG2,x = "GeneRatio",color = "p.adjust",showCategory = 25) +
  facet_grid(rows = vars(category),scales = 'free_y',space = 'free_y') +
  scale_color_gradientn(colors = c('#BF1E27','#FEB466','#F9FCCB','#6296C5','#38489D'))

p1
p2 <- barplot(Enrichment_KEGG2,x = "GeneRatio",color = "p.adjust",showCategory = 25) +
  facet_grid(rows = vars(category),scales = 'free_y',space = 'free_y') +
  scale_fill_gradientn(colors = c('#BF1E27','#FEB466','#F9FCCB','#6296C5','#38489D'))
p2
library(ggplot2)
mytheme <- theme(axis.text.x = element_text(hjust = 0.5,size = 20), 
                 axis.ticks.y = element_blank(), ## 删去y轴刻度线
                 axis.text.y = element_text(size = 20), 
                 axis.title.x = element_text(size = 20), 
                 axis.title.y = element_text(size = 20), 
                 axis.line = element_line(size = 1),
                 plot.margin = unit(c(1,1,1,1), "cm"),#画布边缘距离上(top)、右(right)、下(bottom)、左(left) 
                 plot.title = element_text(hjust = 0.5,size =  22),
                 legend.title = element_blank(), 
                 legend.text = element_text(size = 22), 
                 legend.position = c(0.6,0.2),
                 legend.background = element_rect(fill = 'transparent'))
Enrichment_KEGG2 <- Enrichment_KEGG2@result
Enrichment_KEGG2 <- Enrichment_KEGG2[order(Enrichment_KEGG2$pvalue),]
dt <- Enrichment_KEGG2[order(Enrichment_KEGG2$Count, decreasing = T),][1:20,]
dt <- dt[order(dt$pvalue),]
dt$Description <- factor(dt$Description,levels = rev(dt$Description)) 
dt
range(-log10(as.numeric(dt$pvalue)))#[1] 8.027827 39.300064
range(as.numeric(dt$Count))#[1]18 46

library(stringr)
p  <- ggplot(dt, aes(x = Description, y = -log10(as.numeric(pvalue)),fill = factor(category))) +
  geom_bar(stat = "identity", width = 0.8) + #绘制条形图
  scale_fill_manual(values = c("#BEBADA","#B3DE69","#E6AB02","#80B1D3","red")) +
  coord_flip() +
  labs(title = NULL,x = NULL,y = bquote(~-Log[10]~italic("P-value"))) +
  geom_text(aes(label = GeneRatio),size = 2) +
  geom_hline(yintercept = -log10(as.numeric(0.05)),color = 'grey',size = 1,lty = 'dashed') + #添加虚线
  scale_y_continuous(expand = c(0,0),breaks = seq(0,45,5), limits = c(0,45,5)) + 
  scale_x_discrete(labels = function(dat) str_wrap(dat,width = 10)) +
  theme_bw() + mytheme
p
dev.off()

#####机器学习#####

write.table(file="GSE128708.txt",exp,sep="\t",quote=F,col.names=T,row.names=T)
####LASSO####
library(survival)
library(glmnet)
library(ggplot2)
library(ggsci)
library(patchwork)
library(limma)

inputFile="GSE128708.txt"       #输入文件
C="Control"                        #正常组样本名称

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
data=avereps(rt)
data=t(data)
data=data[,read.table("disease.txt", header=F, sep="\t", check.names=F)[,1]]
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
x=as.matrix(data)

#控制组放置最前
afcon=sum(sample[,1]==C)
group=c(rep("0",afcon),rep("1",nrow(data)-afcon))
group=as.matrix(group)
rownames(group)=rownames(data)
y=as.matrix(group[,1])

set.seed(12345)
cvfit = cv.glmnet(x, y,family = "binomial", nlambda=100, alpha=1,nfolds = 10) #这里alpha=1为LASSO回归，如果等于0就是岭回归，10乘交叉验证
#参数 family 规定了回归模型的类型：
#family="gaussian" 适用于一维连续因变量（univariate）
#family="mgaussian" 适用于多维连续因变量（multivariate）
#family="poisson" 适用于非负次数因变量（count）
#family="binomial" 适用于二元离散因变量（binary）
#family="multinomial" 适用于多元离散因变量（category）
#我们这里结局指标是2分类变量，所以使用binomial

fit <- glmnet(x,y,family = "binomial")
cvfit$lambda.min

#提取信息及预测风险
coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene, Coef=actCoef)
write.table(geneCoef, file="geneCoef.xls", sep="\t", quote=F, row.names=F)
write.table(file="lassoset.txt",lassoGene,sep="\t",quote=F,col.names=F,row.names=F) #文件名
#######################简单作图########################################
pdf("lasso.pdf",height = 6,width = 12)
layout(matrix(c(1,1,2,2), 2, 2, byrow = F))   #两行两列，图一占前俩格，图二占后两格，按列排
#pdf("lambda.pdf")
plot(fit,xvar = 'lambda')
#dev.off()
#pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
#dev.off()
dev.off()

####RF####
#引用包
library(randomForest)
library(limma)
library(ggpubr)
set.seed(2267299)
inputFile="GSE128708.txt"       #输入文件
C="Control"                        #正常组样本名称

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)

data=avereps(rt)
data=t(data)
data=data[,read.table("disease.txt", header=F, sep="\t", check.names=F)[,1]]
sample=read.table("sample.txt",sep="\t",header=F,check.names=F,row.names = 1)
data=data[rownames(sample),]
colnames(data)=gsub("-", "afaf", colnames(data))
#控制组放置最前
afcon=sum(sample[,1]==C)
group=c(rep("Control",afcon),rep("COPD",nrow(data)-afcon))

#随机森林树
rf=randomForest(as.factor(group)~., data=data, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()

#找出误差最小的点
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data, ntree=optionTrees)

#查看基因的重要性
#绘制基因的重要性图
importance=importance(x=rf2)
importance=as.data.frame(importance)
#importance$size=gsub("-", "afaf", importance$size)
importance$size=rownames(importance)
importance=importance[,c(2,1)]
names(importance)=c("Gene","importance")
#展示前20个基因的重要性
af=importance[order(importance$importance,decreasing = T),]
af=af[1:14,]
p=ggdotchart(af, x = "Gene", y = "importance",
             color = "importance", # Custom color palette
             sorting = "descending",                       # Sort value in descending order
             add = "segments",                             # Add segments from y = 0 to dots
             add.params = list(color = "lightgray", size = 2), # Change segment color and size
             dot.size = 6,                        # Add mpg values as dot labels
             font.label = list(color = "white", size = 9,
                               vjust = 0.5),               # Adjust label parameters
             ggtheme = theme_bw()         ,               # ggplot2 theme
             rotate=TRUE                                       )#翻转坐标轴 
p1=p+ geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  gradient_color(palette =c(ggsci::pal_npg()(2)[2],ggsci::pal_npg()(2)[1])      ) +#颜色
  grids()   
#保存图片
pdf(file="importance.pdf", width=6, height=6)
print(p1)
dev.off()
#挑选疾病特征基因
rfGenes=importance[order(importance[,"importance"], decreasing = TRUE),]
write.table(rfGenes, file="rfGenes.txt", sep="\t", quote=F, col.names=T, row.names=F)


####SVM####

write.csv(x,file = "deg.csv",row.names=T)
a=read.table("deg.txt", header=T, sep="\t", check.names=F,row.names = 1)

library(tidyverse)
library(glmnet)
source('msvmRFE.R')   #文件夹内自带
library(VennDiagram)
library(sigFeature)
library(e1071)
library(caret)
library(randomForest)
#library(e1071)
#source(msvmRFE.R)
train<-a
input <- train
set.seed(1)
#采用五折交叉验证 (k-fold crossValidation）
svmRFE(input, k = 5, halve.above = 100) #分割数据，分配随机数
nfold = 5
nrows = nrow(input)
folds = rep(1:nfold, len=nrows)[sample(nrows)]
folds = lapply(1:nfold, function(x) which(folds == x))
results = lapply(folds, svmRFE.wrap, input, k=5, halve.above=100) #特征选择
top.features = WriteFeatures(results, input, save=F) #查看主要变量
head(top.features)

#把SVM-REF找到的特征保存到文件
write.csv(top.features,"feature_svm.csv")

# 运行时间主要取决于选择变量的个数，一般的电脑还是不要选择太多变量
# 选前48个变量进行SVM模型构建，体验一下

featsweep = lapply(1:14, FeatSweep.wrap, results, input) #48个变量
featsweep

#load("featsweep.RData")
# 选前300个变量进行SVM模型构建，然后导入已经运行好的结果
#featsweep = lapply(1:300, FeatSweep.wrap, results, input) #300个变量
save(featsweep,file = "featsweep.RData")
#画图
no.info = min(prop.table(table(input[,1])))
errors = sapply(featsweep, function(x) ifelse(is.null(x), NA, x$error))

dev.new(width=2, height=4, bg='white')
pdf("B_svm-error.pdf",width = 10,height = 6)
PlotErrors(errors, no.info=no.info) #查看错误率
dev.off()

dev.new(width=2, height=4, bg='white')
pdf("B_svm-accuracy.pdf",width = 10,height = 6)
Plotaccuracy(1-errors,no.info=no.info) #查看准确率
dev.off()
# 图中红色圆圈所在的位置，即错误率最低点
which.min(errors) 
top<-top.features[1:which.min(errors), "FeatureName"]
write.csv(top,"top.csv")

###XGboost###
library(xgboost)
library(caret)
library(tidyverse)
library(readr)
library(VIM)
library(caret)
library(rpart)
library(rpart.plot)
library(Metrics)
library(stringr)
library(rpart)
library(tibble)
library(bitops)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(tidyverse)
library(limma)
library(pheatmap)
library(visNetwork)
library(ggpol)
library(ggplot2)
library(sparkline)
library(dplyr)
library(tidyverse)
library(caret)
library(DALEX)
library(gbm)

set.seed(99919)
data<-a
group=data$group
data=data[,-(1)]
# Fitting model(用caret实现)
TrainControl <- trainControl( method = "repeatedcv", number = 10, repeats = 4)
model<- train(x=data,y=as.factor(group),  method = "xgbTree", trControl = TrainControl,verbose = FALSE)
plot(varImp(model))
importance <- varImp(model)
head(importance)
important <- as.data.frame(importance$importance) 
xgb<-important
varimpdf <- data.frame(var = row.names(xgb),
                       impor = xgb[,1])
ggplot(varimpdf,aes(x = reorder(var,-impor), y = impor))+
  geom_col(colour = "lightblue",fill = "lightblue")+
  labs(title="Feature gene importance (XGBoost)", x="",y = "importance")+
  theme(plot.title = element_text(size=12,hjust=0.5))+
  theme(axis.text.x = element_text(size = 3))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 50,vjust = 0.85,hjust = 0.75))
write.table(file="xgb.txt",xgb,sep="\t",quote=F,col.names=T,row.names=T)
###GBM###
dev.off()
library(readr)
library(VIM)
library(caret)
library(rpart)
library(rpart.plot)
library(Metrics)
library(stringr)
library(rpart)
library(tibble)
library(bitops)
library(rattle)
library(rpart.plot)
library(RColorBrewer)
library(tidyverse)
library(limma)
library(pheatmap)
library(visNetwork)
library(ggpol)
library(ggplot2)
library(sparkline)
library(dplyr)
library(tidyverse)
library(caret)
library(DALEX)
library(gbm)
library(caret)
data<-a
group=data$group
data=data[,-(1)]
set.seed(1)
metric <- "RMSE"
myControl <- trainControl(method="cv", number=5)
# Fitting model
fitControl <- trainControl( method = "repeatedcv", number = 5, repeats = 5)
fit <- train(x=data,y=as.factor(group),  method = "gbm", trControl = fitControl,verbose = FALSE)
#绘制基因重要性梯度图
importances <- varImp(fit)
importances
importance <- as.data.frame(importances$importance)
#输入完上面那串代码后，显示的结果就是GBM结果，将他们复制出来。
#重要性筛选区域设置多少可以自己定，我这里是只要重要性不为0都可以
#删除为重要性为0的gene后重新导入
GBM<-importance
varimpdf <- data.frame(var = row.names(GBM),
                       impor = GBM[,1])

ggplot(varimpdf,aes(x = reorder(var,-impor), y = impor))+
  geom_col(colour = "lightblue",fill = "lightblue")+
  labs(title="Feature gene importance (Gradient Boosting Machine)", x="",y = "importance")+
  theme(plot.title = element_text(size=12,hjust=0.5))+
  theme(axis.text.x = element_text(size = 5))+
  theme(axis.text.y = element_text(size = 12))+
  theme(axis.text.x = element_text(angle = 50,vjust = 0.85,hjust = 0.75))

write.csv(file="GBM.csv",GBM,quote=F,row.names=T)


install.packages(
  "lightgbm"
  , type = "both"
  , repos = "https://cran.r-project.org"
)
#####LIGHTGBM#####
library(lightgbm)
library(rsample)
a=read.table("deg.txt", header=T, sep="\t", check.names=F,row.names = 1)
train<-a
input <- train
set.seed(12345)
# 数据加载
df <- a
data_train <- df
y <- data_train$group
x <- data.matrix(data_train[,-1])
dtrain <- lgb.Dataset(data = x,label = y)
params <-  list(    num_leaves = 4L,  
                    learning_rate = 0.05,   
                    max_depth=6,   
                    feature_fraction = 0.9,     
                    bagging_fraction = 0.8,   
                    bagging_freq = 5,  
                    objective = "binary")
fit <- lgb.train(  data = dtrain,  
                   params = params,  
                   nrounds = 100L,  
                   verbose = 0)
lgb.importance(fit) 
lgb.importance(fit) %>% lgb.plot.importance()

####ROC####
#因子型
##ROC
inputFile="GSE128708.txt"       #表达矩阵
hub="hub.txt"        #核心基因
#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
exp=rt
hubgenes=c("AKR1B10","VPS13D")
hubgenes_expression<-exp[match(hubgenes,rownames (exp)),]

library(pROC)

hubgenes_expression <- data.matrix(hubgenes_expression)
roc1<- roc(group_list, hubgenes_expression[1,])
roc1<- roc(controls=hubgenes_expression[1,][group_list=="Control"],
           cases=hubgenes_expression[1,][group_list=="COPD"])
roc2<- roc(group_list, hubgenes_expression[2,])
roc2<- roc(controls=hubgenes_expression[2,][group_list=="Control"],
           cases=hubgenes_expression[2,][group_list=="COPD"])
roc3<- roc(group_list, hubgenes_expression[3,])
roc3<- roc(controls=hubgenes_expression[3,][group_list=="Control"],
           cases=hubgenes_expression[3,][group_list=="COPD"])
round(auc(roc1),3)##AUC
round(ci(roc1),3)##95%CI
round(auc(roc2),3)##AUC
round(ci(roc2),3)##95%CI
round(auc(roc3),3)##AUC
round(ci(roc3),3)##95%CI

plot(roc1,
     add = FALSE,
     col = 'red', 
     legacy.axes = TRUE,
     main="AKR1B10",
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.x = 0.5,
     print.auc.y = 0.5)


plot(roc3,
     add = FALSE,
     col = 'red', 
     legacy.axes = TRUE,
     main="VPS13D",
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.x = 0.5,
     print.auc.y = 0.5)

inputFile="10006.txt"       #表达矩阵
hub="hub.txt"        #核心基因
#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F,row.names = 1)
rt=as.matrix(rt)
exp=rt
hubgenes=c("AKR1B10","HTATIP2","VPS13D")
hubgenes_expression<-exp[match(hubgenes,rownames (exp)),]

pdata2=read.table("10006group.txt",header=T,row.names = 1)
pdata=pdata2
#设置参考水平,
group_list <- ifelse(str_detect(pdata$group, "Control"), "Control",
                     "COPD")
#因子型
group_list = factor(group_list,
                    levels = c("Control","COPD"))
group_list
table(group_list)
library(pROC)
hubgenes_expression <- data.matrix(hubgenes_expression)
roc1<- roc(group_list, hubgenes_expression[1,])
roc1<- roc(controls=hubgenes_expression[1,][group_list=="Control"],
           cases=hubgenes_expression[1,][group_list=="COPD"])
roc2<- roc(group_list, hubgenes_expression[2,])
roc2<- roc(controls=hubgenes_expression[2,][group_list=="Control"],
           cases=hubgenes_expression[2,][group_list=="COPD"])
roc3<- roc(group_list, hubgenes_expression[3,])
roc3<- roc(controls=hubgenes_expression[3,][group_list=="Control"],
           cases=hubgenes_expression[3,][group_list=="COPD"])
round(auc(roc1),3)##AUC
round(ci(roc1),3)##95%CI
round(auc(roc2),3)##AUC
round(ci(roc2),3)##95%CI
round(auc(roc3),3)##AUC
round(ci(roc3),3)##95%CI

plot(roc1,
     add = FALSE,
     col = 'red', 
     legacy.axes = TRUE,
     main="AKR1B10",
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.x = 0.5,
     print.auc.y = 0.5)


plot(roc3,
     add = FALSE,
     col = 'red', 
     legacy.axes = TRUE,
     main="VPS13D",
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.x = 0.5,
     print.auc.y = 0.5)


#ROC2
##ROC2
inputFile="20257.txt"       #表达矩阵
hub="hub.txt" 
rt=read.table(inputFile, header=T, sep="\t", check.names=F,row.names = 1)
rt=as.matrix(rt)
rt=t(rt)
exp=rt
hubgenes=c("AKR1B10","HTATIP2","VPS13D")
hubgenes_expression<-exp[match(hubgenes,rownames (exp)),]
pdata2=read.table("20257group.txt",header=T,row.names = 1)
pdata=pdata2
#设置参考水平,
group_list <- ifelse(str_detect(pdata$group, "Control"), "Control",
                     "COPD")
#因子型
group_list = factor(group_list,
                    levels = c("Control","COPD"))
group_list
table(group_list)


library(pROC)

hubgenes_expression <- data.matrix(hubgenes_expression)
roc1<- roc(group_list, hubgenes_expression[1,])
roc1<- roc(controls=hubgenes_expression[1,][group_list=="Control"],
           cases=hubgenes_expression[1,][group_list=="COPD"])
roc2<- roc(group_list, hubgenes_expression[2,])
roc2<- roc(controls=hubgenes_expression[2,][group_list=="Control"],
           cases=hubgenes_expression[2,][group_list=="COPD"])
roc3<- roc(group_list, hubgenes_expression[3,])
roc3<- roc(controls=hubgenes_expression[3,][group_list=="Control"],
           cases=hubgenes_expression[3,][group_list=="COPD"])
round(auc(roc1),3)##AUC
round(ci(roc1),3)##95%CI
round(auc(roc2),3)##AUC
round(ci(roc2),3)##95%CI
round(auc(roc3),3)##AUC
round(ci(roc3),3)##95%CI

plot(roc1,
     add = FALSE,
     col = 'red', 
     legacy.axes = TRUE,
     main="AKR1B10",
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.x = 0.5,
     print.auc.y = 0.5)

plot(roc3,
     add = FALSE,
     col = 'red', 
     legacy.axes = TRUE,
     main="VPS13D",
     xlab = "1-Specificity",
     print.auc =TRUE,
     print.auc.x = 0.5,
     print.auc.y = 0.5)

write.table(t(hubgenes_expression), file="hubgenes_expression.txt",sep="\t",quote=F)
write.table(t(hubgenes_expression), file="hubgenes_expression2.txt",sep="\t",quote=F)

#引用包
library(dplyr)
library(pROC)
library(ggplot2)
library(survival)
library(regplot)
library(rms)
library(ggsci)
library(survminer)
library(timeROC)
library(ggDCA)
library(limma)

inputFile="GSE128708.txt"       #表达矩阵
hub="hub.txt"        #核心基因

#读取输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)
sample=read.table("sample.txt",sep="\t",header=F,check.names=F)
colnames(sample)=c("ID","Type")
data=data[sample$ID,]
aSAH1=data[,read.table(hub, header=F, sep="\t", check.names=F)[,1]]
aSAH=cbind(sample,aSAH1)
#简单看一下ROC曲线AUC的情况
aflist=roc(Type~AKR1B10+VPS13D, data = aSAH)
g3 <- ggroc(aflist, size = 1.2,alpha=.6,)
g5=g3+ggsci::scale_color_lancet()
print(g5)
#诺曼图，高度比：8比10
############逻辑回归模型
dd <- datadist(aSAH)
options(datadist="dd")
fit <- lrm(formula = Type ~ AKR1B10+VPS13D, data =aSAH)
print(fit)
coef=as.data.frame(fit$coefficients)[-1,,drop=F]
coefout=cbind(ID=rownames(coef),coef)
write.table(coefout,file="coefficients.txt",sep="\t",quote=F,row.names = F)
#绘图
pdf(file="nomogram.pdf", width=9, height=7.5)
plot(nomogram(fit,fun.at = seq(0.05,0.95,0.05)),funlabel = "nomogram model")
dev.off()

plot(regplot(fit,plots=c("density","boxes"), observation=T, title="Prediction Nomogram", clickable=T, points=TRUE,droplines=TRUE))

nomoscore=predict(fit, data=t(aSAH))
aSAH$nomoscore=nomoscore
write.table(aSAH,file="nomoscore.txt",sep="\t",quote=F,row.names = F)
write.table(aSAH1,file="HUBEXP.txt",sep="\t",quote=F,row.names = T)
#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
inputFile="HUBEXP2.txt"       #输入文件
rt=read.table(inputFile, header=T, sep="\t", check.names=F,row.names = 1)
library(tidyverse)
library(rms)
colnames(rt)
ddist <- datadist(rt)
options(datadist = "ddist")   #使用函数datadist()将数据打包

###4. Logstic回归
fit <- lrm(status~AKR1B10+VPS13D, data=rt, x=T, y=T)

#利用lrm(）函数对模型进行拟合
fit   #查看模型拟合结果

###5. 构建Nomogram
nom<- nomogram(fit, fun=plogis,
               fun.at=c(0.001,0.03,0.4,0.9,0.99),
               lp=F, funlabel="Risk of COPD")
plot(nom)
Cindex <- rcorrcens(status~predict(fit), data = rt)
Cindex

library(pROC)   #绘制ROC曲线
##绘制ROC曲线并给出阈值和ROC曲线下面积
gfit <- roc(status~predict(fit), data = rt)
plot(gfit,
     print.auc=TRUE,   #输出AUC
     main = "ROC CURVE",   #设置图形的标题
     col= "blue",   #曲线颜色
     print.thres.col="black",   #cut-off值字体的颜色
     identity.col="blue",   #对角线颜色
     identity.lty=1,identity.lwd=1,legacy.axes=TRUE)
###5. 绘制Calibration plot
cal <- calibrate(fit, method="boot", B=1000)
plot(cal,
     xlab="Nomogram-predicted probability of nonadherence",
     ylab="Actual diagnosed nonadherence (proportion)",
     sub=F)# 加载必要的库

library(rms)
formula <- as.formula("status~AKR1B10+VPS13D")
logistic_model <- lrm(formula, data = rt, x = TRUE, y = TRUE)
summary(logistic_model)
nom <- nomogram(logistic_model, fun = function(x) 1/(1 + exp(-x)), 
                funlabel = "Predicted Probability",
                fun.at=c(0.001,0.03,0.4,0.9,0.99),lp=F)
plot(nom)

inputFile="test10006.txt" 
#读取输入文件
test1=read.table(inputFile, header=T, sep="\t", check.names=F,row.names = 1)
external_predictions <- predict(fit, test1, type = "fitted")  
# 假设 external_data 是外部验证集
roc_external <- roc(test1$status, external_predictions)
auc_external <- auc(roc_external)
plot(roc_external,
     print.auc=TRUE,   #输出AUC
     main = "ROC CURVE",   #设置图形的标题
     col= "red",   #曲线颜色
     print.thres.col="black",   #cut-off值字体的颜色
     identity.col="red",   #对角线颜色
     identity.lty=1,identity.lwd=1,legacy.axes=TRUE)

inputFile="test11784.txt" 
#读取输入文件
test1=read.table(inputFile, header=T, sep="\t", check.names=F,row.names = 1)
external_predictions <- predict(fit, test1, type = "fitted")  
# 假设 external_data 是外部验证集
roc_external <- roc(test1$status, external_predictions)
auc_external <- auc(roc_external)
plot(roc_external,
     print.auc=TRUE,   #输出AUC
     main = "ROC CURVE",   #设置图形的标题
     col= "red",   #曲线颜色
     print.thres.col="black",   #cut-off值字体的颜色
     identity.col="red",   #对角线颜色
     identity.lty=1,identity.lwd=1,legacy.axes=TRUE)

inputFile="test20257.txt" 
#读取输入文件
test1=read.table(inputFile, header=T, sep="\t", check.names=F,row.names = 1)
external_predictions <- predict(fit, test1, type = "fitted")  
# 假设 external_data 是外部验证集
roc_external <- roc(test1$status, external_predictions)
auc_external <- auc(roc_external)
plot(roc_external,
     print.auc=TRUE,   #输出AUC
     main = "ROC CURVE",   #设置图形的标题
     col= "red",   #曲线颜色
     print.thres.col="black",   #cut-off值字体的颜色
     identity.col="red",   #对角线颜色
     identity.lty=1,identity.lwd=1,legacy.axes=TRUE)

setwd("C:/Users/Desktop/COPD/ssGSEA")
load("GSE128708_group3.rda")
exp3=read.table("1.norexp_GSE128708.txt", header=T, sep="\t", check.names=F,row.names = 1)
exprSet=exp3
library(GSEABase)
library(GSVA)
library(dplyr)
library(tibble)
library(tidyverse)
library(ggpubr)
geneSet <- read.csv("geneSet.csv",header = T) 
geneSet <- geneSet %>%
  column_to_rownames("X1")%>%t()
a <- geneSet
set <- colnames(a)
l <- list()
#i <- "Activated CD8 T cell"#
for (i in set) {
  x <-  as.character(a[,i])
  x <- x[nchar(x)!=0]
  x <-  as.character(x)
  l[[i]] <-x
}
exprSet=as.matrix(exprSet)
#3. 开始进行ssGSEA
ssgsea<- gsva(exprSet, l, method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
write.table(file="ssGSEA.txt",ssgsea,sep="\t",quote=F,col.names=T,row.names=T)
ssgsea.1 <- ssgsea
#for (i in colnames(ssgsea)) {
#i <- colnames(ssgsea)[1]
#ssgsea.1[,i] <- (ssgsea[,i] -min(ssgsea[,i]))/(max(ssgsea[,i] )-min(ssgsea[,i] ))

#}
apply(ssgsea.1[,1:6], 2, range)
ssgsea.2 <- ssgsea.1 %>% t()%>% as.data.frame()
####
pdata2=read.table("sample.txt",header=F,row.names = 1)
pdata=pdata2
#设置参考水平,
group_list <- ifelse(str_detect(pdata$V2, "Control"), "Control",
                     "COPD")
#因子型
group_list = factor(group_list,
                    levels = c("Control","COPD"))
group_list
table(group_list)
library(pheatmap)
sample=read.table("sample.txt",sep="\t",header=T,check.names=F,row.names = 1)
library(stringr)
annotation_col = data.frame(sample)
annotation_col = data.frame(group_list)
rownames(annotation_col)=colnames(exprSet)
colnames(annotation_col)[1]<-"group"
p=pheatmap(
  ssgsea.1,
  border_color = NA,
  cluster_rows = T,cluster_cols = F,
  color = colorRampPalette(colors = c("blue","white","tomato"))(100),
  labels_row = NULL,
  annotation_col = annotation_col,
  clustering_method = "ward.D2",
  fontsize_col = 3,
  cutree_cols = 2,
  show_rownames = T,
  show_colnames = F,
)
####


group=read.table("group.txt",sep="\t",header=T,check.names=F)
y=ssgsea.2

y2=t(ssgsea)
y=y2

data <- cbind(y,group)
data <-data[,c(29,30,1:28)]
data=pivot_longer(data=data,
                  cols = 3:30,
                  names_to = "celltype",
                  values_to = "proportion")
pdf("box4444.pdf", width = 15, height = 8)
ggboxplot(data = data ,
          x ="celltype",#籍形图中的分组变量
          y ="proportion",#会制箱形图的响应变量。
          combine = TRUE,#是否将数据合并为一个箱形图
          merge = FALSE,#是香将相同值的分组合并。
          color ="black",#箱形图边框的颜色。
          fill ="group",#形图填充色
          palette = c("blue","red"),
          title = NULL,#图形标题。
          xlab ="ssGSEA",#x 标
          ylab ="Expression",#y 抽标签
          bxp.errorbar = FALSE,#是否在箱形图中会制误差务
          bxp.errorbar.width = 0.2,#误差务密度。
          facet.by = NULL,#基于哪些变量进行分面
          panel.labs = NULL,#分面的
          short.panel.labs = TRUE,
          linetype ="solid",#线
          size = NULL,#产图形大小奔箱形图的宽度
          widthnotch= FALSE,#产是否在箱形
          outlier.shape = 20,#异常
          select = NULL,#要绘制的变量
          remove = NULL,#产不要绘制的变最。
          order = NULL,#籍形图的排序方式。",#如何会制误差，可以是
          error.plot = "pointrange",
          label = NULL,#产要添加的标签
          font.label = list(size = 12,color ="black"),
          label.select = NULL,#产要添加标签的数据点
          repel = TRUE,
          label.rectangle = TRUE, 
          ggtheme = theme_pubr())+ 
  theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 1)) + 
  stat_compare_means(label = "p.signif",method = "wilcox.test",ref.group= ".all.",hide.ns = F,symnum.args=list(cutpoints = c(0,0.001,0.01,0.05,1),symbols = c("***","**","*","ns")))
dev.off()

sig_gene <- c( "AKR1B10","VPS13D")


library(psych)

x <- t(exprSet)

x <- x[,sig_gene]

y <- t(ssgsea.1)

library(psych)

d <- corr.test(x,y,use="complete",method = 'spearman')

r <- d$r
p <- d$p
write.table(file="REA.txt",r,sep="\t",quote=F,col.names=T,row.names=T)
write.table(file="PVA.txt",p,sep="\t",quote=F,col.names=T,row.names=T)
write.table(file="DATA.txt",data,sep="\t",quote=F,col.names=T,row.names=T)

library(ggcorrplot)
pdf("box4.pdf", width = 15, height = 8)
ggcorrplot(t(d$r), 
           
           show.legend = T, 
           
           digits = 2,  sig.level = 0.05,
           
           insig = 'blank',lab = T)+coord_flip() 

library(pheatmap)

library(reshape2)



if (!is.null(p)){
  
  ssmt <- p< 0.001
  
  p[ssmt] <-'***'
  
  smt <- p >0.001& p < 0.01
  
  p[ssmt] <-'**'
  
  smt <- p >0.01& p <0.05
  
  p[smt] <- '*'
  
  p[!ssmt&!smt]<- ''
  
} else {
  
  p <- F
  
}
mycol<-colorRampPalette(c("blue","white","tomato"))(100)

pheatmap(r,scale = "none",cluster_row = F, cluster_col = F, border=NA,
         
         display_numbers = p,fontsize_number = 12, number_color = "white",
         
         cellwidth = 20, cellheight =20,color=mycol)
dev.off()
#单基因GSEA##
remotes::install_github("loukesio/ltc_palettes")
library(ltc)
setwd("C://Users//COPD//GSEA1") 
library(tidyverse)
load("MLG.rda")
exp3=read.table("1.norexp_GSE128708.txt", header=T, sep="\t", check.names=F,row.names = 1)
rt3=read.table("sample.txt", header=F, sep="\t", check.names=F,row.names = 1)
coldat1 <- exp3["AKR1B10",-c(1:21)] %>% t() %>% as.data.frame()
coldat1$group <- ifelse(coldat1$AKR1B10 < median(coldat1$AKR1B10) ,"low","high")
coldat1 <- coldat1 %>% arrange(AKR1B10)
exp <- exp3[,rownames(coldat1)]
identical(rownames(coldat1),colnames(exp))
group_list <- factor(coldat1$group,levels = c("low","high"))

#先做差异分析##
library(limma)
design <- model.matrix(~group_list)
### 比较矩阵命名
design
### 2.线性模型拟合
fit <- lmFit(exp,design)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
deg <- topTable(fit2, coef = 2, number = Inf)
DEG = na.omit(deg) #differently expressed genes
write.table(DEG, file='DEG1.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
#进行注释 fold change
logFC_cutoff <- 0.585
type1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
write.table(DEG, file='DEG2.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
table(DEG$change)
DEG <- DEG %>% rownames_to_column("Gene")
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
#单基因GSEA##
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
##1.2 按照logFC值对基因进行排序
###1: 提取logFC值，并储存在一个向量中
geneList = DEG[,2]
###2: 对geneList进行命名
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
###3: 根据logFC值降序排列
geneList = sort(geneList, decreasing = TRUE)

msigdb_GMTs <- "msigdb_v7.0_GMTs"
##指定用于GSEA富集分析的gmt文件
msigdb <- "c2.cp.kegg.v7.0.entrez.gmt"
##读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))
C2<-GSEA(geneList,TERM2GENE = kegmt, pvalueCutoff = 0.5) #GSEA分析
#转换成数据框
C2_df <- as.data.frame(C2)
save(C2, file = "C2.Rdata")
load("C2.Rdata")

write.table(C2_df,file="KEGG AKR1B10.txt",sep = "\t",row.names = T,col.names = NA,quote = F)
gseaplot2(C2, geneSetID = c("KEGG_METABOLISM_OF_XENOBIOTICS_BY_CYTOCHROME_P450",
                            "KEGG_STEROID_HORMONE_BIOSYNTHESIS",
                            "KEGG_GLUTATHIONE_METABOLISM",
                            "KEGG_TRYPTOPHAN_METABOLISM",
                            "KEGG_ARACHIDONIC_ACID_METABOLISM"
), subplots = 1:2,color = c("red","yellow","blue","green","purple"))

####2akr1b10####
coldat1 <- exp3["AKR1B10",-c(1:21)] %>% t() %>% as.data.frame()
coldat1$group <- ifelse(coldat1$AKR1B10 < median(coldat1$AKR1B10) ,"low","high")
coldat1 <- coldat1 %>% arrange(AKR1B10)
exp <- exp3[,rownames(coldat1)]
identical(rownames(coldat1),colnames(exp))
group_list <- factor(coldat1$group,levels = c("low","high"))
#先做差异分析##
library(limma)
design <- model.matrix(~group_list)
### 比较矩阵命名
design
### 2.线性模型拟合
fit <- lmFit(exp,design)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
deg <- topTable(fit2, coef = 2, number = Inf)
DEG = na.omit(deg) #differently expressed genes
write.table(DEG, file='DEG1.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
#进行注释 fold change
logFC_cutoff <- 0.585
type1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
write.table(DEG, file='DEG2.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
table(DEG$change)
DEG <- DEG %>% rownames_to_column("Gene")
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
#单基因GSEA##
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
##1.2 按照logFC值对基因进行排序
###1: 提取logFC值，并储存在一个向量中
geneList = DEG[,2]
###2: 对geneList进行命名
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
###3: 根据logFC值降序排列
geneList = sort(geneList, decreasing = TRUE)

msigdb_GMTs <- "msigdb_v7.0_GMTs"
##指定用于GSEA富集分析的gmt文件
msigdb <- "h.all.v2023.2.Hs.entrez.gmt"
##读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))
C2<-GSEA(geneList,TERM2GENE = kegmt, pvalueCutoff = 0.25) #GSEA分析
#转换成数据框
C2_df <- as.data.frame(C2)
save(C2, file = "C2.Rdata")
load("C2.Rdata")

write.table(C2_df,file="HALL AKR1B10-1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

gseaplot2(C2, geneSetID = c("HALLMARK_XENOBIOTIC_METABOLISM",
                            "HALLMARK_FATTY_ACID_METABOLISM",
                            "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                            "HALLMARK_MTORC1_SIGNALING",
                            "HALLMARK_ADIPOGENESIS"
), subplots = 1:2,color = c("red","yellow","blue","green","purple"))
####3akr1b10####
coldat1 <- exp3["AKR1B10",-c(1:21)] %>% t() %>% as.data.frame()

coldat1$group <- ifelse(coldat1$AKR1B10 < median(coldat1$AKR1B10) ,"low","high")
coldat1 <- coldat1 %>% arrange(AKR1B10)

exp <- exp3[,rownames(coldat1)]
identical(rownames(coldat1),colnames(exp))

group_list <- factor(coldat1$group,levels = c("low","high"))

#先做差异分析##
library(limma)
design <- model.matrix(~group_list)
### 比较矩阵命名
design

### 2.线性模型拟合
fit <- lmFit(exp,design)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
deg <- topTable(fit2, coef = 2, number = Inf)

DEG = na.omit(deg) #differently expressed genes
write.table(DEG, file='DEG1.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
#进行注释 fold change

logFC_cutoff <- 0.585
type1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
write.table(DEG, file='DEG2.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
table(DEG$change)
DEG <- DEG %>% rownames_to_column("Gene")
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)
#单基因GSEA##
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
##1.2 按照logFC值对基因进行排序
###1: 提取logFC值，并储存在一个向量中
geneList = DEG[,2]
###2: 对geneList进行命名
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
###3: 根据logFC值降序排列
geneList = sort(geneList, decreasing = TRUE)

msigdb_GMTs <- "msigdb_v7.0_GMTs"
##指定用于GSEA富集分析的gmt文件
msigdb <- "c2.cp.reactome.v7.0.entrez.gmt"
##读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))
C2<-GSEA(geneList,TERM2GENE = kegmt, pvalueCutoff = 0.25) #GSEA分析
#转换成数据框
C2_df <- as.data.frame(C2)
save(C2, file = "C2.Rdata")
load("C2.Rdata")

write.table(C2_df,file="REAT AKR1B10.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

gseaplot2(C2, geneSetID = c("HALLMARK_FATTY_ACID_METABOLISM",
                            "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
                            "HALLMARK_P53_PATHWAY",
                            "HALLMARK_GLYCOLYSIS",
                            "HALLMARK_IL2_STAT5_SIGNALING"
), subplots = 1:2,color = c("red","yellow","blue","green","purple"))




####################################################
coldat1 <- exp3["VPS13D",-c(1:21)] %>% t() %>% as.data.frame()

coldat1$group <- ifelse(coldat1$VPS13D < median(coldat1$VPS13D) ,"low","high")
coldat1 <- coldat1 %>% arrange(VPS13D)

exp <- exp3[,rownames(coldat1)]
identical(rownames(coldat1),colnames(exp))

group_list <- factor(coldat1$group,levels = c("low","high"))

#先做差异分析##
library(limma)
design <- model.matrix(~group_list)
### 比较矩阵命名
design

### 2.线性模型拟合
fit <- lmFit(exp,design)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
deg <- topTable(fit2, coef = 2, number = Inf)

DEG = na.omit(deg) #differently expressed genes
write.table(DEG, file='DEG1.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
#进行注释 fold change

logFC_cutoff <- 0.585
type1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
write.table(DEG, file='DEG2.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
table(DEG$change)

DEG <- DEG %>% rownames_to_column("Gene")

library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)



#单基因GSEA##

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
##1.2 按照logFC值对基因进行排序
###1: 提取logFC值，并储存在一个向量中
geneList = DEG[,2]
###2: 对geneList进行命名
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
###3: 根据logFC值降序排列
geneList = sort(geneList, decreasing = TRUE)

msigdb_GMTs <- "msigdb_v7.0_GMTs"

##指定用于GSEA富集分析的gmt文件
msigdb <- "c2.cp.kegg.v7.0.entrez.gmt"

##读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))
C2<-GSEA(geneList,TERM2GENE = kegmt, pvalueCutoff = 0.25) #GSEA分析
#转换成数据框
C2_df <- as.data.frame(C2)
save(C2, file = "C2.Rdata")
load("C2.Rdata")

write.table(C2_df,file="KEGG VPS13D-1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

gseaplot2(C2, geneSetID = c("KEGG_ALLOGRAFT_REJECTION",
                            
                            "KEGG_PRIMARY_IMMUNODEFICIENCY",
                            "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION",
                            "KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
                            "KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY"
), subplots = 1:2,color = c("red","yellow","blue","green","purple"))

###################VPS13D#####################
coldat1 <- exp3["VPS13D",-c(1:21)] %>% t() %>% as.data.frame()

coldat1$group <- ifelse(coldat1$VPS13D < median(coldat1$VPS13D) ,"low","high")
coldat1 <- coldat1 %>% arrange(VPS13D)

exp <- exp3[,rownames(coldat1)]
identical(rownames(coldat1),colnames(exp))

group_list <- factor(coldat1$group,levels = c("low","high"))

#先做差异分析##
library(limma)
design <- model.matrix(~group_list)
### 比较矩阵命名
design

### 2.线性模型拟合
fit <- lmFit(exp,design)
### 3.贝叶斯检验
fit2 <- eBayes(fit)
deg <- topTable(fit2, coef = 2, number = Inf)

DEG = na.omit(deg) #differently expressed genes
write.table(DEG, file='DEG1.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
#进行注释 fold change

logFC_cutoff <- 0.585
type1 = (DEG$P.Value < 0.05)&(DEG$logFC < -logFC_cutoff)
type2 = (DEG$P.Value < 0.05)&(DEG$logFC > logFC_cutoff)
DEG$change = ifelse(type1,"DOWN",ifelse(type2,"UP","NOT"))
write.table(DEG, file='DEG2.txt',sep = "\t",row.names = T,col.names = NA,quote = F)
table(DEG$change)

DEG <- DEG %>% rownames_to_column("Gene")

library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)



#单基因GSEA##

genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')

DEG <- inner_join(DEG,genelist,by=c("Gene"="SYMBOL"))
##1.2 按照logFC值对基因进行排序
###1: 提取logFC值，并储存在一个向量中
geneList = DEG[,2]
###2: 对geneList进行命名
names(geneList) = as.character(DEG[,'ENTREZID'])
head(geneList)
###3: 根据logFC值降序排列
geneList = sort(geneList, decreasing = TRUE)

msigdb_GMTs <- "msigdb_v7.0_GMTs"

##指定用于GSEA富集分析的gmt文件
msigdb <- "h.all.v2023.2.Hs.entrez.gmt"

##读取上面指定的gmt文件
kegmt <- read.gmt(file.path(msigdb_GMTs,msigdb))
C2<-GSEA(geneList,TERM2GENE = kegmt, pvalueCutoff = 0.25) #GSEA分析
#转换成数据框
C2_df <- as.data.frame(C2)
save(C2, file = "C2.Rdata")
load("C2.Rdata")

write.table(C2_df,file="HALL VPS13D-1.txt",sep = "\t",row.names = T,col.names = NA,quote = F)

gseaplot2(C2, geneSetID = c("HALLMARK_INTERFERON_GAMMA_RESPONSE",
                            "HALLMARK_ALLOGRAFT_REJECTION",
                            "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                            "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"
                            
                            
                            
), subplots = 1:2,color = c("red","yellow","blue","green","purple"))




