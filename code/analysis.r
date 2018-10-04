#Packages
library(DistMap)#for ground truth
library(useful)#for corner
library(nnet)#needed function to compute position of max
#library(devtools)
#install_github("marouenbg/seurat")
#Install my version of seurat where I fixed a few things (l15 rbind in function inside FitGeneK)
library(Seurat)#Loading my version of Seurat
#require(scales)#for hue palette
#library(XLConnect) #for reading xls

#change working directory
setwd("/home/marouen/dreamChallenge/data/singleCellData")

##load data
#Raw data
raw.data = read.table("dge_raw.txt",
                      sep = "\t",
                      row.names = NULL,
                      stringsAsFactors = F,
                      quote = "")

raw.data.genes = raw.data$V1
raw.data$V1 = NULL
print(grep("'",raw.data.genes,value = T,fixed = T))
raw.data.genes = gsub("'","",raw.data.genes,fixed = T)
raw.data = as.matrix(raw.data)
rownames(raw.data) = raw.data.genes

#Normalized data
normalized.data = read.table("dge_normalized.txt",
                             sep = "\t",
                             row.names = NULL,
                             stringsAsFactors = F,
                             quote = "")

normalized.data.genes = normalized.data$row.names
normalized.data$row.names = NULL
print(grep("'",normalized.data.genes,value = T,fixed = T))
normalized.data.genes = gsub("'","",normalized.data.genes,fixed = T)
normalized.data = as.matrix(normalized.data)
rownames(normalized.data) = normalized.data.genes
stopifnot(all(normalized.data.genes == raw.data.genes))

#In situ
setwd('../refDB')
insitu.matrix = read.table("binarized_bdtnp.csv", sep = ",",header = T)
insitu.genes_orig <- colnames(insitu.matrix)
missingGenes = insitu.genes_orig[which(!insitu.genes_orig %in% normalized.data.genes)]
print(missingGenes)
insitu.genes = gsub(".","-",insitu.genes_orig,fixed = T)
insitu.genes = gsub("-spl-","(spl)",insitu.genes,fixed = T)
stopifnot(all(insitu.genes %in% raw.data.genes))
insitu.matrix = as.matrix(insitu.matrix)
colnames(insitu.matrix) = insitu.genes
#Reduce set of genes
res = dim(insitu.matrix)
indGenes=1:res[2]
nGenes=length(indGenes)
insitu.matrix = insitu.matrix[,indGenes]

#Retrieve cell coordinates
geometry = read.csv("geometry.txt",sep = " ")
colnames(geometry) = c("x","y","z")

#Distmap compute ground truth
dm = new("DistMap",
         raw.data=raw.data,
         data=normalized.data,
         insitu.matrix=insitu.matrix,
         geometry=as.matrix(geometry))
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))
dm <- mapCells(dm)
#Mcc-ordered position matrix descending
#posMat is the ground truth matrix
posMat = matrix(0L, nrow = dim(insitu.matrix)[1], ncol = dim(raw.data)[2])
for(i in 1:dim(raw.data)[2]){
  posMat[,i]=order(-dm@mcc.scores[,i])#highest MCC is more likely
}

#Call Seurat 
#Found duplciate gene betaCOP at row index 386, 387 
#but they have different column values
#meaning that they are actually different genes but were mislabled
#1.rename betaCOP
#ind=which(duplicated(rownames(raw.data))==TRUE)
#rownames(raw.data)[ind]="DebugbetaCOP"
#2.eliminate duplicate gene
ind=which(duplicated(rownames(raw.data))==TRUE)
raw.data=raw.data[-ind,]
#confirm that values are differnet
#all(raw.data[,386]==raw.data[,387])
#Fix ambiguous genes for regex later
rownames(raw.data)[which(rownames(raw.data)=="Blimp-1")]="DebugBlimp1"
rownames(raw.data)[which(rownames(raw.data)=="E(spl)m5-HLH")]="DebugEsplm5HLH"
colnames(insitu.matrix)[which(colnames(insitu.matrix)=="Blimp-1")]="DebugBlimp1"
colnames(insitu.matrix)[which(colnames(insitu.matrix)=="E(spl)m5-HLH")]="DebugEsplm5HLH"

#conversion table of points between the actual repersentation
#and a projection compatible with seurat where the first cell
#is in the top left corner
convTable  = matrix(0L, nrow = dim(geometry)[1]+1, ncol = 2)
convTable[,1]=1:3040 #first column is new system
convTable[1:3039,2]=order(-geometry[,3],geometry[,1])#sort by decreasing 3, and increasing 1
convTable[3040,2]=3040

#Call Seurat and some hacks to prevent errors
#meta.data <- data.frame(rep(1, ncol(raw.data)))
#ident = factor(rep(1, ncol(raw.data)))
colnames(raw.data)=colnames(normalized.data)#Cell names
nbt=CreateSeuratObject("seurat",raw.data=raw.data)
#names(nbt@ident) = colnames(normalized.data)

#Scale Data
nbt=NormalizeData(nbt)
nbt=ScaleData(nbt)
nbt=FindVariableGenes(nbt,mean.function = ExpMean,dispersion.function = LogVMR)

#PCA
nbt=RunPCA(nbt, do.print=FALSE, pcs.compute = 40)#because there were more than 20 sig PCs, PCA on variable genes
#my_color_palette <- CustomPalette(low = "red", high = "green")
PCAPlot(nbt, 1, 2, pt.size = 2)#, cols.use=my_color_palette) #+ scale_color_manual(values = my_color_palette)

#Determine significant PCs
nbt=JackStraw(nbt,num.replicate = 200, num.pc = 40, prop.freq=0.025)#should be 1000
JackStrawPlot(nbt,PCs = 1:40)
#=>There are 24 significant PCs p<0.05
sigPCs=24
nbt <- ProjectPCA(nbt, do.print = FALSE,do.center=FALSE)
genes.sig <- PCASigGenes(nbt,pcs.use = 1:sigPCs, pval.cut = 1e-2, use.full = FALSE)#this should be TRUE but documentation says it is ok

#tSNE
nbt=RunTSNE(nbt,dims.use = 1:sigPCs,max_iter=2000)
TSNEPlot(nbt,pt.size = 1)

#Clustering
nbt=DBClustDimension(nbt,1,2,reduction.use = "tsne",G.use = 2.2,set.ident = TRUE)#DBScan
#nbt=KClustDimension(nbt,1,2,reduction.use = "tsne",set.ident = TRUE, k.use = 10)#KNN
nClusters=length(unique(factor(nbt@ident)))#1st cluster is un-assigned cells
TSNEPlot(nbt,pt.size = 1)

#Plot cells against each other
par(mfrow=c(2,2))
CellPlot(nbt,nbt@cell.names[1],nbt@cell.names[2],do.ident = FALSE)
CellPlot(nbt,nbt@cell.names[3],nbt@cell.names[4],do.ident = FALSE)

#Add in situ genes
#I project the spheroid into rectangle
#I added 1 empty location to make 3039+1=76*40 or 80*38 or 95*32 (Leading dimension is 76= ncols)
insitu.matrix=rbind(insitu.matrix,integer(dim(insitu.matrix)[2]))
insitu.genes=colnames(insitu.matrix)
#reorder bins
insitu.matrix=insitu.matrix[convTable[,2],]
nbt@spatial@insitu.matrix =as.data.frame(insitu.matrix)
lasso.genes.use=unique(c(genes.sig,nbt@var.genes))

#Impute gene expression
nbt <- AddImputedScore(nbt, genes.use=lasso.genes.use,genes.fit=insitu.genes, do.print=FALSE, s.use=40, gram=FALSE)

#Continue with Seurat mapping
#i="rho" i="Blimp-1"
for(i in rev(insitu.genes)){
  nbt=FitGeneK(nbt,i,do.plot=FALSE,do.k = 2,start.pct=mean(nbt@spatial@insitu.matrix[,i]),num.iter = 1)
} 

#show an example mixture model
#par(mfrow=c(2,2))
#nbt_temp=FitGeneK(nbt,"zfh1",do.plot=TRUE,do.k = 2,start.pct=mean(nbt@spatial@insitu.matrix[,"zfh1"]))
#nbt_temp=FitGeneK(nbt,"zen2",do.plot=TRUE,do.k = 2,start.pct=mean(nbt@spatial@insitu.matrix[,"zen2"]))
#nbt_temp=FitGeneK(nbt,"twi",do.plot=TRUE,do.k = 2,start.pct=mean(nbt@spatial@insitu.matrix[,"twi"]))
#nbt_temp=FitGeneK(nbt,"tsh",do.plot=TRUE,do.k = 2,start.pct=mean(nbt@spatial@insitu.matrix[,"tsh"]))

#Map cells
nbt <- InitialMapping(nbt)

#Refine mapping
num.pc=3; num.genes=6
genes.use=PCTopGenes(nbt,pc.use = 1:num.pc,num.genes = num.genes,use.full = TRUE,do.balanced = TRUE)

#impute values for those genes
new.imputed=genes.use[!genes.use%in%rownames(nbt@imputed)]
lasso.genes.use=unique(c(nbt@var.genes,PCASigGenes(nbt,pcs.use = c(1,2,3), pval.cut = 1e-2, use.full = FALSE)))
#use.full should be TRUE but it throws error, the documentation says FALSE is ok
nbt <- AddImputedScore(nbt, genes.use=lasso.genes.use,genes.fit=new.imputed, do.print=FALSE, s.use=40, gram=FALSE)

#Actual refinment
nbt <- RefinedMapping(nbt,genes.use)

#identify positions
corner(nbt@spatial@final.prob)

#calculate centroids for each cell (we do not need this step because here 1bin=1cell)
nbt.centroids=GetCentroids(nbt); colnames(nbt.centroids)=c("bin.tier","bin.dv")
corner(nbt.centroids)

#final mapping and convert back to original coordinates
mapMat = matrix(0L, nrow = dim(nbt@spatial@final.prob)[1], ncol = dim(raw.data)[2])
for(i in 1:dim(raw.data)[2]){
  mapMat[,i]=convTable[1,order(nbt@spatial@final.prob[,i])]#p values are from the smallest
}

#remove 3040

#write accuracy test


