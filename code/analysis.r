#Packages
library(DistMap)#for ground truth
library(useful)#for corner
library(nnet)#needed function to compute position of max
library(devtools)
install_github("marouenbg/seurat")
#Install my version of seurat where I fixed a few things (l15 rbind in function inside FitGeneK)
library(Seurat)#Loading my version of Seurat
#require(scales)#for hue palette
#library(XLConnect) #for reading xls

accuracyTest <- function(posMat,geometry,nbt,convTable,raw.data,numberGenes,mapping,varMethod){
  #final mapping and convert back to original coordinates
  mapMat = matrix(0L, nrow = dim(nbt@spatial@final.prob)[1], ncol = dim(raw.data)[2])
  for(i in 1:dim(raw.data)[2]){
    mapMat[,i]=convTable[order(nbt@spatial@final.prob[,i]),1]#p values are from the smallest
  }
  
  #remove 3040 vv
  mapMat = mapMat[-3040,]
  
  #write accuracy test
  #1. counting if max MCC cell is in top 10
  count=0
  for(i in 1:dim(posMat)[2]){
    if (posMat[1,i] %in% mapMat[1:10,i]){
      count=count+1
    } 
  }
  perc=count/dim(posMat)[2]#percentage of cells identified corrcetly in top 10
  #2. average distance between cells and the max MCC cell
  avDistance=0
  for(i in 1:dim(posMat)[2]){
    avDistance=avDistance+dist(rbind(geometry[posMat[1,i],], geometry[mapMat[1,i],]))
  }
  avDistance=avDistance/dim(posMat)[2]
  
  #write results to csv file
  write.csv(t(mapMat)[,1:10], paste("submission",numberGenes,mapping,varMethod,".csv",sep=""), row.names = 1:dim(raw.data)[2], sep=",")
  return(list(perc,avDistance))
}

#Parameters
mapping='projection'
binary='seurat'

#change working directory
setwd("../data/singleCellData")

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
setwd("../refDB")
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
# and binary DistMap as well
if (binary=='genDistmap'){
  rownames(dm@binarized.data)[which(rownames(dm@binarized.data)=="Blimp-1")]="DebugBlimp1"
  rownames(dm@binarized.data)[which(rownames(dm@binarized.data)=="E(spl)m5-HLH")]="DebugEsplm5HLH"
}
  
#conversion table of points between the actual repersentation
#and a projection compatible with seurat where the first cell
#is in the top left corner
convTable  = matrix(0L, nrow = dim(geometry)[1]+1, ncol = 2)
convTable[,1]=1:3040 #first column is new system
if(mapping=='original'){
    convTable[,2]=1:3040
}else if (mapping=='projection'){
    convTable[1:3039,2]=order(-geometry[,3],geometry[,1])#sort by decreasing 3, and increasing 1
    convTable[3040,2]=3040
}

#Call Seurat and some hacks to prevent errors
colnames(raw.data)=colnames(normalized.data)#Cell names
nbt=CreateSeuratObject("seurat",raw.data=raw.data)

#Scale Data
nbt=NormalizeData(nbt)
nbt=ScaleData(nbt)
nbt=FindVariableGenes(nbt,mean.function = ExpMean,dispersion.function = LogVMR, y.cutoff=1, x.low.cutoff = 0.1,
x.high.cutoff = 8)

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

insitu.genes.save=insitu.genes
insitu.matrix.save=insitu.matrix
hvg.info.save=nbt@hvg.info
for(numberGenes in c(20,40,60)){
  for(varMethod in c("dispersion","geneMean")){
    insitu.genes=insitu.genes.save
    insitu.matrix=insitu.matrix.save
    nbt@hvg.info=hvg.info.save
    #select genes based on variability
    #1-dispersion/ 2 gene expression mean
    if(varMethod=='dispersion'){
      #the dispersion is default order
    }else if(varMethod=='geneMean'){
      nbt@hvg.info=nbt@hvg.info[order(-nbt@hvg.info$gene.mean),]
    }
    #intersect and keep order of priority
    insitu.genes=intersect(rownames(nbt@hvg.info),insitu.genes)
    if(numberGenes==20){
      insitu.genes=insitu.genes[1:20]
    }else if(numberGenes==40){
      insitu.genes=insitu.genes[1:40]
    }else if(numberGenes==60){
      insitu.genes=insitu.genes[1:60]
    }
    insitu.matrix=insitu.matrix[,insitu.genes]
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
    
    #use binarization of DistMap instead of mixture models
    if (binary=='genDistmap'){
      dmBinary = t(dm@binarized.data)+1
      write.table(dmBinary, file = "/home/marouen/seurat/R/dmBinary", sep = " ")
    }
    
    #Continue with Seurat mapping
    #i="rho" i="Blimp-1"
    for(i in rev(insitu.genes)){
      nbt=FitGeneK(nbt,i,do.plot=FALSE,do.k = 2,start.pct=mean(nbt@spatial@insitu.matrix[,i]),num.iter = 1)
    } 
    
    #Map cells
    mapping='Initial'
    nbt <- InitialMapping(nbt)
    resultsBeforeRefine=accuracyTest(posMat,geometry,nbt,convTable,raw.data,numberGenes,mapping,varMethod)
    print(resultsBeforeRefine)
    
    #Refine mapping
    num.pc=2; num.genes=4
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
    
    mapping='refined'
    resultsAfterRefine=accuracyTest(posMat,geometry,nbt,convTable,raw.data,numberGenes,mapping,varMethod)
    print(resultsAfterRefine)
  }
}
