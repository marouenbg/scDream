#Packages
library(DistMap)

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

#Geometry
geometry = read.csv("geometry.txt",sep = " ")
colnames(geometry) = c("x","y","z")

##Run DistMap
dm = new("DistMap",
         raw.data=raw.data,
         data=normalized.data,
         insitu.matrix=insitu.matrix,
         geometry=as.matrix(geometry))
dm <- binarizeSingleCellData(dm, seq(0.15, 0.5, 0.01))

#write.table(dm@binarized.data,file = "binarizedData_distMap.csv",sep = ",",row.names = T,col.names = T)
dm <- mapCells(dm)

#find non unique max
count=0
dimMcc=dim(dm@mcc.scores)
for (i in 1:dimMcc[2]){
  if (sort(dm@mcc.scores[,i])[dimMcc[1]]==sort(dm@mcc.scores[,i])[dimMcc[1]-2]){
    count=count+1
  }
}

#Questions
#what is DV
#Why there arent one max MCC



