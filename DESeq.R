


#library(tximport)
oDir=getwd()
sfDir="/media/sven/data/hayo/counts/"
setwd(sfDir)

files=list.files(pattern = "*.sf",path=)
#tximport(files=files,type = "salmon")

#files=/media/sven/data/hayo/counts/





All <- lapply(files,function(i){
  read.csv(i, header=TRUE, sep="\t")
})
counts=round(data.frame(cbind(All[[1]]$NumReads,All[[2]]$NumReads,All[[3]]$NumReads,All[[4]]$NumReads)))
head(counts)
rownames(counts)=All[[1]]$Name
colnames(counts)=files

library(DESeq2)
coldata=data.frame(sampleName=factor(colnames(counts)),treatment=c("med21","med21","wt","wt"))
coldata
dds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = coldata,
                                 design = ~ treatment)
vsd <- vst(dds)
plotPCA(vsd, "treatment")



dds <- DESeq(dds)
res <- results(dds)
summary(res)
res.05 <- results(dds, alpha=.05)
head(res.05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)
resLFC1.05=resLFC1[which(resLFC1$padj<0.1),]
resLFC1.05Order=resLFC1.05[order(resLFC1.05$padj),]

plotMA(res, ylim=c(-5,5))
head(resLFC1.05Order)

resLFC1.05OrderbaseMean=resLFC1.05[order(resLFC1.05$baseMean),]
resLFC1.05OrderbaseMean
