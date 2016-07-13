library("airway")
dir <- system.file("extdata", package="airway", mustWork=TRUE)
list.files(dir)
csvfile <- file.path(dir,"sample_table.csv")
sampleTable <- read.csv(csvfile,row.names=1)
sampleTable
filenames <- file.path(dir, paste0(sampleTable$Run, "_subset.bam"))
file.exists(filenames)
library("Rsamtools")
bamfiles <- BamFileList(filenames, yieldSize=2000000)

library("GenomicFeatures")
gtffile <- file.path(dir,"Homo_sapiens.GRCh37.75_subset.gtf")
txdb <- makeTxDbFromGFF(gtffile, format="gtf")
ebg <- exonsBy(txdb, by="gene")

library("GenomicAlignments")
library("BiocParallel")

register(MulticoreParam(2))
se <- summarizeOverlaps(features=ebg, 
                        reads=bamfiles,
                        mode="Union",
                        singleEnd=FALSE,
                        ignore.strand=TRUE,
                        fragments=TRUE )
colData(se) <- DataFrame(sampleTable)

data("airway")
se <- airway

se$dex <- relevel(se$dex, "untrt")
se$dex


data("airway")
se <- airway
se <- se[ rowSums(assay(se)) >= 5, ]
se$dex <- relevel(se$dex, "untrt")

library("DESeq2")
dds <- DESeqDataSet(se, design = ~ cell + dex)
countdata <- assay(se)
vsd <- vst(dds)
plotPCA(vsd, "dex")
dds <- DESeq(dds)
res <- results(dds)
table(res)

resLFC1 <- results(dds, lfcThreshold=1)
res.05 <- results(dds, alpha=.05)
table(res.05$padj < .05)
table(resLFC1$padj < 0.1)

resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)


library("AnnotationDbi")
library("Homo.sapiens")
columns(Homo.sapiens)

res$symbol <- mapIds(Homo.sapiens,
                     keys=row.names(res),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
y$genes$symbol <- res$symbol

res$entrez <- mapIds(Homo.sapiens,
                     keys=row.names(res),
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
y$genes$entrez <- res$entrez

res$symbol <- mapIds(Homo.sapiens,
                     keys=row.names(res),
                     column="GENENAME",
                     keytype="ENSEMBL",
                     multiVals="first")
y$genes$symbol <- res$symbol

resOrdered <- res[order(res$padj),]
head(resOrdered)

resOrderedDF <- as.data.frame(resOrdered)[seq_len(100),]
write.csv(resOrderedDF, file="results.csv")

library("Glimma")
glMDPlot(lrt, 
         counts=y$counts, 
         anno=y$genes, 
         groups=y$samples$dex, 
         samples=colnames(y),
         status=tt.all$table$FDR < 0.1,
         id.column="gene.id")
