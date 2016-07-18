library(biomaRt)


mart <- useMart(biomart = "ensembl", dataset = "scerevisiae_gene_ensembl")
results <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart)
head(results)
head(listAttributes(mart),n=200)
?getBM
length(unique(results$ensembl_gene_id))
length(unique(results$ensembl_transcript_id)
)
a=c("a","a","b")
length(unique(a))
