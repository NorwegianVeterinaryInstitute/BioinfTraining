## To install cummeRbund and DESeq2 (do it once)
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

getwd()
setwd('../Desktop/course_data/DESeq2/')

library('DESeq2')

data <- read.delim('../featureCounts/count_gene', skip=1, sep="\t")
dim(data)
head(data)
colnames(data)

count <- data[7:12]

colnames(count) <- c('Con1_Rep1', 'Con1_Rep2', 'Con1_Rep3', 'Con2_Rep1', 'Con2_Rep2', 'Con2_Rep3')
rownames(count) <- data$Geneid

condition <- data.frame(c('con1', 'con1', 'con1', 'con2', 'con2', 'con2'))
colnames(condition) <- 'group'
rownames(condition) <- colnames(count)

dds <- DESeqDataSetFromMatrix(count, condition, design = ~group)
dds <- DESeq(dds)
res <- results(dds)
summary(res)

summary(res, alpha=0.05)
res_05 <- subset(res, padj <= 0.05)

write.table(res_05, 'DESeq2.txt', quote=F, row.names=F, sep='\t')

plotMA(res)
