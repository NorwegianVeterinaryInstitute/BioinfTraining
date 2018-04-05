## To install cummeRbund and DESeq2 (do it once)
# source("https://bioconductor.org/biocLite.R")
# biocLite("cummeRbund")

# Get working directory
getwd()

# change working directory
setwd('../Desktop/course_data/cummeRbund/')
getwd()

# Load library cummeRbund
library('cummeRbund')

# Load cuffdiff data
cuff_short <- readCufflinks('../cuffdiff_short_output')
cuff_short

cuff_long <- readCufflinks('../cuffdiff_long_output')
cuff_long

disp_short_genes <- dispersionPlot(genes(cuff_short))
disp_short_isoforms <- dispersionPlot(isoforms(cuff_short))

disp_short_genes
disp_short_isoforms

scv_short_genes <- fpkmSCVPlot(genes(cuff_short))
scv_short_isoforms <- fpkmSCVPlot(isoforms(cuff_short))

scv_short_genes
scv_short_isoforms

dendro_short_genes <- csDendro(genes(cuff_short))
dendro_short_genes_rep <- csDendro(genes(cuff_short), replicates=T)

dendro_short_isoforms <- csDendro(isoforms(cuff_short))
dendro_short_isoforms_rep <- csDendro(isoforms(cuff_short), replicates=T)

MDSplot(genes(cuff_short))
plotPCA(genes(cuff_short))

MDSplot(genes(cuff_short), replicates=T)


#Significant number of genes
# sigMatrix(genes(cuff_short)) - This will not work
# alpha is q_value value and 0.05 is default in cummeRbund
sigMatrix(cuff_short, level='genes')
sigMatrix(cuff_short, level='isoforms')


diffData_short_genes <- diffData(genes(cuff_short))
dim(diffData_short_genes)
head(diffData_short_genes)

fpkm_short_genes <- fpkmMatrix(genes(cuff_short))
head(fpkm_short_genes)

fpkm_short_genes_rep <- repFpkmMatrix(genes(cuff_short))
head(fpkm_short_genes_rep)

# Subset for significantly DE genes and isoforms
DE_short_genes_sig <- subset(diffData_short_genes, significant == "yes")
DE_short_genes_sig <- subset(diffData_short_genes, q_value <= 0.05)

dim(DE_short_genes_sig)
head(DE_short_genes_sig)

# Write significantly DE genes output to file
write.table(DE_short_genes_sig, "DE_short_genes_sig.txt", quote=F, row.names=F, sep='\t')


myGeneIds <- DE_short_genes_sig$gene_id
myGene <- getGenes(cuff_short, myGeneIds)
myGene

csHeatmap(myGene)
csHeatmap(myGene, cluster="both", replicates=T)



myGeneIds <- DE_short_genes_sig$gene_id
myGene <- getGenes(cuff_short, myGeneIds[1:5])
myGene

csHeatmap(myGene)
expressionBarplot(myGene)
expressionBarplot(myGene, replicates=T)

myGeneIds <- "FBgn0000064"
myGene <- getGenes(cuff_short, myGeneIds)
myGene
expressionBarplot(isoforms(myGene))

expressionPlot(isoforms(myGene))


# Subset again for up- and down- regulated genes
DE_short_genes_sig_UP <- subset(DE_short_genes_sig, log2_fold_change > 0)
DE_short_genes_sig_DOWN <- subset(DE_short_genes_sig, log2_fold_change < 0)
