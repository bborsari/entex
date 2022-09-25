.libPaths("/nfs/users/rg/bborsari/software/R-3.5.2/library")


# this script performs DE analysis 
# between personal and ref genome mappings
# using 5 RNA-seq experiments
# with available replicates
# for GM12878



#************
# LIBRARIES *
#************

library(DESeq2)
library(ggplot2)
library(reshape2)
library(edgeR)
library(data.table)
suppressPackageStartupMessages(library(rtracklayer))
options(stringsAsFactors =F)


#********
# BEGIN *
#********

# 0. set working directory
setwd("/no_backup/rg/bborsari/projects/ENCODE/EN-TEx/DE.analysis/tasks/DESeq2_reference_diploid/GM12878/")


# 1. read diploid and ref expression matrices
gm.diploid <- fread('GM12878.diploid.expected_count.tsv', h=T, data.table = F)
rownames(gm.diploid) <- gm.diploid$gene_id
gm.diploid$gene_id <- NULL

gm.reference <- fread('GM12878.expected_count.tsv', h=T, data.table = F)
rownames(gm.reference) <- gm.reference$gene_id
gm.reference$gene_id <- NULL

stopifnot(identical(rownames(gm.diploid), rownames(gm.reference)))


# 2. set NAs to zero
gm.reference[is.na(gm.reference)] <- 0 
gm.diploid[is.na(gm.diploid)] <- 0


# 3. import gencode v19 annotation
genes <- import('/users/rg/projects/references/Annotation/H.sapiens/gencode19/gencode.v19.annotation.gtf')


# 4. import gene types
# this file contains classification of ensembl types into
# long vs. short genes 
# + protein-coding vs. non-coding 
gene_types <- read.table('/no_backup/rg/bborsari/projects/ENCODE/EN-TEx/DE.analysis/super_gene_type.tsv', sep="\t",h=F)


# 6. retrieve long genes from gencode v19
long_genes <- unique(as.data.frame(genes[genes$gene_type %in% gene_types[gene_types$V2=='long','V1'],])[,'gene_id'])


# 7. retrieve mt genes 
mtGenes <- as.data.frame(genes[seqnames(genes) %in% c('chrM') & genes$type=='gene', ])[,'gene_id']


# 8. merge reference and diploid matrices 
# into a unique matrix
gm <- cbind(gm.reference[rownames(gm.diploid),], gm.diploid)


# 9. keep only long genes + remove chrM genes
gm <- gm[long_genes,]
gm <- gm[grep('ENSG', rownames(gm)),]
gm <- gm[-which(rownames(gm) %in% mtGenes),]


# 10. set parameters for DE analysis
minCov=1
fc=2
pvalue=0.05
pseudo=1


# 11. DeSeq2 analysis

## 11.1. filter genes for minimum coverage
gm <- gm[rowMeans(gm) > minCov, ]

# 11.2. create "c" metadata dataframe including batch info
gm.c <- data.frame(type = rep(c("reference", "diploid"), each=10))
rownames(gm.c) <- c(colnames(gm.reference), colnames(gm.diploid))
gm.c$batch <- gsub("\\_.*","", rownames(gm.c))

gm.c$type <- factor(gm.c$type)
gm.c$batch <- factor(gm.c$batch)

# 11.3. perform DE analysis
gm.dds <- DESeqDataSetFromMatrix( round(gm),
                                  colData = gm.c,
                                  design= ~batch + type )
gm.dds <- DESeq(gm.dds)
resultsNames(gm.dds) # lists the coefficients
res <- results(gm.dds, name = "type_reference_vs_diploid")
write.table(res, file = paste0("res.GM12878.batch.tsv"), sep="\t", quote=F)

