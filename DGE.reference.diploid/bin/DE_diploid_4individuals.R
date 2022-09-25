.libPaths("/nfs/users/rg/bborsari/software/R-3.5.2/library")

# Original script by Anna Vlasova
# Modified by Beatrice Borsari


# this script performs DE analysis 
# between personal and ref genome mappings
# for each of the 4 entex individuals
# using tissues as replicates


#************
# LIBRARIES *
#************

library(DESeq2)
library(ggplot2)
library(reshape2)
library(edgeR)
suppressPackageStartupMessages(library(rtracklayer))
options(stringsAsFactors =F)


#********
# BEGIN *
#********

setwd("/no_backup/rg/bborsari/projects/ENCODE/EN-TEx/DE.analysis")
outfolder <- "/no_backup/rg/bborsari/projects/ENCODE/EN-TEx/DE.analysis/tasks/DESeq2_reference_diploid/entex.4indivs"


# 1. read diploid and ref expression matrices
m.diploid <- read.csv('EN-TEx/EnTex.genes.97samples.diploid.expected_count.tsv.gz', sep="\t",h=T)
rownames(m.diploid) <- m.diploid$gene_id
m.diploid$gene_id <- NULL

m.reference <- read.csv('EN-TEx/EnTex.genes.97samples.reference.expected_count.tsv.gz', sep="\t",h=T)
colnames(m.reference) <- c(colnames(m.reference)[-1], 'last')
m.reference <- m.reference[,-98]


# 2. set NAs to zero
m.reference[is.na(m.reference)] <- 0 
m.diploid[is.na(m.diploid)] <- 0


# 3. read metadata
meta <- read.csv('EN-TEx/metadata.97samples.reference.diploid.tsv', sep="\t",h=T, comment.char = "")
rownames(meta) <- meta$labExpId


# 4. import gencode v24 annotation
genes <- import('/users/rg/bborsari/public_html/paper_ENTEx/YY/input/gencode.v24.primary_assembly.annotation.gtf')


# 5. import gene types
# this file contains classification of ensembl types into
# long vs. short genes 
# + protein-coding vs. non-coding 
gene_types <- read.table('super_gene_type.tsv', sep="\t",h=F)


# 6. retrieve long genes from gencode v24
long_genes <- unique(as.data.frame(genes[genes$gene_type %in% gene_types[gene_types$V2=='long','V1'],])[,'gene_id'])


# 7. retrieve mt genes 
mtGenes <- as.data.frame(genes[seqnames(genes) %in% c('chrM') & genes$type=='gene', ])[,'gene_id']


# 8. merge reference and diploid matrices 
# into a unique matrix
m <- cbind(m.reference[rownames(m.diploid),], m.diploid)
m_spikes <- m[grep('ERC', rownames(m)),]


# 9. keep only long genes + remove chrM genes
m <- m[long_genes,]
m <- m[grep('ENSG', rownames(m)),]
m <- m[-which(rownames(m) %in% mtGenes),]


# 10. set parameters for DE analysis
minCov=1
fc=2
pvalue=0.05
pseudo=1

# 11. pre-processing spike-ins
N <- colSums(m)
m_spikes <- m_spikes[,colnames(m)]
#TMM normalized, default for calcNormFactors
spikeFactor <- calcNormFactors(m_spikes, lib.size=N)
names(spikeFactor) <- colnames(m_spikes)


# 12. DE analysis: all refs vs. all diploid
indList = c("ENC002", "ENC003", "ENC001", "ENC004")
genesDE_perInd = lapply(indList, function(ind){
  print(ind)
  m_ind <- m[,meta[meta$donorId==ind,'labExpId']]
  m_ind <- m_ind[rowMeans(m_ind)>minCov,]
  c <- meta[colnames(m_ind), c('type'), drop=FALSE]
  rownames(c) <- colnames(m_ind)
  c$type <- factor(c$type)
  
  dds <- DESeqDataSetFromMatrix( round(m_ind),
                                 colData = c,
                                 design= ~type)
  
  sizeFactors(dds) <- spikeFactor[colnames(m_ind)]
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  write.table(res, file = paste0(outfolder, "/", ind, "/res.", ind,".tsv"),
              sep="\t", quote=F)
  
  save(dds, file=paste0(outfolder, "/", ind, paste0("dds.DESeq2.reference_vs_diploid.", ind,".Rdata")))
  
})

