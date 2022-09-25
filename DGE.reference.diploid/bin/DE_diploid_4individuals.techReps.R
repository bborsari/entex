.libPaths("/nfs/users/rg/bborsari/software/R-3.5.2/library")


# this script performs DE analysis between personal and ref genome mappings
# for 7 RNA-seq experiments 
# with available replicates


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
outfolder <- "/no_backup/rg/bborsari/projects/ENCODE/EN-TEx/DE.analysis/tasks/DESeq2_reference_diploid/techReps/"


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


# 11. list experiments with available replicates

loe <- list()
loe[[1]] <- c("enc001-ENCSR276MMH-adrenal_gland", "ENCSR276MMH_1", "ENCSR276MMH_2", "ENCSR276MMH_1_d", "ENCSR276MMH_2_d")
loe[[2]] <- c("enc001-ENCSR429EWK-thoracic_aorta", "ENCSR429EWK_1", "ENCSR429EWK_2", "ENCSR429EWK_1_d",	"ENCSR429EWK_2_d")
loe[[3]] <- c("enc001-ENCSR532LJV-thyroid_gland",	"ENCSR532LJV_1", "ENCSR532LJV_2", "ENCSR532LJV_1_d", "ENCSR532LJV_2_d")
loe[[4]] <- c("enc001-ENCSR853BNH-gastrocnemius_medialis", "ENCSR853BNH_1", "ENCSR853BNH_2", "ENCSR853BNH_1_d",	"ENCSR853BNH_2_d")
loe[[5]] <- c("enc002-ENCSR023ZXN-thyroid_gland", "ENCSR023ZXN_1_d",	"ENCSR023ZXN_2_d",	"ENCSR023ZXN_1",	"ENCSR023ZXN_2")
loe[[6]] <- c("enc002-ENCSR954PZB-adrenal_gland",	"ENCSR954PZB_1", "ENCSR954PZB_2", "ENCSR954PZB_1_d", "ENCSR954PZB_2_d")
loe[[7]] <- c("enc003-ENCSR504QMK-ENCSR226KML-liver",	"ENCSR504QMK_1", "ENCSR226KML_1", "ENCSR504QMK_1_d", "ENCSR226KML_1_d")



# 12. perform DE for each replicated experiment
for (i in 1:length(loe)) {
  
  print(i)
  
  m.sub <- m[, colnames(m) %in% loe[[i]][2:5]]
  m.sub <- m.sub[, loe[[i]][2:5]]
  
  
  N <- colSums(m.sub)
  m.sub_spikes <- m_spikes[,colnames(m.sub)]
  #TMM normalized, default for calcNormFactors
  spikeFactor <- calcNormFactors(m.sub_spikes, lib.size=N)
  names(spikeFactor) <- colnames(m.sub_spikes)
  
  m.sub <- m.sub[rowMeans(m.sub)>minCov,]
  c <- meta[colnames(m.sub), c('type'), drop=FALSE]
  rownames(c) <- colnames(m.sub)
  c$type <- factor(c$type)
  
  dds <- DESeqDataSetFromMatrix(round(m.sub),
                                colData = c,
                                design= ~type)
  
  sizeFactors(dds) <- spikeFactor[colnames(m.sub)]
  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)
  res <- results(dds)
  
  id <- loe[[i]][1]
  
  write.table(res, file = paste0(outfolder, "/", id, "/res.tsv" ),
              sep="\t", quote=F)
  
  save(dds, file=paste0(outfolder, "/", id, "/dds.DESeq2.reference_vs_diploid.Rdata"))
  
}

