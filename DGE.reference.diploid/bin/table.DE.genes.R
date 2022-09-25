.libPaths("/nfs/users/rg/bborsari/software/R-3.5.2/library")


#************
# LIBRARIES *
#************

library(ggplot2)
library(dplyr)
library(scales)

#********
# BEGIN *
#********

# 1. set wd
setwd("/no_backup/rg/bborsari/projects/ENCODE/EN-TEx/DE.analysis/tasks/DESeq2_reference_diploid")

#---------------------------------------------------------------------------------
# 2. list of DE genes - 4 indivs separately, tissues as replicates (std analyis) |
#---------------------------------------------------------------------------------

## 2.1. upregulated genes (aka upregulated in the reference)
upreg <- read.delim("entex.4indivs/upreg.tsv", h=F, sep="\t", stringsAsFactors = F)
upreg.ids <- read.delim("entex.4indivs/upreg.geneIds.tsv", h=F, sep="\t", stringsAsFactors = F)
upreg.ids$V2 <- ifelse(upreg.ids$V2 == "", NA, upreg.ids$V2)
upreg <- merge(upreg, upreg.ids, all.x = T, by = "V1")
upreg$V3 <- "upregulated"

## 2.2. downregulated genes (aka downregulated in the diploid)
downreg <- read.delim("entex.4indivs/downreg.tsv", h=F, sep="\t", stringsAsFactors = F)
downreg.ids <- read.delim("entex.4indivs/downreg.geneIds.tsv", h=F, sep="\t", stringsAsFactors = F)
downreg.ids$V2 <- ifelse(downreg.ids$V2 == "", NA, downreg.ids$V2)
downreg <- merge(downreg, downreg.ids, all.x = T, by = "V1")
downreg$V3 <- "downregulated"

DEG <- rbind(upreg, downreg)

# write.table(DEG, 
#             file = "~/public_html/paper_ENTEx/entex.webpage.files/S1/table.DE.genes.tsv",
#             row.names = F, col.names = F, sep="\t", quote=F)




#---------------------------------------------------------------------------------
# 3. list of DE genes - 4 indivs separately, tissues as replicates (std analyis) |
#---------------------------------------------------------------------------------

## 3.1. upregulated genes (aka upregulated in the reference)
techreps.upreg <- read.delim("techReps/upreg.tsv", h=F, sep="\t", stringsAsFactors = F)
techreps.upreg.ids <- read.delim("techReps/upreg.geneIds.tsv", h=F, sep="\t", stringsAsFactors = F)
techreps.upreg.ids$V2 <- ifelse(techreps.upreg.ids$V2 == "", NA, techreps.upreg.ids$V2)
techreps.upreg <- merge(techreps.upreg, techreps.upreg.ids, all.x = T, by = "V1")
techreps.upreg$V3 <- "upregulated"
techreps.upreg$V3 <- ifelse(techreps.upreg$V1 %in% upreg$V1, techreps.upreg$V3, paste0(techreps.upreg$V3, "*"))

## 3.2. downregulated genes (aka downregulated in the diploid)
techreps.downreg <- read.delim("techReps/downreg.tsv", h=F, sep="\t", stringsAsFactors = F)
techreps.downreg.ids <- read.delim("techReps/downreg.geneIds.tsv", h=F, sep="\t", stringsAsFactors = F)
techreps.downreg.ids$V2 <- ifelse(techreps.downreg.ids$V2 == "", NA, techreps.downreg.ids$V2)
techreps.downreg <- merge(techreps.downreg, techreps.downreg.ids, all.x = T, by = "V1")
techreps.downreg$V3 <- "downregulated"
techreps.downreg$V3 <- ifelse(techreps.downreg$V1 %in% downreg$V1, techreps.downreg$V3, paste0(techreps.downreg$V3, "*"))


techreps.DEG <- rbind(techreps.upreg, techreps.downreg)

# write.table(techreps.DEG, 
#             file = "~/public_html/paper_ENTEx/entex.webpage.files/S1/table.DE.genes.techReps.liver.tsv",
#             row.names = F, col.names = F, sep="\t", quote=F)




#---------------------------------------------------------------------------------
# 4. list of DE genes - 4 indivs separately, tissues as replicates (std analyis) |
#---------------------------------------------------------------------------------

## 4.1. upregulated genes (aka upregulated in the reference)
GM12878.upreg <- read.delim("GM12878/upreg.tsv", h=F, sep="\t", stringsAsFactors = F)
GM12878.upreg.ids <- read.delim("GM12878/upreg.geneIds.tsv", h=F, sep="\t", stringsAsFactors = F)
GM12878.upreg.ids$V2 <- ifelse(GM12878.upreg.ids$V2 == "", NA, GM12878.upreg.ids$V2)
GM12878.upreg <- merge(GM12878.upreg, GM12878.upreg.ids, all.x = T, by = "V1")
GM12878.upreg$V3 <- "upregulated"
GM12878.upreg$V3 <- ifelse(((GM12878.upreg$V1 %in% upreg$V1) | (GM12878.upreg$V1 %in% techreps.upreg$V1)), 
                           GM12878.upreg$V3, 
                           paste0(GM12878.upreg$V3, "*"))

## 4.2. downregulated genes (aka downregulated in the diploid)
GM12878.downreg <- read.delim("GM12878/downreg.tsv", h=F, sep="\t", stringsAsFactors = F)
GM12878.downreg.ids <- read.delim("GM12878/downreg.geneIds.tsv", h=F, sep="\t", stringsAsFactors = F)
GM12878.downreg.ids$V2 <- ifelse(GM12878.downreg.ids$V2 == "", NA, GM12878.downreg.ids$V2)
GM12878.downreg <- merge(GM12878.downreg, GM12878.downreg.ids, all.x = T, by = "V1")
GM12878.downreg$V3 <- "downregulated"
GM12878.downreg$V3 <- ifelse(((GM12878.downreg$V1 %in% downreg$V1) | (GM12878.downreg$V1 %in% techreps.downreg$V1)), 
                             GM12878.downreg$V3, 
                             paste0(GM12878.downreg$V3, "*"))


GM12878.DEG <- rbind(GM12878.upreg, GM12878.downreg)

# write.table(GM12878.DEG, 
#             file = "~/public_html/paper_ENTEx/entex.webpage.files/S1/table.DE.genes.GM12878.tsv",
#             row.names = F, col.names = F, sep="\t", quote=F)
