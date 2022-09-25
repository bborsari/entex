#------------------------
# DEG analysis of       |
# each indiv. using     |
# tissues as replicates | 
#------------------------

# 1. DEG analysis 
# ~/software/R-3.5.2/bin/Rscript bin/DE_diploid_4individuals.R

# 2. get gene names of DE genes 
# cat entex.4indivs/ENC00*/res.ENC00* | awk '$3>1 && $NF < 0.1' | cut -f1 | sort -u |cut -d "." -f1 | ~bborsari/software/R-3.4.1/bin/Rscript ~bborsari/scripts/R/biomaRt.R -a "c(\"ensembl_gene_id\", \"hgnc_symbol\")" -f "c(\"ensembl_gene_id\")" -i stdin -o entex.4indivs/upreg.geneIds.tsv
# cat entex.4indivs/ENC00*/res.ENC00* | awk '$3<-1 && $NF < 0.1' | cut -f1 | sort -u |cut -d "." -f1 | ~bborsari/software/R-3.4.1/bin/Rscript ~bborsari/scripts/R/biomaRt.R -a "c(\"ensembl_gene_id\", \"hgnc_symbol\")" -f "c(\"ensembl_gene_id\")" -i stdin -o entex.4indivs/downreg.geneIds.tsv

# 3. get list of all DE genes
# cat entex.4indivs/ENC00*/res.ENC00* | awk '$3>1 && $NF < 0.1' | cut -f1 | sort -u |cut -d "." -f1 > entex.4indivs/upreg.tsv
# cat entex.4indivs/ENC00*/res.ENC00* | awk '$3<-1 && $NF < 0.1' | cut -f1 | sort -u |cut -d "." -f1 > entex.4indivs/downreg.tsv

# 4. list of genes for entex portal obtained with:
# bin/table.DE.genes.R


#-------------------------------
# DEG analyis of 7 experiments |
# 6 exps. w/ tech reps         |
# 1 exp. duplicated            |
#-------------------------------

# 1. DEG analysis 
# ~/software/R-3.5.2/bin/Rscript bin/DE_diploid_4individuals.techReps.R

# 2. get gene names of DE genes
# cat techReps/enc00*/res.tsv | awk '$3>1 && $NF < 0.1' | cut -f1 | sort -u |cut -d "." -f1 | ~bborsari/software/R-3.4.1/bin/Rscript ~bborsari/scripts/R/biomaRt.R -a "c(\"ensembl_gene_id\", \"hgnc_symbol\")" -f "c(\"ensembl_gene_id\")" -i stdin -o techReps/upreg.geneIds.tsv
# cat techReps/enc00*/res.tsv | awk '$3<-1 && $NF < 0.1' | cut -f1 | sort -u |cut -d "." -f1 | ~bborsari/software/R-3.4.1/bin/Rscript ~bborsari/scripts/R/biomaRt.R -a "c(\"ensembl_gene_id\", \"hgnc_symbol\")" -f "c(\"ensembl_gene_id\")" -i stdin -o techReps/downreg.geneIds.tsv

# 3. get list of all DE genes
# cat techReps/enc00*/res.tsv | awk '$3>1 && $NF < 0.1' | cut -f1 | sort -u |cut -d "." -f1 > techReps/upreg.tsv
# cat techReps/enc00*/res.tsv | awk '$3<-1 && $NF < 0.1' | cut -f1 | sort -u |cut -d "." -f1 > techReps/downreg.tsv



#--------------------------
# DEG analysis of GM12878 |
#--------------------------

# 1. DEG analysis
# ~/software/R-3.5.2/bin/Rscript bin/DE_diploid_GM12878.R

# 2. get gene names of DE genes
# awk '$3>1 && $NF < 0.1' GM12878/res.GM12878.batch.tsv | cut -f1 | sort -u |cut -d "." -f1 | ~bborsari/software/R-3.4.1/bin/Rscript ~bborsari/scripts/R/biomaRt.R -a "c(\"ensembl_gene_id\", \"hgnc_symbol\")" -f "c(\"ensembl_gene_id\")" -i stdin -o GM12878/upreg.geneIds.tsv
# awk '$3<-1 && $NF < 0.1' GM12878/res.GM12878.batch.tsv | cut -f1 | sort -u |cut -d "." -f1 | ~bborsari/software/R-3.4.1/bin/Rscript ~bborsari/scripts/R/biomaRt.R -a "c(\"ensembl_gene_id\", \"hgnc_symbol\")" -f "c(\"ensembl_gene_id\")" -i stdin -o GM12878/downreg.geneIds.tsv

# 3. get list of all DE genes
awk '$3>1 && $NF < 0.1' GM12878/res.GM12878.batch.tsv | cut -f1 | sort -u |cut -d "." -f1 > GM12878/upreg.tsv
awk '$3<-1 && $NF < 0.1' GM12878/res.GM12878.batch.tsv | cut -f1 | sort -u |cut -d "." -f1 > GM12878/downreg.tsv
