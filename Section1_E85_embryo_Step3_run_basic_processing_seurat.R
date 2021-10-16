#############################
#### Section1_E85_embryo ####
#############################

#### The scripts used for analyzing E85 embryo data

###########################################
#### Step3_run_basic_processing_Seurat ####
###########################################


work_path = ""

library(Seurat)
library(dplyr)
library(ggplot2)
source("help_code/help_code.R")

res_doublet_DEG = readRDS(paste0(work_path, "/doublet_cluster_2/res_doubelt_DEG.rds"))
res_doublet_DEG = res_doublet_DEG[res_doublet_DEG$doublet_DEG,]

dat = readRDS(paste0(work_path, "/dat.rds"))
count = dat[["gene_count"]]
fd = dat[["df_gene"]]
pd = dat[["df_cell"]]

pd$doublet_DEG = rownames(pd) %in% as.vector(res_doublet_DEG$sample)
pd = pd[,c("sample", "unmatched_rate", "all_exon", "all_intron", "UMI_count", "gene_count", "doublet_score", "detected_doublets", "doublet_cluster", "doublet_DEG", "RT_group")]
print(sum(pd$sample != colnames(count)))
print(sum(rownames(fd) != rownames(count)))
pd$group = "E85_sci"
pd$day = "E8.5"

#### calculate MT_pct and Ribo_pct per cell, for checking cell quality
gene = fd
gene$gene_short_name = gene$gene_name
MT_gene = as.vector(gene[grep("^mt-",gene$gene_short_name),]$gene_id)
Rpl_gene = as.vector(gene[grep("^Rpl",gene$gene_short_name),]$gene_id)
Mrpl_gene = as.vector(gene[grep("^Mrpl",gene$gene_short_name),]$gene_id)
Rps_gene = as.vector(gene[grep("^Rps",gene$gene_short_name),]$gene_id)
Mrps_gene = as.vector(gene[grep("^Mrps",gene$gene_short_name),]$gene_id)
RIBO_gene = c(Rpl_gene, Mrpl_gene, Rps_gene, Mrps_gene)

pd$MT_pct = 100 * Matrix::colSums(count[gene$gene_id %in% MT_gene, ])/Matrix::colSums(count)
pd$RIBO_pct = 100 * Matrix::colSums(count[gene$gene_id %in% RIBO_gene, ])/Matrix::colSums(count)
pd$EXON_pct = 100 * pd$all_exon / (pd$all_exon + pd$all_intron)
pd$log2_umi = log2(pd$UMI_count)

#### cells labeled as doublets (by Scrublet) or from doublet-derived subclusters were filtered out 
pd = pd[!(pd$detected_doublets | pd$doublet_cluster | pd$doublet_DEG),]

#### first we removed cells with exon% > 85%
#### then we identified the peak of the histgram and then set the UMI cutoff as mean(x) +/- 2*sd(x)
### 154,313 cells were retained

x1 = 10.88034; x2 = 15.26966
pd = pd[pd$log2_umi >= x1 & pd$log2_umi <= x2 & EXON_pct <= 85,]
print(dim(pd))

### For each cell, we only retain protein-coding genes, lincRNA genes and pseudogenes
fd = fd[(fd$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & fd$chr %in% paste0("chr", c(1:19, "M")),]
count = count[rownames(fd), rownames(pd)]
obj = CreateSeuratObject(count, assay = "RNA", meta.data = pd, min.cells=0, min.features=0)
saveRDS(obj, paste0(work_path, "/obj_E8.5b_exp.rds"))

obj_processed = doClusterSeurat(obj)
saveRDS(obj, paste0(work_path, "/obj_E8.5b_processed.rds"))



