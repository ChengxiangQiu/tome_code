#############################
#### Section5_pseudobulk ####
#############################

#### Inferring the molecular histories of individual cell types

########################
#### Step3_ordering ####
########################

library(monocle)
library(Seurat)
library(dplyr)
work_path = ""
source("help_code/help_code.R")

count = readRDS(paste0(work_path, "/embryo_1.rds"))
pd = readRDS(paste0(work_path, "/pd_1.rds"))
pd = data.frame(unique(pd[,c("embryo_id", "embryo_sex", "project", "day", "group")]))
rownames(pd) = as.vector(pd$embryo_id)
pd = pd[colnames(count),]
count1 = count; pd1 = pd

count = readRDS(paste0(work_path, "/embryo_2.rds"))
pd = readRDS(paste0(work_path, "/pd_2.rds"))
pd = data.frame(unique(pd[,c("embryo_id", "embryo_sex", "project", "day", "group")]))
rownames(pd) = as.vector(pd$embryo_id)
pd = pd[colnames(count),]
count2 = count; pd2 = pd

count = readRDS(paste0(work_path, "/embryo_3.rds"))
pd = readRDS(paste0(work_path, "/pd_3.rds"))
pd = data.frame(unique(pd[,c("embryo_id", "embryo_sex", "project", "day", "group")]))
rownames(pd) = as.vector(pd$embryo_id)
pd = pd[colnames(count),]
count3 = count; pd3 = pd

count = readRDS(paste0(work_path, "/embryo_4.rds"))
pd = readRDS(paste0(work_path, "/pd_4.rds"))
pd = data.frame(unique(pd[,c("embryo_id", "embryo_sex", "project", "day", "group")]))
rownames(pd) = as.vector(pd$embryo_id)
pd = pd[colnames(count),]
count4 = count; pd4 = pd

### Of note, here we used the gene list by intersecting the Pijuan-Sala's and Beth/Jun's data

gene_overlap = intersect(rownames(count2), rownames(count4))
count1 = count1[rownames(count1) %in% gene_overlap,]
count2 = count2[rownames(count2) %in% gene_overlap,]
count3 = count3[rownames(count3) %in% gene_overlap,]
count4 = count4[rownames(count4) %in% gene_overlap,]

obj1 = CreateSeuratObject(count1, meta.data = pd1)
obj2 = CreateSeuratObject(count2, meta.data = pd2)
obj3 = CreateSeuratObject(count3, meta.data = pd3)
obj4 = CreateSeuratObject(count4, meta.data = pd4)

obj = merge(obj1, list(obj2, obj3, obj4))

count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])
fd = data.frame(gene_ID = rownames(count))
fd = fd %>% left_join(mouse_gene, by = "gene_ID") %>% as.data.frame()
rownames(fd) = as.vector(fd$gene_ID)

pd = new("AnnotatedDataFrame", data = pd)
fd = new("AnnotatedDataFrame", data = fd)
cds = newCellDataSet(count, pd, fd)

cds = estimateSizeFactors(cds)
cds = estimateDispersions(cds)

saveRDS(cds, paste0(work_path, "/cds_embryo.rds"))

#### using monocle/2 for ordering embryos/samples

library(monocle)
library(dplyr)
work_path = ""
cds = readRDS(paste0(work_path, "/cds_embryo.rds"))
genes_ordering = readRDS(paste0(work_path, "/genes_ordering.rds"))

cds = setOrderingFilter(cds, genes_ordering)

cds = reduceDimension(cds, max_components = 2, method = "DDRTree")

cds = orderCells(cds)

cds$day = factor(cds$day, levels = paste0("E", c(3.5, 4.5, 5.25, 5.5, 6.25, seq(6.5, 8.25, 0.25), "8.5a", "8.5b", seq(9.5, 13.5, 1))))
plot_cell_trajectory(cds, color_by = "day")

GM_state <- function(cds){
    if (length(unique(pData(cds)$State)) > 1){
        T0_counts <- table(pData(cds)$State, pData(cds)$day)[,"E6.25"]
        return(as.numeric(names(T0_counts)[which
                                           (T0_counts == max(T0_counts))]))
    } else {
        return (1)
    }
}
cds <- orderCells(cds, root_state = GM_state(cds))
plot_cell_trajectory(cds, color_by = "Pseudotime")
head(pData(cds))
saveRDS(data.frame(pData(cds)), paste0(work_path, "/cds_embryo_pseudotime.rds"))

