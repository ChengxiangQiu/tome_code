#############################
#### Section5_pseudobulk ####
#############################

#### Inferring the molecular histories of individual cell types

################################
#### Step2_identifying_DEGs ####
################################

#### identifying gene list used for ordering

library(monocle3)
library(Seurat)
library(dplyr)
work_path = ""

#### read individual datasets 
count = readRDS(paste0(work_path, "/embryo_1.rds"))

pd = readRDS(paste0(work_path, "/pd_1.rds"))
pd = data.frame(unique(pd[,c("embryo_id", "embryo_sex", "project", "day", "group")]))
rownames(pd) = as.vector(pd$embryo_id)
pd$day = gsub("E","",pd$day)
pd$day = as.vector(as.numeric(pd$day))
pd = pd[colnames(count),]

fd = data.frame(gene_ID = rownames(count))
fd = fd %>% left_join(mouse_gene, by = "gene_ID") %>% as.data.frame()
rownames(fd) = as.vector(fd$gene_ID)

count1 = count
pd1 = pd


count = readRDS(paste0(work_path, "/embryo_2.rds"))

pd = readRDS(paste0(work_path, "/pd_2.rds"))
pd = data.frame(unique(pd[,c("embryo_id", "embryo_sex", "project", "day", "group")]))
rownames(pd) = as.vector(pd$embryo_id)
pd$day = gsub("E","",pd$day)
pd$day[pd$day == "8.5a"] = 8.5
pd$day = as.vector(as.numeric(pd$day))
pd = pd[colnames(count),]

fd = data.frame(gene_ID = rownames(count))
fd = fd %>% left_join(mouse_gene, by = "gene_ID") %>% as.data.frame()
rownames(fd) = as.vector(fd$gene_ID)

count2 = count
pd2 = pd

### 
count = readRDS(paste0(work_path, "/embryo_3.rds"))

pd = readRDS(paste0(work_path, "/pd_3.rds"))
pd = data.frame(unique(pd[,c("embryo_id", "embryo_sex", "project", "day", "group")]))
rownames(pd) = as.vector(pd$embryo_id)
pd$day = gsub("E","",pd$day)
pd$day[pd$day == "8.5b"] = 8.5
pd$day = as.vector(as.numeric(pd$day))
pd = pd[colnames(count),]

fd = data.frame(gene_ID = rownames(count))
fd = fd %>% left_join(mouse_gene, by = "gene_ID") %>% as.data.frame()
rownames(fd) = as.vector(fd$gene_ID)

count3 = count
pd3 = pd

### 
count = readRDS(paste0(work_path, "/embryo_4.rds"))

pd = readRDS(paste0(work_path, "/pd_4.rds"))
pd = data.frame(unique(pd[,c("embryo_id", "embryo_sex", "project", "day", "group")]))
rownames(pd) = as.vector(pd$embryo_id)
pd$day = gsub("E","",pd$day)
pd$day = as.vector(as.numeric(pd$day))
pd = pd[colnames(count),]

fd = data.frame(gene_ID = rownames(count))
fd = fd %>% left_join(mouse_gene, by = "gene_ID") %>% as.data.frame()
rownames(fd) = as.vector(fd$gene_ID)

count4 = count
pd4 = pd



#### first, identifying DEGs across timepoints of Pijuan's
gene_overlap = intersect(rownames(count2), rownames(count4))
count2 = count2[gene_overlap,]

obj = CreateSeuratObject(count2, meta.data = pd2)

count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])
fd = data.frame(gene_ID = rownames(count))
fd = fd %>% left_join(mouse_gene, by = "gene_ID") %>% as.data.frame()
rownames(fd) = as.vector(fd$gene_ID)

cds = new_cell_data_set(count,
                        cell_metadata = pd,
                        gene_metadata = fd)

gene_fits <- fit_models(cds, model_formula_str = "~day")
fit_coefs <- coefficient_table(gene_fits)
res = fit_coefs %>% filter(term != "(Intercept)") %>% arrange(q_value) %>%
    select(gene_ID, gene_short_name, term, p_value, q_value, estimate)
saveRDS(res, paste0(work_path, "/DEG_2.rds"))


### then, identifying DEGs across timepoints of new E8.5 + deeper sequencing
### and adding project as covariant
gene_overlap = intersect(rownames(count2), rownames(count4))
count3 = count3[gene_overlap,]
count4 = count4[gene_overlap,]

obj3 = CreateSeuratObject(count3, meta.data = pd3)
obj4 = CreateSeuratObject(count4, meta.data = pd4)
obj = merge(obj3, obj4)

count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])
fd = data.frame(gene_ID = rownames(count))
fd = fd %>% left_join(mouse_gene, by = "gene_ID") %>% as.data.frame()
rownames(fd) = as.vector(fd$gene_ID)

cds = new_cell_data_set(count,
                        cell_metadata = pd,
                        gene_metadata = fd)

gene_fits <- fit_models(cds, model_formula_str = "~day + project")
fit_coefs <- coefficient_table(gene_fits)
res = fit_coefs %>% filter(term != "(Intercept)") %>% arrange(q_value) %>%
    select(gene_ID, gene_short_name, term, p_value, q_value, estimate)
saveRDS(res, paste0(work_path, "/DEG_3_4_project.rds"))


### finally, identifying DEGs differentially expressed between cells vs. nuclei
gene_overlap = intersect(rownames(count2), rownames(count4))
count2 = count2[gene_overlap,]
count3 = count3[gene_overlap,]
count4 = count4[gene_overlap,]

obj2 = CreateSeuratObject(count2, meta.data = pd2)
obj2$tech = "10x"
obj3 = CreateSeuratObject(count3, meta.data = pd3)
obj3$tech = "sci"
obj4 = CreateSeuratObject(count4, meta.data = pd4)
obj4$tech = "sci"
obj = merge(obj2, list(obj3, obj4))

count = GetAssayData(obj, slot = "counts")
pd = data.frame(obj[[]])
fd = data.frame(gene_ID = rownames(count))
fd = fd %>% left_join(mouse_gene, by = "gene_ID") %>% as.data.frame()
rownames(fd) = as.vector(fd$gene_ID)

cds = new_cell_data_set(count,
                        cell_metadata = pd,
                        gene_metadata = fd)

gene_fits <- fit_models(cds, model_formula_str = "~tech")
fit_coefs <- coefficient_table(gene_fits)
res = fit_coefs %>% filter(term != "(Intercept)") %>% arrange(q_value) %>%
    select(gene_ID, gene_short_name, term, p_value, q_value, estimate)
saveRDS(res, paste0(work_path, "/DEG_2_3_4_noday.rds"))



###################### make DEGs used for ordering #######################

DEG_2 = readRDS(paste0(work_path, "/DEG_2.rds")); table(DEG_2$term)
DEG_2_3_4_noday = readRDS(paste0(work_path, "/DEG_2_3_4_noday.rds")); table(DEG_2_3_4_noday$term)
DEG_3_4_project = readRDS(paste0(work_path, "/DEG_3_4_project.rds")); table(DEG_3_4_project$term)

### Finally, I use this gene list for ordering embryos/samples
### top 3,000 genes associated with timepoints from Pijuan-Sala's data
### plus top 3,000 genes associated with timepoints from Beth's and Jun's data (using project as covariant)
### excluding significant genes which are differentially expressed between Pijuan-Sala's and Beth's + Jun's data (only using technology as independent variant)
### 534 genes are remained
genes_ordering = union(DEG_2$gene_ID[1:3000], 
                       DEG_3_4_project$gene_ID[DEG_3_4_project$term == "day"][1:3000])
exclude_gene = DEG_2_3_4_noday$gene_ID[DEG_2_3_4_noday$p_value < 0.05 & DEG_2_3_4_noday$term == "techsci"]
genes_ordering = genes_ordering[!genes_ordering %in% as.vector(exclude_gene)]
print(length(genes_ordering))

saveRDS(genes_ordering, paste0(work_path, "/genes_ordering.rds"))


