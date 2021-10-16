#############################
#### Section5_pseudobulk ####
#############################

#### Inferring the molecular histories of individual cell types

###################################
#### Step4_specific_trajectory ####
###################################

#### aggregating cells for each cell state for each individual embryo, excluding ExE tissues ####

library(Seurat)
library(dplyr)
library(monocle3)
work_path = ""
exclude_celltype = c("Embryonic visceral endoderm", "Extraembryonic visceral endoderm", "Parietal endoderm", "Extraembryonic ectoderm")

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])

pd = readRDS(paste0(work_path, "/pd_", kk, ".rds"))
pd_sub = pd[!pd$celltype %in% exclude_celltype,]
pd_sub$state_id = paste(pd_sub$day, pd_sub$celltype, pd_sub$embryo_id, sep = "#")
count = readRDS(paste0(work_path, "/count_", kk, ".rds"))

if(sum(rownames(pd) != colnames(count)) != 0){
    print(XXX)
} else {
    state_list = as.vector(unique(pd_sub$state_id))
    exp = NULL
    for(i in 1:length(state_list)){
        print(paste0(i,"/",length(state_list)))
        if(sum(pd_sub$state_id == state_list[i]) == 1){
            exp_tmp = count[,colnames(count) %in% as.vector(rownames(pd_sub)[pd_sub$state_id == state_list[i]])]
            exp = cbind(exp, exp_tmp)
        } else {
            exp_tmp = count[,colnames(count) %in% as.vector(rownames(pd_sub)[pd_sub$state_id == state_list[i]])]
            exp = cbind(exp, Matrix::rowSums(exp_tmp))
        }
    }
    colnames(exp) = state_list
    rownames(exp) = rownames(count)
}

saveRDS(exp, paste0(work_path, "/embryo_state_", kk, ".rds"))


#### step 2, merging datasets 

library(Seurat)
library(dplyr)
library(monocle3)
work_path = ""

count = readRDS(paste0(work_path, "/embryo_state_1.rds"))
tmp = readRDS(paste0(work_path, "/pd_1.rds"))
pd = data.frame(sample_id = colnames(count),
                 day = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][1])),
                 celltype = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][2])),
                 embryo_id = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][3])))
tmp = data.frame(unique(tmp[,c("embryo_id", "embryo_sex", "project", "group")]))
pd = pd %>% left_join(tmp, by = "embryo_id") %>% as.data.frame()
rownames(pd) = as.vector(pd$sample_id)
pd = pd[colnames(count),]
count1 = count; pd1 = pd

count = readRDS(paste0(work_path, "/embryo_state_2.rds"))
tmp = readRDS(paste0(work_path, "/pd_2.rds"))
pd = data.frame(sample_id = colnames(count),
                day = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][1])),
                celltype = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][2])),
                embryo_id = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][3])))
tmp = data.frame(unique(tmp[,c("embryo_id", "embryo_sex", "project", "group")]))
pd = pd %>% left_join(tmp, by = "embryo_id") %>% as.data.frame()
rownames(pd) = as.vector(pd$sample_id)
pd = pd[colnames(count),]
count2 = count; pd2 = pd

count = readRDS(paste0(work_path, "/embryo_state_3.rds"))
tmp = readRDS(paste0(work_path, "/pd_3.rds"))
pd = data.frame(sample_id = colnames(count),
                day = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][1])),
                celltype = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][2])),
                embryo_id = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][3])))
tmp = data.frame(unique(tmp[,c("embryo_id", "embryo_sex", "project", "group")]))
pd = pd %>% left_join(tmp, by = "embryo_id") %>% as.data.frame()
rownames(pd) = as.vector(pd$sample_id)
pd = pd[colnames(count),]
count3 = count; pd3 = pd

count = readRDS(paste0(work_path, "/embryo_state_4.rds"))
tmp = readRDS(paste0(work_path, "/pd_4.rds"))
pd = data.frame(sample_id = colnames(count),
                day = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][1])),
                celltype = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][2])),
                embryo_id = unlist(lapply(colnames(count), function(x) strsplit(x,"[#]")[[1]][3])))
tmp = data.frame(unique(tmp[,c("embryo_id", "embryo_sex", "project", "group")]))
pd = pd %>% left_join(tmp, by = "embryo_id") %>% as.data.frame()
rownames(pd) = as.vector(pd$sample_id)
pd = pd[colnames(count),]
count4 = count; pd4 = pd

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
obj$Anno = paste0(obj$day, ":", obj$celltype)

saveRDS(obj, paste0(work_path, "/obj_embryo_state.rds"))


#### four interested cellular trajectories ###


library(Seurat)
library(ggplot2)
library(gplots)
library(dplyr)
library(randomcoloR)
source("help_code/help_code.R")
work_path = ""

calculateDEG <- function(x, time, batch){
    x_var = apply(x, 1, sd)
    x = x[x_var > quantile(x_var, 0.9),]
    res = NULL
    x = as.matrix(x)
    for(j in 1:nrow(x)){
        print(paste0(j, "/", nrow(x)))
        fit = lm(x[j,] ~ time + batch)
        res = rbind(res, data.frame(esti = summary(fit)$coefficients[2,1], pval = summary(fit)$coefficients[2,4], esti2 = summary(fit)$coefficients[3,1], pval2 = summary(fit)$coefficients[3,4]))
    }
    rownames(res) = rownames(x)
    return(res)
}

### Neural crest (PNS glia)
neural_crest_state = c(paste0("E", seq(6.25, 7.25, 0.25), ":Epiblast"),
                       paste0("E", seq(7.5, 8, 0.25), ":Rostral neuroectoderm"),
                       paste0("E", c(8.25, "8.5a", "8.5b"), ":Neural crest"),
                       paste0("E", seq(9.5, 13.5, 1), ":Neural crest (PNS glia)"))

### Otic epithelium
otic_epithelium_state = c(paste0("E", c(6.25), ":Epiblast"),
            paste0("E", seq(6.5, 7.25, 0.25), ":Primitive streak and adjacent ectoderm"),
            paste0("E", seq(7.5, 8.25, 0.25), ":Surface ectoderm"),
            paste0("E", c("8.5a", "8.5b"), ":Placodal area"),
            paste0("E", seq(9.5, 13.5, 1), ":Otic epithelium"))

### Midgut/Hindgut epithelium
midgut_hindgut_state = c(paste0("E", c(6.25), ":Epiblast"),
            paste0("E", seq(6.5, 6.75, 0.25), ":Primitive streak and adjacent ectoderm"),
            paste0("E", c(7), ":Anterior primitive streak"),
            paste0("E", c(7.25), ":Definitive endoderm"),
            paste0("E", c(seq(7.5, 8.25, 0.25), "8.5a", "8.5b"), ":Gut"),
            paste0("E", seq(9.5, 10.5, 1), ":Gut and lung epithelium"),
            paste0("E", seq(11.5, 13.5, 1), ":Midgut/Hindgut epithelium"))

### Cardiomyocytes
cardiomyocytes_state = c(paste0("E", c(6.25), ":Epiblast"),
            paste0("E", c(6.5), ":Primitive streak and adjacent ectoderm"),
            paste0("E", seq(6.75, 7.25, 0.25), ":Nascent mesoderm"),
            paste0("E", c(7.5), ":Splanchnic mesoderm"),
            paste0("E", c(seq(7.75, 8.25, 0.25), "8.5a", "8.5b", 9.5), ":First heart field"),
            paste0("E", seq(10.5, 13.5, 1), ":Cardiomyocytes"))

state_list = c(neural_crest_state, otic_epithelium_state, midgut_hindgut_state, cardiomyocytes_state)
celltype_list = as.vector(unique(unlist(lapply(as.vector(state_list), function(x) strsplit(x,"[:]")[[1]][2]))))

set.seed(123)
n <- length(celltype_list)
cols <- distinctColorPalette(n)
names(cols) = celltype_list
pie(rep(1, n), col=cols)

p = data.frame(x = sample(1:3, n, replace = T),
            y = sample(1:3, n, replace = T),
            celltype = celltype_list) %>%
    ggplot(aes(x, y)) + 
    geom_point(aes(color = celltype), size = 5) + 
    scale_colour_manual(values = cols) +
    theme_classic(base_size = 15)

cols = readRDS(paste0(work_path, "/scatter_legend.rds"))
                
obj = readRDS(paste0(work_path, "/obj_embryo_state.rds"))
embryo = readRDS(paste0(work_path, "/cds_embryo_pseudotime.rds"))
pd = data.frame(obj[[]]) %>% left_join(embryo[,c("embryo_id", "Pseudotime")], by = "embryo_id")
obj$Pseudotime = as.vector(pd$Pseudotime)

vline_x = max(obj$Pseudotime[obj$project == "pijuan"])

### Neural crest
obj_sub = subset(obj, subset = Anno %in% neural_crest_state)
table(obj_sub$celltype)

count = GetAssayData(object = obj_sub, slot = "counts")
count = t(t(count) / Matrix::colSums(count)) * 100000
count = log(count + 1)

res = calculateDEG(count, obj_sub$Pseudotime, factor(obj_sub$project))
res = res %>% 
    mutate(gene_ID = rownames(res)) %>%
    left_join(mouse_gene[,c("gene_ID", "gene_short_name")], by = "gene_ID") %>%
    mutate(qval = p.adjust(pval, method = "bon")) %>%
    mutate(qval2 = p.adjust(pval2, method = "bon"))
hist(res$esti, 100)
hist(res$esti2, 100)
up1 = mean(res$esti) + 2*sd(res$esti)
down1 = mean(res$esti) - 2*sd(res$esti)
res_sub1 = res[(res$esti < down1 |
                  res$esti > up1) &
                  res$qval < 0.05,]
up2 = mean(res$esti2) + 2*sd(res$esti2)
down2 = mean(res$esti2) - 2*sd(res$esti2)
res_sub2 = res[(res$esti2 < down2 |
                   res$esti2 > up2) &
                   res$qval2 < 0.05,]
res_sub = res[res$gene_short_name %in% as.vector(res_sub1$gene_short_name) & !res$gene_short_name %in% as.vector(res_sub2$gene_short_name),]
dim(res_sub) # 115
res_sub %>% filter(esti > 0) %>% arrange(qval) %>% head()
res_sub %>% filter(esti < 0) %>% arrange(qval) %>% head()

target_gene_list = c("Sox10","Pax3","Ncam1","Lin28b")
df = data.frame(obj_sub[[]])
for (i in 1:length(target_gene_list)){
    gene_i = mouse_gene$gene_ID[mouse_gene$gene_short_name == target_gene_list[i]]
    p = df %>%
        mutate(exp = as.vector(count[gene_i,])) %>%
        ggplot(aes(Pseudotime, exp)) + 
        geom_point(aes(color = celltype)) + 
        geom_smooth() + 
        scale_colour_manual(values = cols) +
        theme_classic(base_size = 15) +
        theme(legend.position="none") +
        geom_vline(xintercept = vline_x, color = "blue") +
        labs(title = target_gene_list[i])
    assign(paste0("p", i), p)
}

ggsave(paste0(work_path, "/neural_crest_pseudobulk_1_add.pdf"), p1 + p2 + p3 + p4, width=10, height=10, dpi = 300, units = "in", useDingbats=FALSE) 

exp = count[as.vector(res_sub$gene_ID), order(obj_sub$Pseudotime)]
rownames(exp) = res_sub$gene_short_name
pdf(paste0(work_path, "/neural_crest_pseudobulk_2_add.pdf"), 10, 10)
g = heatmap.2(as.matrix(exp), col=viridis, scale="row", Rowv = TRUE, Colv = FALSE, key=T, density.info="none", trace="none", cexRow=0.6, cexCol = 0.6)
dev.off()
gene_list = rownames(exp)[rev(g$rowInd)]
write.table(gene_list, paste0(work_path, "/neural_crest_pseudobulk_gene_list_dd.txt"), col.names=F, row.names=F, quote=F, sep="\t")

state_list = colnames(exp)[(g$colInd)]
df = df[state_list,]
df_sub = df[df$project == "pijuan",]
vline_xxx = c(1:nrow(df))[df$embryo_id == df_sub$embryo_id[nrow(df_sub)]]
p = ggplot(df, aes(1:nrow(df),Pseudotime)) + geom_point(aes(color = celltype)) + 
    geom_vline(xintercept = vline_xxx, color = "blue") +
    scale_colour_manual(values = cols) +
    theme_classic(base_size = 15) + theme(legend.position="none")
ggsave(paste0(work_path, "/neural_crest_pseudobulk_3_add.pdf"), p, width=8, height=4, dpi = 300, units = "in", useDingbats=FALSE) 





### Otic epithelium 
obj_sub = subset(obj, subset = Anno %in% otic_epithelium_state)
table(obj_sub$celltype)

count = GetAssayData(object = obj_sub, slot = "counts")
count = t(t(count) / Matrix::colSums(count)) * 100000
count = log(count + 1)

res = calculateDEG(count, obj_sub$Pseudotime, factor(obj_sub$project))
res = res %>% 
    mutate(gene_ID = rownames(res)) %>%
    left_join(mouse_gene[,c("gene_ID", "gene_short_name")], by = "gene_ID") %>%
    mutate(qval = p.adjust(pval, method = "bon")) %>%
    mutate(qval2 = p.adjust(pval2, method = "bon"))
hist(res$esti, 100)
hist(res$esti2, 100)
up1 = mean(res$esti) + 2*sd(res$esti)
down1 = mean(res$esti) - 2*sd(res$esti)
res_sub1 = res[(res$esti < down1 |
                    res$esti > up1) &
                   res$qval < 0.05,]
up2 = mean(res$esti2) + 2*sd(res$esti2)
down2 = mean(res$esti2) - 2*sd(res$esti2)
res_sub2 = res[(res$esti2 < down2 |
                    res$esti2 > up2) &
                   res$qval2 < 0.05,]
res_sub = res[res$gene_short_name %in% as.vector(res_sub1$gene_short_name) & !res$gene_short_name %in% as.vector(res_sub2$gene_short_name),]
dim(res_sub) # 122
res_sub %>% filter(esti > 0) %>% arrange(qval) %>% head()
res_sub %>% filter(esti < 0) %>% arrange(qval) %>% head()

target_gene_list = c("Pax2","Lmx1a","Rbms3","Car4")
df = data.frame(obj_sub[[]])
for (i in 1:length(target_gene_list)){
    gene_i = mouse_gene$gene_ID[mouse_gene$gene_short_name == target_gene_list[i]]
    p = df %>%
        mutate(exp = as.vector(count[gene_i,])) %>%
        ggplot(aes(Pseudotime, exp)) + 
        geom_point(aes(color = celltype)) + 
        geom_smooth() + 
        scale_colour_manual(values = cols) +
        theme_classic(base_size = 15) +
        theme(legend.position="none") +
        geom_vline(xintercept = vline_x, color = "blue") +
        labs(title = target_gene_list[i])
    assign(paste0("p", i), p)
}

ggsave(paste0(work_path, "/otic_epithelium_pseudobulk_1_add.pdf"), p1 + p2 + p3 + p4, width=10, height=10, dpi = 300, units = "in", useDingbats=FALSE) 

exp = count[as.vector(res_sub$gene_ID), order(obj_sub$Pseudotime)]
rownames(exp) = res_sub$gene_short_name
pdf(paste0(work_path, "/otic_epithelium_pseudobulk_2_add.pdf"), 10, 10)
g = heatmap.2(as.matrix(exp), col=viridis, scale="row", Rowv = TRUE, Colv = FALSE, key=T, density.info="none", trace="none", cexRow=0.6, cexCol = 0.6)
dev.off()
gene_list = rownames(exp)[rev(g$rowInd)]
write.table(gene_list, paste0(work_path, "/otic_epithelium_pseudobulk_gene_list_add.txt"), col.names=F, row.names=F, quote=F, sep="\t")

state_list = colnames(exp)[(g$colInd)]
df = df[state_list,]
df_sub = df[df$project == "pijuan",]
vline_xxx = c(1:nrow(df))[df$embryo_id == df_sub$embryo_id[nrow(df_sub)]]
p = ggplot(df, aes(1:nrow(df),Pseudotime)) + geom_point(aes(color = celltype)) + 
    geom_vline(xintercept = vline_xxx, color = "blue") +
    scale_colour_manual(values = cols) +
    theme_classic(base_size = 15) + theme(legend.position="none")
ggsave(paste0(work_path, "/otic_epithelium_pseudobulk_3_add.pdf"), p, width=8, height=4, dpi = 300, units = "in", useDingbats=FALSE) 






### Midgut/Hindgut epithelium 
obj_sub = subset(obj, subset = Anno %in% midgut_hindgut_state)
table(obj_sub$celltype)

count = GetAssayData(object = obj_sub, slot = "counts")
count = t(t(count) / Matrix::colSums(count)) * 100000
count = log(count + 1)

res = calculateDEG(count, obj_sub$Pseudotime, factor(obj_sub$project))
res = res %>% 
    mutate(gene_ID = rownames(res)) %>%
    left_join(mouse_gene[,c("gene_ID", "gene_short_name")], by = "gene_ID") %>%
    mutate(qval = p.adjust(pval, method = "bon")) %>%
    mutate(qval2 = p.adjust(pval2, method = "bon"))
hist(res$esti, 100)
hist(res$esti2, 100)
up1 = mean(res$esti) + 2*sd(res$esti)
down1 = mean(res$esti) - 2*sd(res$esti)
res_sub1 = res[(res$esti < down1 |
                    res$esti > up1) &
                   res$qval < 0.05,]
up2 = mean(res$esti2) + 2*sd(res$esti2)
down2 = mean(res$esti2) - 2*sd(res$esti2)
res_sub2 = res[(res$esti2 < down2 |
                    res$esti2 > up2) &
                   res$qval2 < 0.05,]
res_sub = res[res$gene_short_name %in% as.vector(res_sub1$gene_short_name) & !res$gene_short_name %in% as.vector(res_sub2$gene_short_name),]
dim(res_sub) # 124
res_sub %>% filter(esti > 0) %>% arrange(qval) %>% head()
res_sub %>% filter(esti < 0) %>% arrange(qval) %>% head()

target_gene_list = c("Gata4","Onecut2","Hnf4g","Cdh6")
df = data.frame(obj_sub[[]])
for (i in 1:length(target_gene_list)){
    gene_i = mouse_gene$gene_ID[mouse_gene$gene_short_name == target_gene_list[i]]
    p = df %>%
        mutate(exp = as.vector(count[gene_i,])) %>%
        ggplot(aes(Pseudotime, exp)) + 
        geom_point(aes(color = celltype)) + 
        geom_smooth() + 
        scale_colour_manual(values = cols) +
        theme_classic(base_size = 15) +
        theme(legend.position="none") +
        geom_vline(xintercept = vline_x, color = "blue") +
        labs(title = target_gene_list[i])
    assign(paste0("p", i), p)
}

ggsave(paste0(work_path, "/midgut_hindgut_pseudobulk_1_add.pdf"), p1 + p2 + p3 + p4, width=10, height=10, dpi = 300, units = "in", useDingbats=FALSE) 

exp = count[as.vector(res_sub$gene_ID), order(obj_sub$Pseudotime)]
rownames(exp) = res_sub$gene_short_name
pdf(paste0(work_path, "/midgut_hindgut_pseudobulk_2_add.pdf"), 10, 10)
g = heatmap.2(as.matrix(exp), col=viridis, scale="row", Rowv = TRUE, Colv = FALSE, key=T, density.info="none", trace="none", cexRow=0.6, cexCol = 0.6)
dev.off()
gene_list = rownames(exp)[rev(g$rowInd)]
write.table(gene_list, paste0(work_path, "/midgut_hindgut_pseudobulk_gene_list_add.txt"), col.names=F, row.names=F, quote=F, sep="\t")

state_list = colnames(exp)[(g$colInd)]
df = df[state_list,]
df_sub = df[df$project == "pijuan",]
vline_xxx = c(1:nrow(df))[df$embryo_id == df_sub$embryo_id[nrow(df_sub)]]
p = ggplot(df, aes(1:nrow(df),Pseudotime)) + geom_point(aes(color = celltype)) + 
    geom_vline(xintercept = vline_xxx, color = "blue") +
    scale_colour_manual(values = cols) +
    theme_classic(base_size = 15) + theme(legend.position="none")
ggsave(paste0(work_path, "/midgut_hindgut_pseudobulk_3_add.pdf"), p, width=8, height=4, dpi = 300, units = "in", useDingbats=FALSE) 






### Cardiomyocytes epithelium 
obj_sub = subset(obj, subset = Anno %in% cardiomyocytes_state)
table(obj_sub$celltype)

count = GetAssayData(object = obj_sub, slot = "counts")
count = t(t(count) / Matrix::colSums(count)) * 100000
count = log(count + 1)

res = calculateDEG(count, obj_sub$Pseudotime, factor(obj_sub$project))
res = res %>% 
    mutate(gene_ID = rownames(res)) %>%
    left_join(mouse_gene[,c("gene_ID", "gene_short_name")], by = "gene_ID") %>%
    mutate(qval = p.adjust(pval, method = "bon")) %>%
    mutate(qval2 = p.adjust(pval2, method = "bon"))
hist(res$esti, 100)
hist(res$esti2, 100)
up1 = mean(res$esti) + 2*sd(res$esti)
down1 = mean(res$esti) - 2*sd(res$esti)
res_sub1 = res[(res$esti < down1 |
                    res$esti > up1) &
                   res$qval < 0.05,]
up2 = mean(res$esti2) + 2*sd(res$esti2)
down2 = mean(res$esti2) - 2*sd(res$esti2)
res_sub2 = res[(res$esti2 < down2 |
                    res$esti2 > up2) &
                   res$qval2 < 0.05,]
res_sub = res[res$gene_short_name %in% as.vector(res_sub1$gene_short_name) & !res$gene_short_name %in% as.vector(res_sub2$gene_short_name),]
dim(res_sub) # 85
res_sub %>% filter(esti > 0) %>% arrange(qval) %>% head()
res_sub %>% filter(esti < 0) %>% arrange(qval) %>% head()

target_gene_list = c("Tnnt2","Myl4","Hbb-bt","Dut")
df = data.frame(obj_sub[[]])
for (i in 1:length(target_gene_list)){
    gene_i = mouse_gene$gene_ID[mouse_gene$gene_short_name == target_gene_list[i]]
    p = df %>%
        mutate(exp = as.vector(count[gene_i,])) %>%
        ggplot(aes(Pseudotime, exp)) + 
        geom_point(aes(color = celltype)) + 
        geom_smooth() + 
        scale_colour_manual(values = cols) +
        theme_classic(base_size = 15) +
        theme(legend.position="none") +
        geom_vline(xintercept = vline_x, color = "blue") +
        labs(title = target_gene_list[i])
    assign(paste0("p", i), p)
}

ggsave(paste0(work_path, "/cardiomyocytes_pseudobulk_1_add.pdf"), p1 + p2 + p3 + p4, width=10, height=10, dpi = 300, units = "in", useDingbats=FALSE) 

exp = count[as.vector(res_sub$gene_ID), order(obj_sub$Pseudotime)]
rownames(exp) = res_sub$gene_short_name
pdf(paste0(work_path, "/cardiomyocytes_pseudobulk_2_add.pdf"), 10, 10)
g = heatmap.2(as.matrix(exp), col=viridis, scale="row", Rowv = TRUE, Colv = FALSE, key=T, density.info="none", trace="none", cexRow=0.6, cexCol = 0.6)
dev.off()
gene_list = rownames(exp)[rev(g$rowInd)]
write.table(gene_list, paste0(work_path, "/cardiomyocytes_pseudobulk_gene_list_add.txt"), col.names=F, row.names=F, quote=F, sep="\t")

state_list = colnames(exp)[(g$colInd)]
df = df[state_list,]
df_sub = df[df$project == "pijuan",]
vline_xxx = c(1:nrow(df))[df$embryo_id == df_sub$embryo_id[nrow(df_sub)]]
p = ggplot(df, aes(1:nrow(df),Pseudotime)) + geom_point(aes(color = celltype)) + 
    geom_vline(xintercept = vline_xxx, color = "blue") +
    scale_colour_manual(values = cols) +
    theme_classic(base_size = 15) + theme(legend.position="none")
ggsave(paste0(work_path, "/cardiomyocytes_pseudobulk_3_add.pdf"), p, width=8, height=4, dpi = 300, units = "in", useDingbats=FALSE) 

