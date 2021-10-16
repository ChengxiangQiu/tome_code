#############################
#### Section5_pseudobulk ####
#############################

#### Inferring the molecular histories of individual cell types

##################################
#### Step1_aggregating_counts ####
##################################


#### E9.5 to E13.5, the deeper sequencing of the previous libraries (Cao et al.)

library(Seurat)
library(dplyr)
orig_data_path = ""
work_path = ""
time_point = paste0("E", seq(9.5, 13.5, 1))

cnt = 1; time_i = time_point[cnt]
obj = readRDS(paste0(work_path, "/obj_", time_i, "_exp.rds"))
gene = rownames(obj)

pd_all = NULL
count = NULL

for(cnt in 1:5){
    time_i = time_point[cnt]
    print(time_i)
    
    pd = readRDS(paste0(orig_data_path, "/", time_i, "_pd.rds"))
    anno = readRDS(paste0(work_path, "/revision/anno/obj_", time_i, "_anno.rds"))
    tmp = anno %>%
        select(sample, celltype) %>%
        left_join(pd %>% select(sample, embryo_id, embryo_sex), by = "sample")
    tmp = data.frame(tmp)
    rownames(tmp) = as.vector(tmp$sample)
    tmp$embryo_id = gsub("embryo", "cao", as.vector(tmp$embryo_id))
    tmp$project = "cao"
    tmp = tmp[,c("embryo_id", "embryo_sex", "project", "celltype")]
    tmp$day = time_i
    pd_all = rbind(pd_all, tmp)
    
    obj = readRDS(paste0(work_path, "/obj_", time_i, "_exp.rds"))
    xx = GetAssayData(object = obj, slot = "counts")
    xx = xx[,rownames(tmp)]
    if(sum(rownames(obj) != gene) == 0){
        print("YES!!")
        count = cbind(count, xx)
    } else {
        print("Error!!")
    }
    
}
pd_all$group = "95"
saveRDS(count, paste0(work_path, "/count_4.rds"))
saveRDS(pd_all, paste0(work_path, "/pd_4.rds"))

#### E6.5 to E8.5, from the Pijuan-Sala's data

library(Seurat)
work_path = ""
time_point = paste0("E", c(seq(6.5, 8.25, 0.25), "8.5a"))
gene_list = rownames(readRDS(paste0(work_path, "/obj_", time_point[2], "_exp.rds")))

exp <- NULL
pd <- NULL

for(i in 1:length(time_point)){
    time_i = time_point[i]
    print(time_i)
    
    obj = readRDS(paste0(work_path, "/obj_", time_i, "_exp.rds"))
    anno = readRDS(paste0(work_path, "/obj_", time_i, "_anno.rds"))
    anno = anno[colnames(obj),]
    anno$celltype = unlist(lapply(as.vector(anno$Anno), function(x) strsplit(x,"[:]")[[1]][2]))
    obj$embryo = as.vector(paste0("pijuan_", obj$embryo))
    obj$celltype = as.vector(anno$celltype)
    obj$group = as.vector(anno$group)
    
    if(i==1){
        obj$tmp = !obj$group %in% c("351", "352")
        obj = subset(obj, subset = tmp)
    }

    all <- GetAssayData(object = obj, slot = "counts")
    print(sum(!gene_list %in% rownames(obj)))
    all <- all[gene_list,]
    
    exp = cbind(exp, all)
    tmp = data.frame(obj[[]])
    tmp$project = "pijuan"
    tmp$day = time_i
    pd = rbind(pd, tmp[,c("embryo","project","celltype","day","group")])
}
pd$embryo_sex = NA
pd = pd[,c("embryo", "embryo_sex", "project", "celltype", "day","group")]
colnames(pd) = c("embryo_id", "embryo_sex", "project", "celltype", "day", "group")
saveRDS(exp, paste0(work_path, "/count_2.rds"))
saveRDS(pd, paste0(work_path, "/pd_2.rds"))


#### new E8.5b embryo data

library(dplyr)
library(Matrix)
work_path = ""
dat = readRDS(paste0(orig_data_path, "/dat.rds"))
gene_count = dat$gene_count

anno = readRDS(pasteo(work_path, "/obj_E8.5b_anno.rds"))
gene_count = gene_count[,rownames(anno)]
anno = anno %>%
    left_join(embryo_sex %>% select(RT_group, embryo_id, embryo_sex), by = "RT_group")
rownames(anno) = as.vector(anno$sample)
anno$project = "beth"
anno = anno[,c("embryo_id", "embryo_sex", "project", "celltype", "day")]
anno$group = "852"
saveRDS(gene_count, paste0(work_path, "/count_3.rds"))
saveRDS(anno, paste0(work_path, "/pd_3.rds"))


#### E6.25:Epiblast, including 4 samples

obj = readRDS(paste0(work_path, "/obj_E6.25_exp.rds"))
anno = readRDS(paste0(work_path, "/obj_E6.25_anno.rds"))
obj$celltype = unlist(lapply(as.vector(anno$Anno), function(x) strsplit(x,"[:]")[[1]][2]))
obj$embryo_id = paste0("cheng_", obj$embryo)
obj$embryo_sex = NA
obj$project = "cheng"
obj$day = "E6.25"
obj = subset(obj, subset = celltype == "Epiblast")
obj = subset(obj, subset = embryo_id != "cheng_9") ### cheng_9 only has 2 cells, so we exclude it
exp = GetAssayData(object = obj, slot = "counts")
pd = data.frame(obj[[]])[,c("embryo_id", "embryo_sex", "project", "celltype", "day")]
pd$group = "625"
saveRDS(exp, paste0(work_path, "/count_1.rds"))
saveRDS(pd, paste0(work_path, "/pd_1.rds"))



#### step 2 - aggregating cells for each individual embryo, excluding ExE tissues


library(Seurat)
library(dplyr)
library(monocle3)
work_path = ""
exclude_celltype = c("Embryonic visceral endoderm", "Extraembryonic visceral endoderm", "Parietal endoderm", "Extraembryonic ectoderm")

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])
#### kk = 1:4

pd = readRDS(paste0(work_path, "/pd_", kk, ".rds"))
pd_sub = pd[!pd$celltype %in% exclude_celltype,]
count = readRDS(paste0(work_path, "/count_", kk, ".rds"))

if(sum(rownames(pd) != colnames(count)) != 0){
    print(XXX)
} else {
    embryo_list = as.vector(unique(pd_sub$embryo_id))
    exp = NULL
    for(i in 1:length(embryo_list)){
        print(paste0(i,"/",length(embryo_list)))
        exp_tmp = count[,colnames(count) %in% as.vector(rownames(pd_sub)[pd_sub$embryo_id == embryo_list[i]])]
        exp = cbind(exp, Matrix::rowSums(exp_tmp))
    }
    colnames(exp) = embryo_list
    rownames(exp) = rownames(count)
}

saveRDS(exp, paste0(work_path, "/embryo_", kk, ".rds"))




