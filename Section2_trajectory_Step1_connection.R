#############################
#### Section2_trajectory ####
#############################

#### The scripts used for creating cellular trajectories

##########################
#### Step1_connection ####
##########################

#### Of note, to connect cell states before (and including) E8.5a, 
#### we used the standard workflow - doClusterSeurat function in help_code.R

#### however, because the large number of cells, to connect cell states 
#### after (and including) E8.5b, we used the reciprocal PCA based method shown in below

#### Reference: https://satijalab.org/seurat/archive/v3.0/integration.html
#### under "Standard Workflow", and "Reciprocal PCA" panels

library(Seurat)
library(future)
library(future.apply)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 100000 * 1024^2)
source("help_code/help_code.R")

work_path = ""

time_point = paste0("E", c(3.5, 4.5, 5.25, 5.5, 6.25, seq(6.5, 8.25, 0.25), "8.5a", "8.5b", seq(9.5, 13.5, 1)))

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])
time_1 = time_point[kk]
time_2 = time_point[kk + 1]
print(time_1)
print(time_2)

### adding both time and batch (if necessary) information

obj_1 = readRDS(paste0(work_path, "/obj_", time_1, "_exp.rds"))
obj_1$group = paste0(time_1, "_", obj_1$group)

obj_2 = readRDS(paste0(work_path, "/obj_", time_2, "_exp.rds"))
obj_2$group = paste0(time_2, "_", obj_1$group)

obj = merge(obj_1, obj_2)

if(!time_2 %in% paste0("E", c("8.5b", "9.5", "10.5", "11.5", "12.5", "13.5"))){
    
    obj.integrated = doClusterSeurat(obj)
    
} else {
    
    obj.list <- SplitObject(obj, split.by = "group")
    obj.list <- future_lapply(X = obj.list, FUN = function(x) {
        x <- NormalizeData(x, verbose = FALSE)
        x <- FindVariableFeatures(x, verbose = FALSE)
    })
    
    features <- SelectIntegrationFeatures(object.list = obj.list)
    obj.list <- future_lapply(X = obj.list, FUN = function(x) {
        x <- ScaleData(x, features = features, verbose = FALSE)
        x <- RunPCA(x, features = features, verbose = FALSE)
    })
    
    anchors <- FindIntegrationAnchors(object.list = obj.list, reduction = "rpca", 
                                      dims = 1:50)
    obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
    
    obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
    obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)
    obj.integrated <- RunUMAP(obj.integrated, dims = 1:30, n.components = 3, min.dist = 0.75)
    
}

saveRDS(obj.integrated, paste0(work_path, "/", time_1, "_", time_2, ".rds"))

emb = data.frame(Embeddings(object = obj.integrated, reduction = "umap"))
saveRDS(emb, file=paste0(work_path, "/", time_1, "_", time_2, "_umap3.rds"))


#### running Knn to find ancestor 

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(htmlwidgets)
library(plotly)
library(monocle3)
source("help_code/help_code.R")

work_path = ""

time_point = paste0("E", c(3.5, 4.5, 5.25, 5.5, 6.25, seq(6.5, 8.25, 0.25), "8.5a", "8.5b", seq(9.5, 13.5, 1)))

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])
time_i = time_point[kk]
time_j = time_point[kk+1]
print(time_i)
print(time_j)

emb = readRDS(paste0(work_path, "/", time_i, "_", time_j, "_umap3.rds"))
emb = data.frame(emb)

anno1 = readRDS(paste0(work_path, "/obj_", time_i, "_anno.rds"))
anno1 = anno1[,c("day", "Anno")]
anno1$day = "pre"
anno1$stage = time_i
anno2 = readRDS(paste0(work_path, "/obj_", time_j, "_anno.rds"))
anno2 = anno2[,c("day", "Anno")]
anno2$day = "nex"
anno2$stage = time_j

anno = rbind(anno1, anno2)
if(nrow(emb) != nrow(anno)){
    print("Error!")
    print(xxx)
}
anno = anno[rownames(emb),]

res = createLineage_Knn(emb, anno,  k_neigh = 5) #### createLineage_Knn function was in help_code.R
saveRDS(res, paste0(work_path, "/", time_i, "_", time_j, "_Knn_umap.rds"))

#### testing if results are robust by different k parameter
for(k_i in c(8, 10, 15, 20)){
    res = createLineage_Knn(emb, anno,  k_neigh = k_i) #### createLineage_Knn function was in help_code.R
    saveRDS(res, paste0(work_path, "/", time_i, "_", time_j, "_Knn_umap_k_", k_i, ".rds"))
}

#### creating the median value of the matrix
replication_times=500
dat = res
state_1 = row.names(dat[[1]])
state_2 = names(dat[[1]])
tmp_1 = matrix(NA,nrow(dat[[1]]),ncol(dat[[1]]))
for(i in 1:nrow(dat[[1]])){
    for(j in 1:ncol(dat[[1]])){
        xx = NULL
        for(k in 1:replication_times){
            xx = c(xx, dat[[k]][i,j])
        }
        tmp_1[i,j] = median(xx[!is.na(xx)])
    }
}
tmp_1 = data.frame(tmp_1)
row.names(tmp_1) = state_1
names(tmp_1) = state_2

write.csv(tmp_1, paste0(work_path, "/", time_i, "_", time_j, "_Knn_umap.csv"))


#### repeat this approach but using 30 PCs instead to testing if results are robust

obj = readRDS(paste0(work_path, "/", time_i, "_", time_j, ".rds"))
emb = data.frame(Embeddings(obj, reduction = "pca"))
print(dim(emb))

anno1 = readRDS(paste0(work_path, "/obj_", time_i, "_anno.rds"))
anno1 = anno1[,c("day", "Anno")]
anno1$day = "pre"
anno1$stage = time_i
anno2 = readRDS(paste0(work_path, "/obj_", time_j, "_anno.rds"))
anno2 = anno2[,c("day", "Anno")]
anno2$day = "nex"
anno2$stage = time_j

anno = rbind(anno1, anno2)
if(nrow(emb) != nrow(anno)){
    print("Error!")
    print(xxx)
}
anno = anno[rownames(emb),]

res = createLineage_Knn(emb, anno, reduction = "pca", replication_times = 100) #### createLineage_Knn function was in help_code.R
saveRDS(res, paste0(work_path, "/", time_i, "_", time_j, "_Knn_pca.rds"))


