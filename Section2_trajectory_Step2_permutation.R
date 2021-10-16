#############################
#### Section2_trajectory ####
#############################

#### The scripts used for creating cellular trajectories

###########################
#### Step2_permutation ####
###########################

#### performing a permutation (shuffling the Anno labels for each connection) by 1,000 times

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(htmlwidgets)
library(plotly)
library(monocle3)
library(FNN)
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
print(dim(emb))

if(ncol(emb) != 3){
    print(XXX)
}

anno1 = readRDS(paste0(work_path, "/obj_", time_i, "_anno.rds"))
anno1 = anno1[,c("day", "Anno")]
anno1$day = "pre"
anno1$stage = time_i
anno2 = readRDS(paste0(work_path, "/obj_", time_j, "_anno.rds"))
anno2 = anno2[,c("day", "Anno")]
anno2$day = "nex"
anno2$stage = time_j

permutation_times = 1000
k_neigh = 5

res = list()

for(rep_i in 1:permutation_times){
    
    anno1$state = anno1$Anno[sample(1:nrow(anno1))]
    anno2$state = anno2$Anno[sample(1:nrow(anno2))]
    
    anno = rbind(anno1, anno2)
    if(nrow(emb) != nrow(anno)){
        print("Error!")
        print(xxx)
    }
    pd = anno[rownames(emb),]
    
    emb_sub = emb
    pd_sub = pd
    
    irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",]
    irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",]
    pd_sub1 <- pd_sub[pd_sub$day == "pre",]
    pd_sub2 <- pd_sub[pd_sub$day == "nex",]
    
    pre_state_min = min(table(as.vector(pd_sub1$state)))
    
    if (pre_state_min < k_neigh & pre_state_min >= 3){
        k_neigh = pre_state_min
        print(k_neigh)
    }
    
    if (pre_state_min < 3){
        next
    }
    
    neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index
    
    tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
    for(i in 1:k_neigh){
        tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
    }
    state1 <- names(table(as.vector(pd_sub1$state)))
    state2 <- names(table(as.vector(pd_sub2$state)))
    
    tmp2 <- matrix(NA,length(state2),length(state1))
    for(i in 1:length(state2)){
        x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
        for(j in 1:length(state1)){
            tmp2[i,j] <- sum(x==state1[j])
        }
    }
    tmp2 <- tmp2/apply(tmp2,1,sum)
    tmp2 <- data.frame(tmp2)
    row.names(tmp2) = state2
    names(tmp2) = state1
    
    res[[rep_i]] = tmp2
    
}

saveRDS(res, paste0(work_path, "/", time_i, "_", time_j, "_Knn_umap_permutation.rds"))

x = list()
z = NULL
for(cnt in 1:(length(time_point)-1)){
    time_i = time_point[cnt]
    time_j = time_point[cnt+1]
    
    dat = readRDS(paste0(work_path, "/", time_i, "_", time_j, "_Knn_umap_permutation.rds"))
    
    permutation_times = 1000
    y = NULL
    
    for(i in 1:permutation_times){
        y = c(y, as.vector(as.matrix(dat[[1]])))
    }
    
    x[[cnt]] = y
    z = c(z, y)
    
    print(paste0(time_i, ":", sum(y > 0.2)/length(y)))
}

print(sum(z >= 0.2)/length(z))
### 0.0089

library(ggplot2)
dat = data.frame(edge_weights_by_permutation = z)

p<-ggplot(dat, aes(x=edge_weights_by_permutation)) +
    geom_histogram(position="identity", alpha=0.5, binwidth=0.01) + 
    geom_vline(xintercept = 0.2, colour = "red") +
    theme_classic(base_size = 15)

