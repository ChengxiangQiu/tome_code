########################
#### Section6_keyTF ####
########################

#### Systematic nomination of key transcription factors for cell type specification

#################################
#### Step1_increase_decrease ####
#################################

### increase key TFs

rm(list=ls())
library(Seurat)
library(dplyr)
source("help_code/help_code.R")
path = ""

args = commandArgs(trailingOnly=TRUE)

time_i = as.numeric(args[1])

time_point = paste0("E", c(3.5, 4.5, 5.25, 5.5, 6.25, seq(6.5, 8.25, 0.25), "8.5a", "8.5b", seq(9.5, 13.5, 1)))

mouse_tf = read.table(paste0(path, "help_code/Mus_musculus_TF.txt"), header=T, as.is=TRUE, sep="\t")
mouse_tf <- as.vector(unique(mouse_tf$Ensembl))

edge = readRDS(paste0(path, "/edge_all.rds"))
edge = edge[edge$prob >= 0.2,]
edge = edge[!edge$nex %in% c("E5.25:Extraembryonic ectoderm","E7:Parietal endoderm"),]

###
getKeyTF <- function(time_i, t1 = 0.7, t2 = 0.3, t3 = 0.3){
  
  
  obj_pre = readRDS(paste0(path, "/obj_", time_point[time_i], "_exp.rds"))
  anno_pre = readRDS(paste0(path, "/obj_", time_point[time_i], "_anno.rds"))
  anno_pre = anno_pre[colnames(obj_pre),]
  obj_pre$Anno = as.vector(anno_pre$Anno)
  Idents(obj_pre) = as.vector(anno_pre$Anno)
  
  exp_pre <- GetAssayData(object = obj_pre, slot = "counts")
  exp_pre <- t(t(exp_pre) / colSums(exp_pre)) * 100000
  exp_pre@x <- log(exp_pre@x + 1)
  
  obj_nex = readRDS(paste0(path, "/revision/exp/obj_", time_point[time_i+1], "_exp.rds"))
  anno_nex = readRDS(paste0(path, "/revision/anno/obj_", time_point[time_i+1], "_anno.rds"))
  anno_nex = anno_nex[colnames(obj_nex),]
  obj_nex$Anno = as.vector(anno_nex$Anno)
  Idents(obj_nex) = as.vector(anno_nex$Anno)
  
  obj_nex <- NormalizeData(obj_nex, normalization.method = "LogNormalize", scale.factor = 10000)
  
  obj_conn <- merge(obj_pre, obj_nex)
  obj_conn <- NormalizeData(obj_conn, normalization.method = "LogNormalize", scale.factor = 10000)
  Idents(obj_conn) <- as.vector(obj_conn$Anno)
  
  result <- NULL
  
  pre_state <- as.vector(unique(edge[edge$pre_time == time_point[time_i],]$pre))
  
  for (i in 1:length(pre_state)){
    print(pre_state[i])
    edge_sub <- edge[edge$pre == pre_state[i],]
    if(nrow(edge_sub)!=1 || edge_sub$pre_cell!=edge_sub$nex_cell){
      
      exp_pre_sub <- exp_pre[,obj_pre$Anno == pre_state[i]]
      exp_pre_sub_mean <- rowMeans(exp_pre_sub)
      DEG_1 = rownames(exp_pre)
      
      gene_use = intersect(mouse_tf, DEG_1)
      
      for(j in 1:nrow(edge_sub)){
        if (edge_sub$pre_cell[j]!=edge_sub$nex_cell[j]){
          print(paste0("running nex state:",j,"/",nrow(edge_sub)))
          ### requirement 2: TF shoule be increasing at child state
          gene_use_sub = intersect(gene_use, rownames(obj_conn))
          DEG_2 <- FindMarkers(obj_conn, ident.1 = edge_sub$nex[j], ident.2 = edge_sub$pre[j], features = gene_use_sub, min.pct = 0, logfc.threshold = t2, only.pos = TRUE)
          DEG_2 <- DEG_2[DEG_2$p_val_adj<0.05 & DEG_2$pct.1 > 0.1,]
          
          if(nrow(DEG_2)!=0){
            
            ### requirement 3: TF should be decreasing at sister state
            if(nrow(edge_sub)!=1){
              gene_use_sub = intersect(gene_use, rownames(obj_nex))
              DEG_3 <- FindMarkers(obj_nex, ident.1 = edge_sub$nex[j], ident.2 = c(as.vector(edge_sub$nex[-j])), features = gene_use_sub, min.pct = 0, logfc.threshold = t3, only.pos = TRUE)
              DEG_3 <- DEG_3[DEG_3$p_val_adj<0.05 & DEG_3$pct.1 > 0.1,]
              
              gene_intersect = intersect(rownames(DEG_2), rownames(DEG_3))
              
              DEG_2 = DEG_2[gene_intersect,]
              DEG_3 = DEG_3[gene_intersect,]
              
              
            } else {
              DEG_3 = DEG_2
              DEG_3[,] = NA
            }
            
            DEG = cbind(DEG_2, DEG_3)
          } else {
            DEG = DEG_2
          }
          
          if(nrow(DEG)>0){
            
            DEG$edge =  paste0(edge_sub$nex[j], "_", edge_sub$pre[j])
            DEG$state = edge_sub$nex[j]
            DEG$gene_short_name = mouse_gene[rownames(DEG),]$gene_short_name
            DEG$gene_ID = rownames(DEG)
            rownames(DEG) = NULL
            
            result = rbind(result, DEG)
          }
        }
      } 
    }
  }
  return(result)
}  

res <- getKeyTF(time_i, t1 = 0, t2 = 0, t3 = 0)
saveRDS(res, paste0(path, '/increase/', time_i, '_tf.rds'))


#### decrease

getKeyTF <- function(time_i, t1 = 0.7, t2 = 0.3, t3 = 0.3){
  
  
  obj_pre = readRDS(paste0(path, "/revision/exp/obj_", time_point[time_i], "_exp.rds"))
  anno_pre = readRDS(paste0(path, "/revision/anno/obj_", time_point[time_i], "_anno.rds"))
  anno_pre = anno_pre[colnames(obj_pre),]
  obj_pre$Anno = as.vector(anno_pre$Anno)
  Idents(obj_pre) = as.vector(anno_pre$Anno)
  
  exp_pre <- GetAssayData(object = obj_pre, slot = "counts")
  exp_pre <- t(t(exp_pre) / colSums(exp_pre)) * 100000
  exp_pre@x <- log(exp_pre@x + 1)
  
  obj_nex = readRDS(paste0(path, "/revision/exp/obj_", time_point[time_i+1], "_exp.rds"))
  anno_nex = readRDS(paste0(path, "/revision/anno/obj_", time_point[time_i+1], "_anno.rds"))
  anno_nex = anno_nex[colnames(obj_nex),]
  obj_nex$Anno = as.vector(anno_nex$Anno)
  Idents(obj_nex) = as.vector(anno_nex$Anno)
  
  obj_nex <- NormalizeData(obj_nex, normalization.method = "LogNormalize", scale.factor = 10000)
  
  obj_conn <- merge(obj_pre, obj_nex)
  obj_conn <- NormalizeData(obj_conn, normalization.method = "LogNormalize", scale.factor = 10000)
  Idents(obj_conn) <- as.vector(obj_conn$Anno)
  
  result <- NULL
  
  pre_state <- as.vector(unique(edge[edge$pre_time == time_point[time_i],]$pre))
  
  for (i in 1:length(pre_state)){
    print(pre_state[i])
    edge_sub <- edge[edge$pre == pre_state[i],]
    if(nrow(edge_sub)!=1 || edge_sub$pre_cell!=edge_sub$nex_cell){
      
      exp_pre_sub <- exp_pre[,obj_pre$Anno == pre_state[i]]
      exp_pre_sub_mean <- rowMeans(exp_pre_sub)
      DEG_1 = rownames(exp_pre)
      
      gene_use = intersect(mouse_tf, DEG_1)
      
      for(j in 1:nrow(edge_sub)){
        if (edge_sub$pre_cell[j]!=edge_sub$nex_cell[j]){
          print(paste0("running nex state:",j,"/",nrow(edge_sub)))
          ### requirement 2: TF shoule be decreasing at child state
          gene_use_sub = intersect(gene_use, rownames(obj_conn))
          DEG_2 <- FindMarkers(obj_conn, ident.1 = edge_sub$pre[j], ident.2 = edge_sub$nex[j], features = gene_use_sub, min.pct = 0, logfc.threshold = t2, only.pos = TRUE)
          DEG_2 <- DEG_2[DEG_2$p_val_adj<0.05 & DEG_2$pct.1 > 0.1,]
          
          if(nrow(DEG_2)!=0){
            
            ### requirement 3: TF should be decreasing at sister state
            if(nrow(edge_sub)!=1){
              gene_use_sub = intersect(gene_use, rownames(obj_nex))
              DEG_3 <- FindMarkers(obj_nex, ident.1 = c(as.vector(edge_sub$nex[-j])), ident.2 = edge_sub$nex[j], features = gene_use_sub, min.pct = 0, logfc.threshold = t3, only.pos = TRUE)
              DEG_3 <- DEG_3[DEG_3$p_val_adj<0.05 & DEG_3$pct.1 > 0.1,]
              
              gene_intersect = intersect(rownames(DEG_2), rownames(DEG_3))
              
              DEG_2 = DEG_2[gene_intersect,]
              DEG_3 = DEG_3[gene_intersect,]
              
              
            } else {
              DEG_3 = DEG_2
              DEG_3[,] = NA
            }
            
            DEG = cbind(DEG_2, DEG_3)
          } else {
            DEG = DEG_2
          }
          
          if(nrow(DEG)>0){
            
            DEG$edge =  paste0(edge_sub$nex[j], "_", edge_sub$pre[j])
            DEG$state = edge_sub$nex[j]
            DEG$gene_short_name = mouse_gene[rownames(DEG),]$gene_short_name
            DEG$gene_ID = rownames(DEG)
            rownames(DEG) = NULL
            
            result = rbind(result, DEG)
          }
        }
      } 
    }
  }
  return(result)
}  

###

res <- getKeyTF(time_i, t1 = 0, t2 = 0, t3 = 0)
saveRDS(res, paste0(path, '/decrease/', time_i, '_tf.rds'))

