########################
#### Section6_keyTF ####
########################

#### Systematic nomination of key transcription factors for cell type specification

#################################
#### Step1_summarize_results ####
#################################

source("help_code/help_code.R")
path = ""

time_point = paste0("E", c(3.5, 4.5, 5.25, 5.5, 6.25, seq(6.5, 8.25, 0.25), "8.5a", "8.5b", seq(9.5, 13.5, 1)))

edge = readRDS(paste0(path, "/edge_all.rds"))
rownames(edge) = paste0(edge$nex, "_", edge$pre)

mouse_tf = read.table(paste0(path, "/Mus_musculus_TF.txt"), header=T, as.is=TRUE, sep="\t")
mouse_tf <- as.vector(unique(mouse_tf$Ensembl))

### increase 
result <- NULL

for (time_i in 1:(length(time_point)-1)){
  print(paste0(time_i, "/", length(time_point)))
  res <- readRDS(paste0(path, '/increase/', time_i, '_tf.rds'))
  if (!is.null(res)){
    if(nrow(res)>0){
      res = data.frame(res)
      names(res) = c("p_val_1", "avg_logFC_1", "pct_1_1", "pct_2_1", "p_val_adj_1", "p_val_2", "avg_logFC_2", "pct_1_2", "pct_2_2", "p_val_adj_2", "edge", "state", "gene_short_name", "gene_ID")
      res = res[,c("gene_ID", "avg_logFC_1", "pct_1_1", "pct_2_1", "p_val_adj_1", "avg_logFC_2", "pct_1_2", "pct_2_2", "p_val_adj_2", "edge", "state", "gene_short_name")]
      
      edge_list = as.vector(unique(res$edge))
      for(i in 1:length(edge_list)){
        tmp = res[res$edge == edge_list[i],]
        if(sum(is.na(tmp$pct_1_2))==0){
          x = scale(tmp$avg_logFC_1)
          y = scale(tmp$avg_logFC_2)
          x = (x + y)/2
        } else {
          x = scale(tmp$avg_logFC_1)
        }
        tmp$score = x
        result = rbind(result, tmp)
      }
    }
  }
}

edge_tmp = edge[as.vector(result$edge),]
result$prob = edge_tmp$prob

write.csv(result,paste0(path, "/mouse_keyTF_increase_all.csv"), row.names=F)

#### further only retain results from edges which are corresponding to newly emerged cell types

### decrease
result <- NULL

for (time_i in 1:(length(time_point)-1)){
  print(paste0(time_i, "/", length(time_point)))
  res <- readRDS(paste0(path, '/decrease/', time_i, '_tf.rds'))
  if (!is.null(res)){
    if(nrow(res)>0){
      res = data.frame(res)
      names(res) = c("p_val_1", "avg_logFC_1", "pct_1_1", "pct_2_1", "p_val_adj_1", "p_val_2", "avg_logFC_2", "pct_1_2", "pct_2_2", "p_val_adj_2", "edge", "state", "gene_short_name", "gene_ID")
      res = res[,c("gene_ID", "avg_logFC_1", "pct_1_1", "pct_2_1", "p_val_adj_1", "avg_logFC_2", "pct_1_2", "pct_2_2", "p_val_adj_2", "edge", "state", "gene_short_name")]
      
      edge_list = as.vector(unique(res$edge))
      for(i in 1:length(edge_list)){
        tmp = res[res$edge == edge_list[i],]
        if(sum(is.na(tmp$pct_1_2))==0){
          x = scale(tmp$avg_logFC_1)
          y = scale(tmp$avg_logFC_2)
          x = (x + y)/2
        } else {
          x = scale(tmp$avg_logFC_1)
        }
        tmp$score = x
        result = rbind(result, tmp)
      }
    }
  }
}

edge_tmp = edge[as.vector(result$edge),]
result$prob = edge_tmp$prob

write.csv(result,paste0(path, "/mouse_keyTF_decrease_all.csv"), row.names=F)

#### further only retain results from edges which are corresponding to newly emerged cell types

  
