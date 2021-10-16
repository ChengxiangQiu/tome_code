#############################
#### Section2_trajectory ####
#############################

#### The scripts used for creating cellular trajectories

#################################
#### Step3_summarize_results ####
#################################

library(reshape2)
library(scales)
library(ggplot2)
library(dplyr)
library(gplots)
library(viridis)

work_path = ""

time_point = paste0("E", c(3.5, 4.5, 5.25, 5.5, 6.25, seq(6.5, 8.25, 0.25), "8.5a", "8.5b", seq(9.5, 13.5, 1)))

replication_times=500
res_median_umap = list()
for(time_i in 1:(length(time_point)-1)){
  print(paste0(time_point[time_i], ":", time_point[time_i+1]))
  dat = readRDS(paste0(work_path, "/",time_point[time_i],"_",time_point[time_i+1],"_Knn_umap.rds"))
  state_1 = row.names(dat[[1]])
  state_1 = paste0(time_point[time_i+1], ":", update_name(gsub(paste0(time_point[time_i+1], ":"), "", state_1)))
  state_2 = names(dat[[1]])
  state_2 = paste0(time_point[time_i], ":", update_name(gsub(paste0(time_point[time_i], ":"), "", state_2)))
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
  res_median_umap[[time_i]] = tmp_1
}

dat = NULL
for(i in 1:length(res_median_umap)){
  print(time_point[i])
  dat = rbind(dat, melt(as.matrix(res_median_umap[[i]])))
}

dat = data.frame(dat)
names(dat) = c("nex", "pre", "prob")

#### As described in the manuscript, anchor-based batch correction and integration 
#### of profiles of E8.5 cells from (Pijuan-Sala et al. 2019) (termed “E8.5a”) 
#### and newly generated profiles of E8.5 nuclei (termed “E8.5b”) worked very 
#### well with the exception of primitive erythroid cells, which we suspect may 
#### be due to more extensive differences between cells vs. nuclei in this cell type.
#### we manually set E8.5b:Primitive erythroid cells as the ancestor of E8.5b:Primitive erythroid cells

dat$prob[dat$nex == "E8.5b:Primitive erythroid cells"] = 0.0
dat$prob[dat$nex == "E8.5b:Primitive erythroid cells" & dat$pre == "E8.5a:Primitive erythroid cells"] = 1.0

dat$pre_time = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][1]))
dat$pre_cell = unlist(lapply(as.vector(dat$pre), function(x) strsplit(x,"[:]")[[1]][2]))
dat$nex_time = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][1]))
dat$nex_cell = unlist(lapply(as.vector(dat$nex), function(x) strsplit(x,"[:]")[[1]][2]))

saveRDS(dat, paste0(work_path, "/edge_all.rds"))

### Here we use “cell state” to mean an annotated cluster at a given stage. 

print(paste0("how many edges: ", nrow(dat)))
print(paste0("how many edges (> 0): ", nrow(dat[dat$prob>0,])))
print(paste0("how many edges (> 0.2): ", nrow(dat[dat$prob>=0.2,])))
print(paste0("how many edges (> 0.7): ", nrow(dat[dat$prob>=0.7,])))
print(paste0("how many edges (> 0.8): ", nrow(dat[dat$prob>=0.8,])))
print(paste0("how many nodes: ", length(unique(c(as.vector(dat$pre), as.vector(dat$nex))))))
print(paste0("how many cell types: ", length(unique(c(as.vector(dat$pre_cell), as.vector(dat$nex_cell))))))

#### extract edges with prob > 0.2

x = dat[dat$prob>=0.2,]
x = x[,c("pre","nex","prob")]
print(paste0("how many nodes now: ", length(unique(c(as.vector(x$pre), as.vector(x$nex))))))

#### As described in the manuscript, we introduced 4 “dummy nodes”, corresponding 
#### to morula at E3.0 (as a root for trophectoderm and inner cell mass),
#### trophectoderm at E3.5 and E4.5 (which had been removed at these timepoints 
#### by immunosurgery (Mohammed et al. 2017)) and parietal endoderm at E6.75 
#### (undetected, likely due to undersampling). For technical reasons (see above), 
#### we also introduced an edge between primitive erythroid cells at E8.5a and E8.5b.

x = x[!x$nex %in% c("E5.25:Extraembryonic ectoderm", "E7:Parietal endoderm"),]
dummy = NULL
dummy = rbind(dummy, c("E3:Morula","E3.5:Inner cell mass",1))
dummy = rbind(dummy, c("E3:Morula","E3.5:Trophectoderm",1))
dummy = rbind(dummy, c("E3.5:Trophectoderm","E4.5:Trophectoderm",1))
dummy = rbind(dummy, c("E4.5:Trophectoderm","E5.25:Extraembryonic ectoderm",1))
dummy = rbind(dummy, c("E6.5:Parietal endoderm","E6.75:Parietal endoderm",1))
dummy = rbind(dummy, c("E6.75:Parietal endoderm","E7:Parietal endoderm",1))
dummy = data.frame(dummy)
names(dummy) = c("pre","nex","prob")

res = rbind(x, dummy)

dat_sub = res
dat_sub$pre_cell = unlist(lapply(as.vector(dat_sub$pre), function(x) strsplit(x,"[:]")[[1]][2]))
dat_sub$nex_cell = unlist(lapply(as.vector(dat_sub$nex), function(x) strsplit(x,"[:]")[[1]][2]))
sum(dat_sub$pre_cell == dat_sub$nex_cell)
sum(dat_sub$pre_cell == dat_sub$nex_cell)/nrow(dat_sub)
sum(dat_sub$pre_cell != dat_sub$nex_cell)
sum(dat_sub$pre_cell != dat_sub$nex_cell)/nrow(dat_sub)

#### summary on finally used to create the tree

print(paste0("how many edges: ", nrow(res)))
print(paste0("how many nodes: ", length(unique(c(as.vector(res$pre), as.vector(res$nex))))))

nex_list = as.vector(unique(res$nex))
tree = NULL
for(i in 1:length(nex_list)){
  res_sub = res[res$nex==nex_list[i],]
  if(nrow(res_sub)==1){
    tree = rbind(tree, res_sub)
  } else {
    res_sub = res_sub[order(res_sub$prob, decreasing = TRUE),]
    tree = rbind(tree, res_sub[1,])
  }
}
tree = data.frame(tree)
tree = tree[,c("pre","nex")]

write.table(res, paste0(work_path, "/edge_prob.txt"), row.names = F, col.names = F, quote = F, sep = "\t")
write.table(tree, paste0(work_path, "/mouse_edge.txt"), row.names = F, col.names = F, quote = F, sep = "\t")

#### create heatmap of each pair of adjacent stages

dat = readRDS(paste0(work_path, "/edge_all.rds"))
time_point = paste0("E", c(3.5, 4.5, 5.25, 5.5, 6.25, seq(6.5, 8.25, 0.25), "8.5a", "8.5b", seq(9.5, 13.5, 1)))
mm_cell_type_list = read.table(paste0("help_code/celltype_group.txt"), as.is=T, sep="\t")
mm_cell_type_list = as.vector(mm_cell_type_list$V1)
save_path = ""

kk = 1
df = data.frame(matrix(c(1,1,0,0),2,2))
rownames(df) = c("Epiblast","Hypoblast")
colnames(df) = c("Inner cell mass","Inner cell mass_NA")
pdf(paste0(save_path, "heatmap_",kk,".pdf"),10,10)
print(heatmap.2(as.matrix(df), col=viridis, scale="none", Rowv = FALSE, Colv = FALSE, key=T, density.info="none", trace="none", cexRow=0.6, cexCol = 0.6))
dev.off()
write.table(rownames(df),paste0(save_path, "rownames_",kk,".txt"), row.names=F, col.names=F, sep="\t", quote=F)
write.table(colnames(df),paste0(save_path, "colnames_",kk,".txt"), row.names=F, col.names=F, sep="\t", quote=F)

for(kk in 2:19){
  print(paste0(kk,"/",19))
  dat_sub = dat[dat$pre_time == time_point[kk],c("pre_cell","nex_cell","prob")]
  df = dcast(dat_sub, nex_cell~pre_cell)
  rownames(df) <- df[,1]; df <- df[,-1]
  x = mm_cell_type_list[mm_cell_type_list %in% rownames(df)]
  y = mm_cell_type_list[mm_cell_type_list %in% colnames(df)]
  df = df[x,y]
  pdf(paste0(save_path, "heatmap_",kk,".pdf"),10,10)
  print(heatmap.2(as.matrix(df), col=viridis, scale="none", Rowv = FALSE, Colv = FALSE, key=T, density.info="none", trace="none", cexRow=0.6, cexCol = 0.6))
  dev.off()
  write.table(rownames(df),paste0(save_path, "rownames_",kk,".txt"), row.names=F, col.names=F, sep="\t", quote=F)
  write.table(colnames(df),paste0(save_path, "colnames_",kk,".txt"), row.names=F, col.names=F, sep="\t", quote=F)
}


