################################
#### Section8_align_species ####
################################

#### Systematic comparison of the cellular trajectories of mouse, zebrafish and frog embryogenesis

#################################
#### Step1_aggregating_state ####
#################################


rm(list=ls())
path = ""
library(Seurat)
source("help_code/help_code.R")
library(preprocessCore)

time_point = paste0("E", c(3.5,4.5,5.25,5.5,6.25,seq(6.5,8.25,0.25),"8.5a","8.5b",seq(9.5,13.5)))

exp = list()
cell_num = NULL

### Of note, here we used the gene list by intersecting the Pijuan-Sala's and Beth/Jun's data

obj_1 = readRDS(paste0(path, "/obj_", "E8.5a", "_exp.rds"))
obj_2 = readRDS(paste0(path, "/obj_", "E8.5b", "_exp.rds"))
gene = intersect(rownames(obj_1), rownames(obj_2))

for(time_i in 1:length(time_point)){
  
  print(time_i)
  
  obj = readRDS(paste0(path, "/obj_", time_point[time_i], "_exp.rds"))
  anno = readRDS(paste0(path, "/obj_", time_point[time_i], "_anno.rds"))
  
  if(time_point[time_i] == "E8.5b"){
    anno$project = "dat_354"
  } 
  
  if(time_point[time_i] %in% paste0("E", seq(9.5,13.5))){
    anno$project = "dat_355"
  }
  
  anno = anno[colnames(obj),]
  obj$Anno = anno$Anno
  obj$project = anno$project
  
  if(time_point[time_i]=="E6.5"){
    obj = subset(obj, subset = project == "dat_353")
  }
  
  count = GetAssayData(object = obj, slot = "counts")
  count = t(t(count) / colSums(count)) * 100000
  count@x = log(count@x + 1)
  
  count = count[rownames(count) %in% gene,]
  
  state_list = as.vector(unique(obj$Anno))
  
  for(state_i in 1:length(state_list)){
    print(paste0(time_point[time_i], ":", state_i, "/", length(state_list)))
    
    count_sub = count[,obj$Anno == state_list[state_i]]
    
    if(ncol(count_sub)>1){
      exp[[state_list[state_i]]] = Matrix::rowSums(count_sub)
    } else {
      exp[[state_list[state_i]]] = count_sub
    }
    
    cell_num = rbind(cell_num, data.frame(Anno = state_list[state_i], num = ncol(count_sub)))
  }
  
}


state_list = names(exp)
celltype_list = unlist(lapply(state_list, function(x) strsplit(x,"[:]")[[1]][2]))

dat = NULL
for(i in 1:length(state_list)){
  t = exp[[state_list[i]]][gene]
  print(paste0(state_list[[i]], ":", sum(is.na(t))))
  t[is.na(t)] = 0
  dat = cbind(dat, t)
}

dat = data.frame(dat)
rownames(dat) = gene
colnames(dat) = state_list

saveRDS(dat, paste0(path, "/mouse_state_sum.rds"))
saveRDS(cell_num, paste0(path, "/mouse_cell_num.rds"))



### aggregating cells for every cell state of Zebrafish

time_point <- paste0("hpf", c(3.3, 3.8, 4.3, 4.8, 5.3, 6, 7, 8, 9, 10, 11, 12, 14, 18, 24))

obj_all = readRDS(paste0(path, "/obj.two.rds"))
count_all = GetAssayData(object = obj_all, slot = "counts")

exp = list()
cell_num = NULL

for(time_i in 1:length(time_point)){
  
  print(time_i)
  
  anno = readRDS(paste0(path, "/anno/obj_", time_point[time_i], "_anno.rds"))
  count_sub = count_all[,rownames(anno)]
  obj = CreateSeuratObject(counts = count_sub, min.cells = 0, min.features = 0, meta.data = anno)

  count = GetAssayData(object = obj, slot = "counts")
  count = t(t(count) / colSums(count)) * 100000
  count@x = log(count@x + 1)
  
  state_list = as.vector(unique(obj$Anno))
  
  for(state_i in 1:length(state_list)){
    print(paste0(time_point[time_i], ":", state_i, "/", length(state_list)))
    
    count_sub = count[,obj$Anno == state_list[state_i]]
    
    if(sum(obj$Anno == state_list[state_i])>1){
      exp[[state_list[state_i]]] = Matrix::rowSums(count_sub)
    } else {
      exp[[state_list[state_i]]] = count_sub
    }
    
    cell_num = rbind(cell_num, data.frame(Anno = state_list[state_i], num = sum(obj$Anno == state_list[state_i])))
  }
  
}

gene = NULL
for(i in 1:length(exp)){
  gene = union(gene, names(exp[[i]]))
}

state_list = names(exp)
celltype_list = unlist(lapply(state_list, function(x) strsplit(x,"[:]")[[1]][2]))
state_list = state_list[celltype_list != "Unknown"]

dat = NULL
for(i in 1:length(state_list)){
  t = exp[[state_list[i]]][gene]
  print(paste0(state_list[[i]], ":", sum(is.na(t))))
  t[is.na(t)] = 0
  dat = cbind(dat, t)
}

dat = data.frame(dat)
rownames(dat) = gene
colnames(dat) = state_list

saveRDS(dat, paste0(path, "/zf_state_sum.rds"))
saveRDS(cell_num, paste0(path, "/zf_cell_num.rds"))



### aggregating cells for every cell state of Xenopus

time_point <- paste0("S", c(8,10,11,12,13,14,16,18,20,22))

obj_all <- readRDS("/obj_all.rds")
count_all = GetAssayData(object = obj_all, slot = "counts")

exp = list()
cell_num = NULL

for(time_i in 1:length(time_point)){
  
  print(time_i)
  
  anno = readRDS(paste0(path, "/anno/obj_", time_point[time_i], "_anno.rds"))
  count_sub = count_all[,rownames(anno)]
  obj = CreateSeuratObject(counts = count_sub, min.cells = 0, min.features = 0, meta.data = anno)
  
  count = GetAssayData(object = obj, slot = "counts")
  count = t(t(count) / colSums(count)) * 100000
  count@x = log(count@x + 1)
  
  state_list = as.vector(unique(obj$Anno))
  
  for(state_i in 1:length(state_list)){
    print(paste0(time_point[time_i], ":", state_i, "/", length(state_list)))
    
    count_sub = count[,obj$Anno == state_list[state_i]]
    
    if(ncol(count_sub)>1){
      exp[[state_list[state_i]]] = Matrix::rowSums(count_sub)
    } else {
      exp[[state_list[state_i]]] = count_sub
    }
    
    cell_num = rbind(cell_num, data.frame(Anno = state_list[state_i], num = ncol(count_sub)))
  }
  
}

gene = NULL
for(i in 1:length(exp)){
  gene = union(gene, names(exp[[i]]))
}

state_list = names(exp)
celltype_list = unlist(lapply(state_list, function(x) strsplit(x,"[:]")[[1]][2]))
state_list = state_list[celltype_list != "Unknown"]

dat = NULL
for(i in 1:length(state_list)){
  t = exp[[state_list[i]]][gene]
  print(paste0(state_list[[i]], ":", sum(is.na(t))))
  t[is.na(t)] = 0
  dat = cbind(dat, t)
}

dat = data.frame(dat)
rownames(dat) = gene
colnames(dat) = state_list

saveRDS(dat, paste0(path, "/frog_state_sum.rds"))
saveRDS(cell_num, paste0(path, "/frog_cell_num.rds"))


###################################################################
### generate mean expression for cell states/types for mm/zf/xp ###
###################################################################

work_path = ""

mm_state_sum = readRDS(paste0(work_path, "/","mouse_state_sum.rds"))
mm_cell_num = readRDS(paste0(work_path, "/","mouse_cell_num.rds"))

zf_state_sum = readRDS(paste0(work_path, "/","zf_state_sum.rds"))
zf_cell_num = readRDS(paste0(work_path, "/","zf_cell_num.rds"))

xp_state_sum = readRDS(paste0(work_path, "/","frog_state_sum.rds"))
xp_cell_num = readRDS(paste0(work_path, "/","frog_cell_num.rds"))

generate_cell_state_type <- function(state_sum, cell_num){
  
  cell_num$celltype = unlist(lapply(as.vector(cell_num$Anno), function(x) strsplit(x,"[:]")[[1]][2]))
  rownames(cell_num) = as.vector(cell_num$Anno)
  cell_num = cell_num[colnames(state_sum),]
  
  cell_state_matrix = t(t(as.matrix(state_sum))/as.vector(cell_num$num))
  
  celltype_list = as.vector(unique(cell_num$celltype))
  xx = NULL
  for(i in 1:length(celltype_list)){
    if(sum(cell_num$celltype == celltype_list[i]) > 1){
      yy = rowSums(state_sum[,cell_num$celltype == celltype_list[i]])
    } else {
      yy = state_sum[,cell_num$celltype == celltype_list[i]]
    }
    yy = yy/sum(cell_num$num[cell_num$celltype == celltype_list[i]])
    xx = cbind(xx, yy)
  }
  colnames(xx) = celltype_list
  cell_type_matrix = xx
  
  return(list(cell_state_matrix, cell_type_matrix, cell_num))
}

res = generate_cell_state_type(mm_state_sum, mm_cell_num)
saveRDS(res[[1]], paste0(work_path, "/", "mm_cell_state_matrix.rds"))
saveRDS(res[[2]], paste0(work_path, "/", "mm_cell_type_matrix.rds"))
saveRDS(res[[3]], paste0(work_path, "/", "mm_cell_num.rds"))

res = generate_cell_state_type(zf_state_sum, zf_cell_num)
saveRDS(res[[1]], paste0(work_path, "/", "zf_cell_state_matrix.rds"))
saveRDS(res[[2]], paste0(work_path, "/", "zf_cell_type_matrix.rds"))
saveRDS(res[[3]], paste0(work_path, "/", "zf_cell_num.rds"))

res = generate_cell_state_type(xp_state_sum, xp_cell_num)
saveRDS(res[[1]], paste0(work_path, "/", "xp_cell_state_matrix.rds"))
saveRDS(res[[2]], paste0(work_path, "/", "xp_cell_type_matrix.rds"))
saveRDS(res[[3]], paste0(work_path, "/", "xp_cell_num.rds"))



