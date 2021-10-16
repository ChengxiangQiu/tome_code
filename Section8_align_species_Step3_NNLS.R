################################
#### Section8_align_species ####
################################

#### Systematic comparison of the cellular trajectories of mouse, zebrafish and frog embryogenesis

####################
#### Step3_NNLS ####
####################

rm(list=ls())
work_path = ""
source("help_code/help_code.R")
source("help_code/Cao_NNLS_regression.R")
library(dplyr)
library(tidyr)
library("gplots")
library("reshape2")
library(viridis)

mm_zf = read.table(paste0(work_path, "/help_code/", "mm_zf.txt"), header=T, row.names=1, as.is=T)
mm_xp = read.table(paste0(work_path, "/help_code/", "mm_xp.txt"), header=T, row.names=1, as.is=T)
zf_xp = read.table(paste0(work_path, "/help_code/", "zf_xp.txt"), header=T, row.names=1, as.is=T)

mm_cell_type_list = read.table(paste0(work_path, "/help_code/mouse_celltype_group.txt"), as.is=T, sep="\t")
mm_cell_type_list = as.vector(mm_cell_type_list$V1)

zf_cell_type_list = read.table(paste0(work_path, "/help_code/zf_celltype_list.txt"), as.is=T, sep="\t")
zf_cell_type_list = as.vector(zf_cell_type_list$V1)

xp_cell_type_list = read.table(paste0(work_path, "/help_code/xp_celltype_list.txt"), as.is=T, sep="\t")
xp_cell_type_list = as.vector(xp_cell_type_list$V1)

mm_keyTF = read.table(paste0(work_path, "/help_code/Table_S6.txt"), header=T, as.is=T, sep="\t")
zf_keyTF = read.table(paste0(work_path, "/help_code/Table_S14.txt"), header=T, as.is=T, sep="\t")
xp_keyTF = read.table(paste0(work_path, "/help_code/Table_S18.txt"), header=T, as.is=T, sep="\t")

mm_emergence_time = unique(mm_keyTF[,c("celltype", "emergence_time")])
zf_emergence_time = unique(zf_keyTF[,c("celltype", "emergence_time")])
xp_emergence_time = unique(xp_keyTF[,c("celltype", "emergence_time")])

print(mm_cell_type_list[!mm_cell_type_list %in% mm_emergence_time$celltype])
mm_emergence_time = rbind(mm_emergence_time, c("Inner cell mass", "3.5"))
mm_emergence_time = rbind(mm_emergence_time, c("Morula", "3"))
mm_emergence_time = rbind(mm_emergence_time, c("Trophectoderm", "3.5"))
mm_emergence_time = rbind(mm_emergence_time, c("Extraembryonic ectoderm", "5.25"))

print(zf_cell_type_list[!zf_cell_type_list %in% zf_emergence_time$celltype])
zf_emergence_time = rbind(zf_emergence_time, c("DEL", "hdf3.8"))
zf_emergence_time = rbind(zf_emergence_time, c("blastomere", "hpf3.3"))

print(xp_cell_type_list[!xp_cell_type_list %in% xp_emergence_time$celltype])
xp_emergence_time = rbind(xp_emergence_time, c("optic vesicle", "hpfS22"))
xp_emergence_time = rbind(xp_emergence_time, c("alpha ionocyte", "hpfS22"))
xp_emergence_time = rbind(xp_emergence_time, c("beta ionocyte", "hpfS22"))
xp_emergence_time = rbind(xp_emergence_time, c("blastula", "hpfS8"))
xp_emergence_time = rbind(xp_emergence_time, c("germ cell", "hpfS8"))
xp_emergence_time = rbind(xp_emergence_time, c("gut", "hpfS18"))
xp_emergence_time$emergence_time = gsub("hpf", "", xp_emergence_time$emergence_time)


#################
### mm vs. zf ###
#################

mm_cell_type_matrix = readRDS(paste0(work_path, "/", "mm_cell_type_matrix.rds"))
zf_cell_type_matrix = readRDS(paste0(work_path, "/", "zf_cell_type_matrix.rds"))

mm_cell_type_matrix = mm_cell_type_matrix[mm_zf$mm_gene,]
rownames(mm_cell_type_matrix) = rownames(mm_zf)
zf_cell_type_matrix = zf_cell_type_matrix[mm_zf$zf_gene,]
rownames(zf_cell_type_matrix) = rownames(mm_zf)

### exclude ExE-embryo cell types
mm_cell_type_matrix_sub = mm_cell_type_matrix[,!colnames(mm_cell_type_matrix) %in% c("Inner cell mass",
                                                                                     "Hypoblast",
                                                                                     "Parietal endoderm",
                                                                                     "Extraembryonic ectoderm",
                                                                                     "Visceral endoderm",
                                                                                     "Embryonic visceral endoderm",
                                                                                     "Extraembryonic visceral endoderm")]
zf_cell_type_matrix_sub = zf_cell_type_matrix[,!colnames(zf_cell_type_matrix) %in% c("blastomere","EVL","periderm","forerunner cells")]

conn <- correlation_analysis_bidirection(as.matrix(mm_cell_type_matrix_sub), as.matrix(zf_cell_type_matrix_sub), fold.change = 1.5, top_gene_num = 1200, spec_gene_num = 1200)
saveRDS(conn, paste0(work_path, "/mm_zf_cell_type.rds"))

conn = readRDS(paste0(work_path, "/mm_zf_cell_type.rds"))
conn$beta <- 2*(conn$beta_1+0.001)*(conn$beta_2+0.001)
write.csv(conn[conn$beta != 2e-6,], paste0(work_path, "/mm_zf_cell_type.csv"))
dat <- conn[,c("source","target","beta")]
dat <- dcast(dat, source~target)
rownames(dat) <- dat[,1]; dat <- dat[,-1]

# top5 vs. top5
x = conn %>% filter(beta > 1e-4) %>% group_by(source) %>% slice_max(order = beta, n = 5) %>% as.data.frame()
y = conn %>% filter(beta > 1e-4) %>% group_by(target) %>% slice_max(order = beta, n = 5) %>% as.data.frame()

res = NULL
x_celltype = as.vector(unique(x$source))
for(i in 1:length(x_celltype)){
  xx = as.vector(x$target[x$source==x_celltype[i]])
  for(j in 1:length(xx)){
    yy = as.vector(y$source[y$target == xx[j]])
    if (x_celltype[i] %in% yy){
      if(x_celltype[i] %in% mm_keyTF$celltype){
        emergence_time_A = paste0("E", mm_keyTF$emergence_time[mm_keyTF$celltype == x_celltype[i]][1])
        keyTF_A = as.vector(mm_keyTF$gene_ID[mm_keyTF$celltype == x_celltype[i]])
      } else {
        emergence_time_A = NA
        keyTF_A = NA
      }
      
      if(xx[j] %in% zf_keyTF$celltype){
        emergence_time_B = zf_keyTF$emergence_time[zf_keyTF$celltype == xx[j]][1]
        keyTF_B = as.vector(zf_keyTF$gene_ID[zf_keyTF$celltype == xx[j]])
      } else {
        emergence_time_B = NA
        keyTF_B = NA
      }
      
      share_keyTF = mm_zf_2[mm_zf_2$zf_gene %in% keyTF_B & mm_zf_2$mm_gene %in% keyTF_A,]
      A_overlap = sort(as.vector(unique(mm_keyTF$key_TF[mm_keyTF$gene_ID %in% as.vector(share_keyTF$mm_gene)])))
      B_overlap = sort(as.vector(unique(zf_keyTF$key_TF[zf_keyTF$gene_ID %in% as.vector(share_keyTF$zf_gene)])))
      
      A_specific = as.vector(unique(mm_keyTF$key_TF[mm_keyTF$gene_ID %in% keyTF_A]))
      A_specific = sort(A_specific[!A_specific %in% A_overlap])
      B_specific = as.vector(unique(zf_keyTF$key_TF[zf_keyTF$gene_ID %in% keyTF_B]))
      B_specific = sort(B_specific[!B_specific %in% B_overlap])
      
      tmp = data.frame(datasetA = x_celltype[i], 
                       datasetB = xx[j], 
                       beta = x$beta[x$source == x_celltype[i] & x$target == xx[j]], 
                       emergence_time_A = emergence_time_A, 
                       emergence_time_B = emergence_time_B, 
                       keyTF_overlap_A = paste(A_overlap, collapse = ", "),
                       keyTF_overlap_B = paste(B_overlap, collapse = ", "),
                       keyTF_specific_A = paste(A_specific, collapse = ", "),
                       keyTF_specific_B = paste(B_specific, collapse = ", "))
      res = rbind(res, tmp)
    }
  }
}
write.csv(res, paste0(work_path, "/mm_zf_cell_type_2.csv"),row.names = F)






#################
### mm vs. xp ###
#################

mm_cell_type_matrix = readRDS(paste0(work_path, "/", "mm_cell_type_matrix.rds"))
xp_cell_type_matrix = readRDS(paste0(work_path, "/", "xp_cell_type_matrix.rds"))

mm_cell_type_matrix = mm_cell_type_matrix[mm_xp$mm_gene,]
rownames(mm_cell_type_matrix) = rownames(mm_xp)
xp_cell_type_matrix = xp_cell_type_matrix[mm_xp$xp_gene,]
rownames(xp_cell_type_matrix) = rownames(mm_xp)

### exclude ExE-embryo cell types
mm_cell_type_matrix_sub = mm_cell_type_matrix[,!colnames(mm_cell_type_matrix) %in% c("Inner cell mass",
                                                                                     "Hypoblast",
                                                                                     "Parietal endoderm",
                                                                                     "Extraembryonic ectoderm",
                                                                                     "Visceral endoderm",
                                                                                     "Embryonic visceral endoderm",
                                                                                     "Extraembryonic visceral endoderm")]
xp_cell_type_matrix_sub = xp_cell_type_matrix

conn <- correlation_analysis_bidirection(as.matrix(mm_cell_type_matrix_sub), as.matrix(xp_cell_type_matrix_sub), fold.change = 1.5, top_gene_num = 1200, spec_gene_num = 1200)
saveRDS(conn, paste0(work_path, "/mm_xp_cell_type.rds"))

conn = readRDS(paste0(work_path, "/mm_xp_cell_type.rds"))
conn$beta <- 2*(conn$beta_1+0.001)*(conn$beta_2+0.001)
conn$comparing = "mm vs. xp"
write.csv(conn[conn$beta != 2e-6,], paste0(work_path, "/mm_xp_cell_type.csv"))

# top5 vs. top5
x = conn %>% filter(beta > 1e-4) %>% group_by(source) %>% slice_max(order = beta, n = 5) %>% as.data.frame()
y = conn %>% filter(beta > 1e-4) %>% group_by(target) %>% slice_max(order = beta, n = 5) %>% as.data.frame()

res = NULL
x_celltype = as.vector(unique(x$source))
for(i in 1:length(x_celltype)){
  xx = as.vector(x$target[x$source==x_celltype[i]])
  for(j in 1:length(xx)){
    yy = as.vector(y$source[y$target == xx[j]])
    if (x_celltype[i] %in% yy){
      if(x_celltype[i] %in% mm_keyTF$celltype){
        emergence_time_A = paste0("E", mm_keyTF$emergence_time[mm_keyTF$celltype == x_celltype[i]][1])
        keyTF_A = as.vector(mm_keyTF$gene_ID[mm_keyTF$celltype == x_celltype[i]])
      } else {
        emergence_time_A = NA
        keyTF_A = NA
      }
      
      if(xx[j] %in% xp_keyTF$celltype){
        emergence_time_B = xp_keyTF$emergence_time[xp_keyTF$celltype == xx[j]][1]
        keyTF_B = as.vector(xp_keyTF$key_TF[xp_keyTF$celltype == xx[j]])
      } else {
        emergence_time_B = NA
        keyTF_B = NA
      }
      
      share_keyTF = mm_xp_2[mm_xp_2$xp_gene %in% keyTF_B & mm_xp_2$mm_gene %in% keyTF_A,]
      A_overlap = sort(as.vector(unique(mm_keyTF$key_TF[mm_keyTF$gene_ID %in% as.vector(share_keyTF$mm_gene)])))
      B_overlap = sort(as.vector(unique(xp_keyTF$key_TF[xp_keyTF$key_TF %in% as.vector(share_keyTF$xp_gene)])))
      
      A_specific = as.vector(unique(mm_keyTF$key_TF[mm_keyTF$gene_ID %in% keyTF_A]))
      A_specific = sort(A_specific[!A_specific %in% A_overlap])
      B_specific = as.vector(unique(xp_keyTF$key_TF[xp_keyTF$key_TF %in% keyTF_B]))
      B_specific = sort(B_specific[!B_specific %in% B_overlap])
      
      tmp = data.frame(datasetA = x_celltype[i], 
                       datasetB = xx[j], 
                       beta = x$beta[x$source == x_celltype[i] & x$target == xx[j]], 
                       emergence_time_A = emergence_time_A, 
                       emergence_time_B = emergence_time_B, 
                       keyTF_overlap_A = paste(A_overlap, collapse = ", "),
                       keyTF_overlap_B = paste(B_overlap, collapse = ", "),
                       keyTF_specific_A = paste(A_specific, collapse = ", "),
                       keyTF_specific_B = paste(B_specific, collapse = ", "))
      res = rbind(res, tmp)
    }
  }
}
write.csv(res, paste0(work_path, "/mm_xp_cell_type_2.csv"),row.names = F)



#################
### zf vs. xp ###
#################

zf_cell_type_matrix = readRDS(paste0(work_path, "/", "zf_cell_type_matrix.rds"))
xp_cell_type_matrix = readRDS(paste0(work_path, "/", "xp_cell_type_matrix.rds"))

zf_cell_type_matrix = zf_cell_type_matrix[zf_xp$zf_gene,]
rownames(zf_cell_type_matrix) = rownames(zf_xp)
xp_cell_type_matrix = xp_cell_type_matrix[zf_xp$xp_gene,]
rownames(xp_cell_type_matrix) = rownames(zf_xp)

### exclude ExE-embryo cell types
xp_cell_type_matrix_sub = xp_cell_type_matrix
zf_cell_type_matrix_sub = zf_cell_type_matrix[,!colnames(zf_cell_type_matrix) %in% c("blastomere","EVL","periderm","forerunner cells")]

conn <- correlation_analysis_bidirection(as.matrix(zf_cell_type_matrix_sub), as.matrix(xp_cell_type_matrix_sub), fold.change = 1.5, top_gene_num = 1200, spec_gene_num = 1200)
saveRDS(conn, paste0(work_path, "/zf_xp_cell_type.rds"))

conn = readRDS(paste0(work_path, "/zf_xp_cell_type.rds"))
conn$beta <- 2*(conn$beta_1+0.001)*(conn$beta_2+0.001)
conn$comparing = "zf vs. xp"
write.csv(conn[conn$beta != 2e-6,], paste0(work_path, "/zf_xp_cell_type.csv"))

# top5 vs. top5
x = conn %>% filter(beta > 1e-4) %>% group_by(source) %>% slice_max(order = beta, n = 5) %>% as.data.frame()
y = conn %>% filter(beta > 1e-4) %>% group_by(target) %>% slice_max(order = beta, n = 5) %>% as.data.frame()

res = NULL
x_celltype = as.vector(unique(x$source))
for(i in 1:length(x_celltype)){
  xx = as.vector(x$target[x$source==x_celltype[i]])
  for(j in 1:length(xx)){
    yy = as.vector(y$source[y$target == xx[j]])
    if (x_celltype[i] %in% yy){
      if(x_celltype[i] %in% zf_keyTF$celltype){
        emergence_time_A = paste0("E", zf_keyTF$emergence_time[zf_keyTF$celltype == x_celltype[i]][1])
        keyTF_A = as.vector(zf_keyTF$gene_ID[zf_keyTF$celltype == x_celltype[i]])
      } else {
        emergence_time_A = NA
        keyTF_A = NA
      }
      
      if(xx[j] %in% xp_keyTF$celltype){
        emergence_time_B = xp_keyTF$emergence_time[xp_keyTF$celltype == xx[j]][1]
        keyTF_B = as.vector(xp_keyTF$key_TF[xp_keyTF$celltype == xx[j]])
      } else {
        emergence_time_B = NA
        keyTF_B = NA
      }
      
      share_keyTF = zf_xp[zf_xp$xp_gene %in% keyTF_B & zf_xp$zf_gene %in% keyTF_A,]
      A_overlap = sort(as.vector(unique(zf_keyTF$key_TF[zf_keyTF$gene_ID %in% as.vector(share_keyTF$zf_gene)])))
      B_overlap = sort(as.vector(unique(xp_keyTF$key_TF[xp_keyTF$key_TF %in% as.vector(share_keyTF$xp_gene)])))
      
      A_specific = as.vector(unique(zf_keyTF$key_TF[zf_keyTF$gene_ID %in% keyTF_A]))
      A_specific = sort(A_specific[!A_specific %in% A_overlap])
      B_specific = as.vector(unique(xp_keyTF$key_TF[xp_keyTF$key_TF %in% keyTF_B]))
      B_specific = sort(B_specific[!B_specific %in% B_overlap])
      
      tmp = data.frame(datasetA = x_celltype[i], 
                       datasetB = xx[j], 
                       beta = x$beta[x$source == x_celltype[i] & x$target == xx[j]], 
                       emergence_time_A = emergence_time_A, 
                       emergence_time_B = emergence_time_B, 
                       keyTF_overlap_A = paste(A_overlap, collapse = ", "),
                       keyTF_overlap_B = paste(B_overlap, collapse = ", "),
                       keyTF_specific_A = paste(A_specific, collapse = ", "),
                       keyTF_specific_B = paste(B_specific, collapse = ", "))
      res = rbind(res, tmp)
    }
  }
}
write.csv(res, paste0(work_path, "/zf_xp_cell_type_2.csv"),row.names = F)



