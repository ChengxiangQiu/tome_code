################################
#### Section8_align_species ####
################################

#### Systematic comparison of the cellular trajectories of mouse, zebrafish and frog embryogenesis

##################################
#### Step4_overlapping_keyTFs ####
##################################

rm(list=ls())
library(dplyr)
library(reshape2)
library("gplots")
work_path = ""

zf_xp = read.table(paste0(work_path, "/help_code/", "zf_xp.txt"), header=T, row.names=1, as.is=T)
mm_zf_2 = read.table(paste0(work_path, "/help_code/", "mm_zf_2.txt"), header=T, row.names=1, as.is=T)
mm_xp_2 = read.table(paste0(work_path, "/help_code/", "mm_xp_2.txt"), header=T, row.names=1, as.is=T)
mm_zf = mm_zf_2
mm_xp = mm_xp_2

mm = read.table(paste0(work_path, "/help_code/Table_S6.txt"), header=T, as.is=T, sep="\t")
mm$ID = as.vector(mm$gene_ID)
zf = read.table(paste0(work_path, "/help_code/Table_S14.txt"), header=T, as.is=T, sep="\t")
zf$ID = as.vector(zf$gene_ID)
xp = read.table(paste0(work_path, "/help_code/Table_S18.txt"), header=T, as.is=T, sep="\t")
xp$ID = as.vector(xp$key_TF)

mm = unique(mm[,c("celltype", "ID")])
zf = unique(zf[,c("celltype", "ID")])
xp = unique(xp[,c("celltype", "ID")])

mm = mm[!mm$celltype %in% c("Inner cell mass",
                            "Hypoblast",
                            "Parietal endoderm",
                            "Extraembryonic ectoderm",
                            "Visceral endoderm",
                            "Embryonic visceral endoderm",
                            "Extraembryonic visceral endoderm"),]
zf = zf[!zf$celltype %in% c("blastomere","EVL","periderm","forerunner cells"),]

mm_celltype = as.vector(unique(mm$celltype))
zf_celltype = as.vector(unique(zf$celltype))
xp_celltype = as.vector(unique(xp$celltype))

mm_keyTF = as.vector(unique(mm$ID))
zf_keyTF = as.vector(unique(zf$ID))
xp_keyTF = as.vector(unique(xp$ID))

mm_keyTF_num = table(mm$celltype)
zf_keyTF_num = table(zf$celltype)
xp_keyTF_num = table(xp$celltype)

### obs of mm vs. zf ###################################################################################
mm_TF_score = table(as.vector(mm$ID))/length(mm_celltype)
zf_TF_score = table(as.vector(zf$ID))/length(zf_celltype)

obs = matrix(NA, length(mm_celltype), length(zf_celltype))

for(i in 1:length(mm_celltype)){
  for(j in 1:length(zf_celltype)){
    
    x1 = as.vector(mm$ID[mm$celltype == mm_celltype[i]])
    x2 = as.vector(zf$ID[zf$celltype == zf_celltype[j]])
    
    overlap = mm_zf[mm_zf$mm_gene %in% x1 & mm_zf$zf_gene %in%x2,]
    if(nrow(overlap)>0){
      y = mm_TF_score[overlap$mm_gene] * zf_TF_score[overlap$zf_gene]
      obs[i, j] = prod(y) * length(x1) * length(x2)
    }
  }
}


######### permutation ##################
perm = list()
perm_num = 10000

for(cnt in 1:perm_num){
  
  print(paste0(cnt, "/", perm_num))
  
  mm_sim = data.frame(celltype = rep(names(mm_keyTF_num), times = as.vector(mm_keyTF_num)),
                      ID = unlist(lapply(as.vector(mm_keyTF_num), function(x) sample(mm_keyTF, x, replace = F))))
  zf_sim = data.frame(celltype = rep(names(zf_keyTF_num), times = as.vector(zf_keyTF_num)),
                      ID = unlist(lapply(as.vector(zf_keyTF_num), function(x) sample(zf_keyTF, x, replace = F))))
  
  
  mm_TF_score_sim = table(as.vector(mm_sim$ID))/length(mm_celltype)
  zf_TF_score_sim = table(as.vector(zf_sim$ID))/length(zf_celltype)
  
  perm_cnt = matrix(NA, length(mm_celltype), length(zf_celltype))
  
  for(i in 1:length(mm_celltype)){
    for(j in 1:length(zf_celltype)){
      
      x1 = as.vector(mm_sim$ID[mm_sim$celltype == mm_celltype[i]])
      x2 = as.vector(zf_sim$ID[zf_sim$celltype == zf_celltype[j]])
      
      overlap = mm_zf[mm_zf$mm_gene %in% x1 & mm_zf$zf_gene %in%x2,]
      if(nrow(overlap)>0){
        y = mm_TF_score[overlap$mm_gene] * zf_TF_score[overlap$zf_gene]
        perm_cnt[i, j] = prod(y) * length(x1) * length(x2)
      }
    }
  }
  
  perm[[cnt]] = perm_cnt
  
}

saveRDS(list(obs = obs, perm = perm), paste0(work_path, "/mm_zf_perm.rds"))


mm_zf_obs = obs
mm_zf_perm = perm
mm_zf_obs_sig = matrix(NA, length(mm_celltype), length(zf_celltype))
mm_zf_perm_num = matrix(NA, length(mm_celltype), length(zf_celltype))

for(i in 1:length(mm_celltype)){
  for(j in 1:length(zf_celltype)){
    print(paste0(i, "/", j))
    if(!is.na(mm_zf_obs[i, j])){
      tmp = unlist(lapply(mm_zf_perm, function(x) x[i, j]))
      tmp = tmp[!is.na(tmp)]
      mm_zf_perm_num[i, j] = length(tmp)
      mm_zf_obs_sig[i, j] = sum(tmp < mm_zf_obs[i, j])
    }
  }
}

rownames(mm_zf_obs_sig) = rownames(mm_zf_perm_num) = mm_celltype
colnames(mm_zf_obs_sig) = colnames(mm_zf_perm_num) = zf_celltype

mm_zf_obs_sig = melt(mm_zf_obs_sig)
mm_zf_perm_num = melt(mm_zf_perm_num)

sum(mm_zf_obs_sig$Var1 != mm_zf_perm_num$Var1)
sum(mm_zf_obs_sig$Var2 != mm_zf_perm_num$Var2)

dat = data.frame(mm = mm_zf_obs_sig$Var1,
                 zf = mm_zf_obs_sig$Var2,
                 obs = mm_zf_obs_sig$value,
                 perm = mm_zf_perm_num$value)

dat$mm = as.vector(dat$mm)
dat$zf = as.vector(dat$zf)

dat = dat[!is.na(dat$perm) & dat$perm > 1000, ]
dat$sig = dat$obs / dat$perm
dat$sig_fdr = p.adjust(dat$sig, method = "fdr", n = nrow(dat))

result = dat %>% filter(sig < 0.01) %>% arrange(mm, sig) %>% as.data.frame()

tmp = NULL
for(i in 1:nrow(result)){
  if(result$mm[i] %in% mm$celltype){
    keyTF_A = as.vector(mm$ID[mm$celltype == result$mm[i]])
  } else {
    keyTF_A = NA
  }
  
  if(result$zf[i] %in% zf$celltype){
    keyTF_B = as.vector(zf$ID[zf$celltype == result$zf[i]])
  } else {
    keyTF_B = NA
  }
  
  share_keyTF = mm_zf_convert[mm_zf_convert$zf_gene %in% keyTF_B & mm_zf_convert$mm_gene %in% keyTF_A,]
  A_overlap = sort(as.vector(unique(mm_keyTF$key_TF[mm_keyTF$gene_ID %in% as.vector(share_keyTF$mm_gene)])))
  B_overlap = sort(as.vector(unique(zf_keyTF$key_TF[zf_keyTF$gene_ID %in% as.vector(share_keyTF$zf_gene)])))
  
  A_specific = as.vector(unique(mm_keyTF$key_TF[mm_keyTF$gene_ID %in% keyTF_A]))
  A_specific = sort(A_specific[!A_specific %in% A_overlap])
  B_specific = as.vector(unique(zf_keyTF$key_TF[zf_keyTF$gene_ID %in% keyTF_B]))
  B_specific = sort(B_specific[!B_specific %in% B_overlap])
  
  if (result$mm[i] %in% mm_emergence_time$celltype){
    time_A = as.vector(mm_emergence_time$emergence_time[mm_emergence_time$celltype == result$mm[i]])
  } else {
    time_A = NA
  }
  
  if (result$zf[i] %in% zf_emergence_time$celltype){
    time_B = as.vector(zf_emergence_time$emergence_time[zf_emergence_time$celltype == result$zf[i]])
  } else {
    time_B = NA
  }
  
  tmp = rbind(tmp, data.frame(time_A = time_A,
                              time_B = time_B,
                              keyTF_overlap_A = paste(A_overlap, collapse = ", "),
                              keyTF_overlap_B = paste(B_overlap, collapse = ", "),
                              keyTF_specific_A = paste(A_specific, collapse = ", "),
                              keyTF_specific_B = paste(B_specific, collapse = ", "))
  )
}

result = cbind(result, tmp)
write.table(result, paste0(work_path, "mm_zf_align_keyTF.txt"), row.names=F, quote=F, sep="\t")






