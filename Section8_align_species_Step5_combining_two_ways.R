################################
#### Section8_align_species ####
################################

#### Systematic comparison of the cellular trajectories of mouse, zebrafish and frog embryogenesis

##################################
#### Step5_combining_two_ways ####
##################################

library(tidyverse)
library(readxl)
library(dplyr)

work_path = ""

nnls = read_excel(paste0(work_path, "/help_code/Table.S21-23.xlsx"), sheet = "S21")
nnls = nnls[,c(1:4)]
names(nnls) = c("comparing", "datasetA", "datasetB", "comment")
nnls = nnls[nnls$comment == "Yes",]

mm_zf_nnls = nnls %>% filter(comparing == "Mouse vs. Zebrafish") %>% select("datasetA", "datasetB")
mm_xp_nnls = nnls %>% filter(comparing == "Mouse vs. Frog") %>% select("datasetA", "datasetB")
zf_xp_nnls = nnls %>% filter(comparing == "Zebrafish vs. Frog") %>% select("datasetA", "datasetB")

keyTF = read_excel(paste0(work_path, "/help_code/Table.S21-23.xlsx"), sheet = "S22")
keyTF = keyTF[,c(1:4)]
names(keyTF) = c("comparing", "datasetA", "datasetB", "comment")
keyTF = keyTF[keyTF$comment == "Yes",]

mm_zf_keyTF = keyTF %>% filter(comparing == "Mouse vs. Zebrafish") %>% select("datasetA", "datasetB")
mm_xp_keyTF = keyTF %>% filter(comparing == "Mouse vs. Frog") %>% select("datasetA", "datasetB")
zf_xp_keyTF = keyTF %>% filter(comparing == "Zebrafish vs. Frog") %>% select("datasetA", "datasetB")

mm_zf_nnls$method = "nnls"
mm_zf_keyTF$method = "keyTF"
mm_zf = rbind(mm_zf_nnls, mm_zf_keyTF)
mm_zf$ID = paste0(mm_zf$datasetA, "_", mm_zf$datasetB)
mm_zf_table = table(mm_zf$ID)
mm_zf_1 = mm_zf[mm_zf$ID %in% names(mm_zf_table)[mm_zf_table == 1],]
mm_zf_2 = mm_zf[mm_zf$ID %in% names(mm_zf_table)[mm_zf_table == 2],]
mm_zf_2 = unique(mm_zf_2[,c("datasetA", "datasetB")])
mm_zf_2$method = "nnls, keyTF"
mm_zf = rbind(mm_zf_1[,c("datasetA", "datasetB", "method")], mm_zf_2)
mm_zf$A_vs_B = "mouse vs. zebrafish"

mm_xp_nnls$method = "nnls"
mm_xp_keyTF$method = "keyTF"
mm_xp = rbind(mm_xp_nnls, mm_xp_keyTF)
mm_xp$ID = paste0(mm_xp$datasetA, "_", mm_xp$datasetB)
mm_xp_table = table(mm_xp$ID)
mm_xp_1 = mm_xp[mm_xp$ID %in% names(mm_xp_table)[mm_xp_table == 1],]
mm_xp_2 = mm_xp[mm_xp$ID %in% names(mm_xp_table)[mm_xp_table == 2],]
mm_xp_2 = unique(mm_xp_2[,c("datasetA", "datasetB")])
mm_xp_2$method = "nnls, keyTF"
mm_xp = rbind(mm_xp_1[,c("datasetA", "datasetB", "method")], mm_xp_2)
mm_xp$A_vs_B = "mouse vs. frog"

zf_xp_nnls$method = "nnls"
zf_xp_keyTF$method = "keyTF"
zf_xp = rbind(zf_xp_nnls, zf_xp_keyTF)
zf_xp$ID = paste0(zf_xp$datasetA, "_", zf_xp$datasetB)
zf_xp_table = table(zf_xp$ID)
zf_xp_1 = zf_xp[zf_xp$ID %in% names(zf_xp_table)[zf_xp_table == 1],]
zf_xp_2 = zf_xp[zf_xp$ID %in% names(zf_xp_table)[zf_xp_table == 2],]
zf_xp_2 = unique(zf_xp_2[,c("datasetA", "datasetB")])
zf_xp_2$method = "nnls, keyTF"
zf_xp = rbind(zf_xp_1[,c("datasetA", "datasetB", "method")], zf_xp_2)
zf_xp$A_vs_B = "zebrafish vs. frog"

dat = rbind(mm_zf, mm_xp, zf_xp)
dat$A_vs_B = factor(dat$A_vs_B, levels = c("mouse vs. zebrafish", "mouse vs. frog", "zebrafish vs. frog"))
dat = dat %>% arrange(A_vs_B, datasetA, method)
dat = dat[,c("A_vs_B", "datasetA", "datasetB", "method")]
names(dat) = c("A vs. B", "species A", "species B", "aligning method")
write.table(dat, paste0(work_path, "/corr_cell_types.txt"), row.names=F, quote=F, sep="\t")


###########################################
### make a plot ###########################
###########################################

# Libraries
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)

mm_zf$datasetA = paste0("mm:", mm_zf$datasetA)
mm_zf$datasetB = paste0("zf:", mm_zf$datasetB)
mm_xp$datasetA = paste0("mm:", mm_xp$datasetA)
mm_xp$datasetB = paste0("xp:", mm_xp$datasetB)
zf_xp$datasetA = paste0("zf:", zf_xp$datasetA)
zf_xp$datasetB = paste0("xp:", zf_xp$datasetB)

dat = rbind(mm_zf, mm_xp, zf_xp)
dat = dat[,c("datasetA", "datasetB", "method")]
names(dat) = c("A","B","method")

nodes = as.vector(unique(c(dat$A, dat$B)))
Nodes = data.frame(name = nodes,
                   group = unlist(lapply(nodes, function(x) strsplit(x,"[:]")[[1]][1])))
Nodes$ID = unlist(lapply(as.vector(Nodes$name), function(x) strsplit(x,"[:]")[[1]][2]))
print(table(Nodes$group))
Links = dat[,c("A","B","method")]

write.table(Nodes, "Nodes.txt", row.names=F, sep="\t", quote=F)
write.table(Links, "Links.txt", row.names=F, sep="\t", quote=F)

