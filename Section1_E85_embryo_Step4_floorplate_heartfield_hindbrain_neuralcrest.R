#############################
#### Section1_E85_embryo ####
#############################

#### The scripts used for analyzing E85 embryo data

###########################################################
#### Step4_floorplate_heartfield_hindbrain_neuralcrest ####
###########################################################

library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
source("help_code/help_code.R")

obj = readRDS("obj_E8.5b_processed.rds")
anno = readRDS("seurat_object_E8.5b.rds")
anno = data.frame(anno[[]])
anno$Anno = as.vector(anno$cell_state)
sum(rownames(anno) != colnames(obj))
obj$Anno = as.vector(anno$Anno)
obj$celltype = unlist(lapply(as.vector(anno$Anno), function(x) strsplit(x,"[:]")[[1]][2]))

#### add information for individual embryos with somite counts 
embryo_info = readRDS("help_code/embryo_info.rds")
anno = anno %>%
    left_join(embryo_info, by = "RT_group")
obj$somite_stage = as.vector(anno$somite_stage)
obj$somite_number = as.vector(anno$somite_number)
obj$embryo_id = as.vector(anno$embryo_id)
obj$embryo_sex = as.vector(anno$embryo_sex)

#### Part1: anterior and posterior floor plates

obj_sub = subset(obj, subset = celltype %in% c("Anterior floor plate", "Posterior floor plate"))
obj_sub = doClusterSeurat(obj_sub)
obj_sub$somite_number = as.numeric(obj_sub$somite_number)
obj_sub$somite_stage = factor(obj_sub$somite_stage, levels = paste0("stage_", c(0, 2:12)))

cbind(obj_sub[[]], Embeddings(obj_sub, reduction = "umap")) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = somite_stage)) +
    geom_point(size=1) + 
    theme_classic(base_size = 10) +
    scale_color_viridis_d() +
    theme(legend.position="none")

#### Part2: first and second heart field

obj_sub = subset(obj, subset = celltype %in% c("First heart field", "Second heart field"))
obj_sub = doClusterSeurat(obj_sub)
obj_sub$somite_number = as.numeric(obj_sub$somite_number)
obj_sub$somite_stage = factor(obj_sub$somite_stage, levels = paste0("stage_", c(0, 2:12)))

cbind(obj_sub[[]], Embeddings(obj_sub, reduction = "umap")) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = somite_stage)) +
    geom_point(size=1) + 
    theme_classic(base_size = 10) +
    scale_color_viridis_d() +
    theme(legend.position="none")

#### Part3: segmentations of hindbrain and neural crest

obj_sub = subset(obj, subset = celltype %in% c("Forebrain/midbrain",
                                               "Hindbrain",
                                               "Spinal cord",
                                               "Neural crest"))
obj_sub = doClusterSeurat(obj_sub)
obj_sub$somite_number = as.numeric(obj_sub$somite_number)
obj_sub$somite_stage = factor(obj_sub$somite_stage, levels = paste0("stage_", c(0, 2:12)))

cbind(obj_sub[[]], Embeddings(obj_sub, reduction = "umap")) %>%
    ggplot(aes(UMAP_1, UMAP_2, color = somite_stage)) +
    geom_point(size=1) + 
    theme_classic(base_size = 10) +
    scale_color_viridis_d() +
    theme(legend.position="none")







