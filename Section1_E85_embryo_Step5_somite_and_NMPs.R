#############################
#### Section1_E85_embryo ####
#############################

#### The scripts used for analyzing E85 embryo data

###############################
#### Step5_somite_and_NMPs ####
###############################

library(Seurat)
library(monocle3)
library(dplyr)
library(htmlwidgets)
library(plotly)
source("help_code/help_code.R")

obj = readRDS("obj_E8.5b_processed.rds")
anno = readRDS("obj_E8.5b_anno.rds")
sum(rownames(anno) != colnames(obj))
obj$Anno = as.vector(anno$Anno)
obj$celltype = unlist(lapply(as.vector(anno$Anno), function(x) strsplit(x,"[:]")[[1]][2]))

if(ncol(Embeddings(obj, reduction = "umap"))!=3){
    obj <- RunUMAP(object = obj, reduction = "pca", dims = 1:30, min.dist = 0.1, n.components = 3)
}

#### add information for individual embryos with somite counts 
embryo_info = readRDS("help_code/embryo_info.rds")
anno = anno %>%
    left_join(embryo_info, by = "RT_group")
obj$somite_stage = as.vector(anno$somite_stage)
obj$somite_number = as.vector(anno$somite_number)
obj$embryo_id = as.vector(anno$embryo_id)
obj$embryo_sex = as.vector(anno$embryo_sex)

pd = data.frame(obj[[]], Embeddings(obj, reduction = "umap"))
somite_stage_level = paste0("stage_", c(0, 2:12))
pd$somite_stage = factor(pd$somite_stage, levels = somite_stage_level)
pd$somite_number = as.numeric(pd$somite_number)

### 3D UMAP, with coloring cells by somite counts
fig = plot_ly(pd, x=~UMAP_1, y=~UMAP_2, z=~UMAP_3, size = I(30), color = ~somite_number, mode = "markers", marker = list(colors='Viridis'))
fig = fig %>% layout(
    scene = list(xaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 yaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE),
                 zaxis = list(title = '', autorange = TRUE, showgrid = FALSE, zeroline = FALSE, showline = FALSE, autotick = TRUE, ticks = '', showticklabels = FALSE)))


#### performing Knn followed by averaging the somite counts of the nearest 5 cells 

library(FNN)

k_neigh = 5
emb = as.matrix(pd[,c("UMAP_1", "UMAP_2", "UMAP_3")])
neighbors <- get.knnx(emb, emb, k = k_neigh + 1)$nn.index

tmp = matrix(NA, nrow(emb), k_neigh)
for(i in 1:k_neigh){
    tmp[,i] = pd$somite_number[neighbors[,(i+1)]]
}

pd$somite_number_smooth = apply(tmp, 1, mean)

celltype_list = table(pd$celltype)
celltype_list = names(celltype_list)[celltype_list > 100]
cor_res = rep(NA, length(celltype_list))
for(i in 1:length(celltype_list)){
    pd_sub = pd %>% filter(celltype == celltype_list[i])
    fit = cor.test(pd_sub$somite_number_smooth, pd_sub$somite_number)
    cor_res[i] = fit$estimate
}

df = data.frame(celltype = celltype_list, corr = cor_res) %>%
    arrange(desc(corr))
df$celltype = factor(df$celltype, levels = as.vector(df$celltype))
df %>%
    ggplot(aes(x=celltype, y=corr)) +
    geom_bar(stat="identity") +
    theme_classic(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    theme(legend.position="none") 

#### focusing on NMPs 

obj_sub = subset(obj, subset = celltype %in% c("Neuromesodermal progenitors"))
obj_sub = doClusterSeurat(obj_sub)
obj_sub$somite_number = as.numeric(obj_sub$somite_number)
obj_sub$somite_stage = factor(obj_sub$somite_stage, levels = paste0("stage_", c(0, 2:12)))

#### Visualization of PCA 

pca <- obj_sub[["pca"]]
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / sum(eigValues)
print(varExplained)

p1 = cbind(obj_sub[[]], Embeddings(obj_sub, reduction = "pca")[,c(1,2)]) %>%
    ggplot(aes(PC_1, PC_2, color = somite_stage)) +
    geom_point(size=0.5) + 
    theme_classic(base_size = 10) +
    scale_color_viridis_d() +
    theme(legend.position="none") +
    labs(x = "PC 1 (23.7%)", y = "PC 2 (15.1%)") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

p2 = cbind(obj_sub[[]], Embeddings(obj_sub, reduction = "pca")[,c(1,3)]) %>%
    ggplot(aes(PC_1, PC_3, color = somite_stage)) +
    geom_point(size=0.5) + 
    theme_classic(base_size = 10) +
    scale_color_viridis_d() +
    theme(legend.position="none") +
    labs(x = "PC 1 (23.7%)", y = "PC 3 (8.4%)") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

p3 = cbind(obj_sub[[]], Embeddings(obj_sub, reduction = "pca")[,c(2,3)]) %>%
    ggplot(aes(PC_2, PC_3, color = somite_stage)) +
    geom_point(size=0.5) + 
    theme_classic(base_size = 10) +
    scale_color_viridis_d() +
    theme(legend.position="none") +
    labs(x = "PC 2 (15.1%)", y = "PC 3 (8.4%)") +
    theme(axis.text.x = element_text(color="black"), axis.text.y = element_text(color="black"))

df = cbind(obj_sub[[]], Embeddings(obj_sub, reduction = "pca")[,1:3])
fig = plot_ly(df, x=~PC_1, y=~PC_2, z=~PC_3, size = I(30), color = ~somite_number, mode = "markers", marker = list(colors='Viridis')) %>%
    layout(scene =list(xaxis = list(title = "PC 1 (23.7%)"),
                       yaxis = list(title = "PC 2 (15.1%)"),
                       zaxis = list(title = "PC 3 (8.4%)")))
