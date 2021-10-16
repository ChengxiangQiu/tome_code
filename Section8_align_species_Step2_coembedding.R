################################
#### Section8_align_species ####
################################

#### Systematic comparison of the cellular trajectories of mouse, zebrafish and frog embryogenesis

###########################
#### Step2_coembedding ####
###########################

library(Seurat)
library(dplyr)
work_path = ""
source("help_code/help_code.R")

gene_overlap = read.table("helpc_code/mm_zf_xp.txt", header=T, row.names=1, as.is=T)

mm_cell_type_list = read.table("helpc_code/mouse_celltype_group.txt", as.is=T, sep="\t")
mm_cell_type_list = as.vector(mm_cell_type_list$V1)

zf_cell_type_list = read.table("help_code/zf_celltype_list.txt", as.is=T, sep="\t")
zf_cell_type_list = as.vector(zf_cell_type_list$V1)

xp_cell_type_list = read.table("help_code/xp_celltype_list.txt", as.is=T, sep="\t")
xp_cell_type_list = as.vector(xp_cell_type_list$V1)


time_point = paste0("E", c(3.5,4.5,5.25,5.5,6.25,seq(6.5,8.25,0.25),"8.5a","8.5b",seq(9.5,13.5)))
zf_time_point = paste0("hpf", c(3.3, 3.8, 4.3, 4.8, 5.3, 6, 7, 8, 9, 10, 11, 12, 14, 18, 24))
xp_time_point = paste0("S", c(8,10,11,12,13,14,16,18,20,22))

mm_cell_num = readRDS(paste0(work_path, "/", "mm_cell_num.rds"))
mm_cell_num$time = unlist(lapply(as.vector(mm_cell_num$Anno), function(x) strsplit(x,"[:]")[[1]][1]))
mm_cell_num$species = "mm"
mm_cell_num$dataset = "mm_1"
mm_cell_num$dataset[mm_cell_num$time %in% paste0("E", "8.5b")] = "mm_2"
mm_cell_num$dataset[mm_cell_num$time %in% paste0("E",seq(9.5,13.5))] = "mm_2"
mm_cell_num = mm_cell_num[!mm_cell_num$celltype %in% c("Inner cell mass",
                                                       "Hypoblast",
                                                       "Parietal endoderm",
                                                       "Extraembryonic ectoderm",
                                                       "Visceral endoderm",
                                                       "Embryonic visceral endoderm",
                                                       "Extraembryonic visceral endoderm"),]

zf_cell_num = readRDS(paste0(work_path, "/", "zf_cell_num.rds"))
zf_cell_num$time = unlist(lapply(as.vector(zf_cell_num$Anno), function(x) strsplit(x,"[:]")[[1]][1]))
zf_cell_num$species = "zf"
zf_cell_num$dataset = "zf"
zf_cell_num = zf_cell_num[!zf_cell_num$celltype %in% c("blastomere","EVL","periderm","forerunner cells"),]

xp_cell_num = readRDS(paste0(work_path, "/", "xp_cell_num.rds"))
xp_cell_num$time = unlist(lapply(as.vector(xp_cell_num$Anno), function(x) strsplit(x,"[:]")[[1]][1]))
xp_cell_num$species = "xp"
xp_cell_num$dataset = "xp"

mm_cell_state_matrix = readRDS(paste0(work_path, "/", "mm_cell_state_matrix.rds"))
gene_overlap_sub = gene_overlap[!is.na(gene_overlap$mm_gene),]
mm_cell_state_matrix = mm_cell_state_matrix[as.vector(gene_overlap_sub$mm_gene), as.vector(mm_cell_num$Anno)]
rownames(mm_cell_state_matrix) = rownames(gene_overlap_sub)

zf_cell_state_matrix = readRDS(paste0(work_path, "/", "zf_cell_state_matrix.rds"))
gene_overlap_sub = gene_overlap[!is.na(gene_overlap$zf_gene),]
zf_cell_state_matrix = zf_cell_state_matrix[as.vector(gene_overlap_sub$zf_gene), as.vector(zf_cell_num$Anno)]
rownames(zf_cell_state_matrix) = rownames(gene_overlap_sub)

xp_cell_state_matrix = readRDS(paste0(work_path, "/", "xp_cell_state_matrix.rds"))
gene_overlap_sub = gene_overlap[!is.na(gene_overlap$xp_gene),]
xp_cell_state_matrix = xp_cell_state_matrix[as.vector(gene_overlap_sub$xp_gene), as.vector(xp_cell_num$Anno)]
rownames(xp_cell_state_matrix) = rownames(gene_overlap_sub)

### 
obj_mm = CreateSeuratObject(as.matrix(mm_cell_state_matrix), meta.data = mm_cell_num)
obj_zf = CreateSeuratObject(as.matrix(zf_cell_state_matrix), meta.data = zf_cell_num)
obj_xp = CreateSeuratObject(as.matrix(xp_cell_state_matrix), meta.data = xp_cell_num)

obj = merge(obj_mm, list(obj_zf, obj_xp))
obj$group = obj$dataset
table(obj$group)
count = GetAssayData(object = obj, slot = "counts")

#obj = doClusterSeurat(obj, k.filter = 100)

obj.list <- SplitObject(object = obj, split.by = "group")

for (i in 1:length(x = obj.list)) {
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
  obj.list[[i]] <- SetAssayData(object = obj.list[[i]], slot = "data", new.data = count[,obj$group==names(obj.list)[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], 
                                        selection.method = "vst", nfeatures = 5000, verbose = FALSE)
}

reference.list <- obj.list[names(table(obj$group))]
obj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30, k.filter = 100)

obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:30)

DefaultAssay(object = obj.integrated) <- "integrated"

obj <- obj.integrated 
obj <- ScaleData(object = obj, verbose = FALSE)
obj <- RunPCA(object = obj, npcs = 30, verbose = FALSE)
obj <- RunUMAP(object = obj, reduction = "pca", dims = 1:30, min.dist = 0.6)
obj <- RunTSNE(object = obj, reduction = "pca", dims = 1:30)
obj <- FindNeighbors(object = obj, reduction = "umap", dims = 1:2)
obj <- FindClusters(object = obj, resolution = 10)
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend()
DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1, group.by = "group") + NoLegend()

### make the scatter plot ###
df = data.frame(obj[[]])
df$umap_1 = Embeddings(object = obj, reduction = "umap")[,1]
df$umap_2 = Embeddings(object = obj, reduction = "umap")[,2]
df$pca_1 = Embeddings(object = obj, reduction = "pca")[,1]
df$pca_2 = Embeddings(object = obj, reduction = "pca")[,2]

mm_cell_type = mm_cell_type_list[mm_cell_type_list %in% as.vector(df$celltype[df$species=="mm"])]
mm_cell_type = paste0("mm:",mm_cell_type)
zf_cell_type = zf_cell_type_list[zf_cell_type_list %in% as.vector(df$celltype[df$species=="zf"])]
zf_cell_type = paste0("zf:",zf_cell_type)
xp_cell_type = xp_cell_type_list[xp_cell_type_list %in% as.vector(df$celltype[df$species=="xp"])]
xp_cell_type = paste0("xp:",xp_cell_type)
label = data.frame(celltype = c(mm_cell_type, zf_cell_type, xp_cell_type), 
                   id = c(1:length(mm_cell_type), 1:length(zf_cell_type), 1:length(xp_cell_type)))
label$xx = unlist(lapply(as.vector(label$celltype), function(x) strsplit(x,"[:]")[[1]][2]))
rownames(label) = as.vector(label$celltype)
df$species_celltype = paste0(df$species, ":", df$celltype)
label_sub = label[as.vector(df$species_celltype),]
df$label = label_sub$id

library(ggplot2)
cols <- c("mm" = "#FF9999", "zf" = "lightblue", "xp" = "lightgreen")
ggplot() +
  geom_point(data=df, aes(x=umap_1, y=umap_2, fill=species), size = 3, shape=21, color="black") +
  geom_text(data=df, aes(x=umap_1, y=umap_2, label=label), size = 2) +
  theme_classic(base_size = 10) +
  scale_fill_manual(values = cols)

ggplot() +
  geom_point(data=df, aes(x=pca_1, y=pca_2, fill=species), size = 3.5, shape=21, color="black") +
  geom_text(data=df, aes(x=pca_1, y=pca_2, label=label), size = 2.5) +
  theme_classic(base_size = 10)


