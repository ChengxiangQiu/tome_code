#####################################
#### Section4_spatial_inferrence ####
#####################################

#### Inference of the approximate spatial locations of cell states during mouse gastrulation

##################################
#### Step1_prepare_CIBETSORTx ####
##################################

#### GEO-seq data was downloaded from the paper
#### https://www.nature.com/articles/s41586-019-1469-8

dat = NULL

E60 = read.table("E6.0.txt", header=T, sep="\t", as.is=T, row.names=1)
gene = rownames(E60)
dat = cbind(dat, as.matrix(E60))

E55 = read.table("E5.5.txt", header=T, sep="\t", as.is=T)
dat = cbind(dat, as.matrix(E55[,-1]))

E65 = read.table("E6.5.txt", header=T, sep="\t", as.is=T)
dat = cbind(dat, as.matrix(E65[,-1]))

E70 = read.table("E7.0.txt", header=T, sep="\t", as.is=T)
dat = cbind(dat, as.matrix(E70[,-1]))

E75 = read.table("E7.5.txt", header=T, sep="\t", as.is=T)
dat = cbind(dat, as.matrix(E75[,-1]))

#### It has been confirmed by the author, that there is one data point missing
dat = data.frame(dat)
rownames(dat) = gene
dat$X3A.E7.0 = (dat$X2A.E7.0 + dat$X4A.E7.0)/2

x = data.frame(ID = names(dat), name = gsub("X", "", names(dat)))
x$ID = as.vector(x$ID); x$name = as.vector(x$name)
x$name[29:58] = paste0(x$name[29:58], ".E6.5")

corn = read.table("/CornPlot/example.txt", row.names=1, header=T)
y = data.frame(ID = names(corn), name = gsub("X", "", names(corn)))
y$ID = as.vector(y$ID); y$name = as.vector(y$name)

x$name[19:28] = y$name[1:10]

sum(!x$name %in% y$name)

rownames(x) = as.vector(x$name)
x = x[as.vector(y$name),]
sum(rownames(x) != y$name)

y$add = as.vector(x$ID)
dat = dat[,as.vector(y$add)]
names(dat) = as.vector(y$name)

y$time = c(rep("E5.5",10), rep("E6",18), rep("E6.5",30), rep("E7",73), rep("E7.5",83))
dat[1:5,y$time=="E5.5"][,1:3]
y = y[,c("ID","name","time")]

saveRDS(list(exp = dat, pd = y), "dat.rds")



#### prepare FPKMs for single-cell reference signatures 


rm(list=ls())
library(Seurat)
library(DESeq2)          ### used for normalization
library(GenomicFeatures) ### used for annotation of genes
library(dplyr)
source("help_code/help_code.R")

mouse_gene = read.table("help_code/mouse.v12.geneID.txt", header = T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)


#### E5.5 

time_i = "E5.5" 
obj = readRDS(paste0("obj_", time_i, "_exp.rds"))
anno = readRDS(paste0("obj_", time_i, "_anno.rds"))
anno = anno[colnames(obj),]
anno$celltype = unlist(lapply(as.vector(anno$Anno), function(x) strsplit(x,"[:]")[[1]][2]))
Idents(obj) = as.vector(anno$celltype)

#### only keey 352 dataset (Cheng et al.)

obj$keep = anno$group != 351

obj = subset(obj, keep)

count = GetAssayData(object = obj, slot = "counts")
anno = anno[colnames(obj), c("celltype", "group")]
anno$celltype = gsub(" ", "_", anno$celltype)
table(anno$celltype)
sum(rownames(anno) != colnames(count))

mouse_gene_sub = mouse_gene[rownames(mouse_gene) %in% rownames(count),]
count = count[rownames(mouse_gene_sub),]
rownames(count) = as.vector(mouse_gene_sub$ID)

### Generate DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = anno,
                              design = ~ celltype)

rowRanges(dds) <- ebg[rownames(count)]

#### obtain normalized FPKM value, and output this dataset for further use

fpkm.size <- fpkm(dds, robust = TRUE)

mouse_gene_sub$rowSums = Matrix::rowSums(fpkm.size)
mouse_gene_sub = mouse_gene_sub[mouse_gene_sub$gene_short_name %in% rownames(exp),]

#### exclude duplicated names 

tmp = mouse_gene_sub %>% group_by(gene_short_name) %>% slice_max(order=rowSums, n=1) %>% as.data.frame()
x = table(tmp$gene_short_name)
tmp[tmp$gene_short_name %in% rownames(x)[x>1],]
tmp = tmp[!tmp$gene_short_name %in% rownames(x)[x>1],]

fpkm = fpkm.size[as.vector(tmp$ID),]
rownames(fpkm) = as.vector(tmp$gene_short_name)
colnames(fpkm) = as.vector(anno$celltype)

fpkm = round(fpkm, digits = 2)
write.table(fpkm, "E55_sig.txt", quote=F, sep="\t")

sum(!rownames(fpkm) %in% rownames(exp))
exp_mix = exp[rownames(fpkm), pd$time == time_i]
write.table(exp_mix, "E55_mix.txt", quote=F, sep="\t")




#### E6.25

time_i = "E6.25"
obj = readRDS(paste0("obj_", time_i, "_exp.rds"))
anno = readRDS(paste0("obj_", time_i, "_anno.rds"))
anno = anno[colnames(obj),]
anno$group = anno$embryo
anno$celltype = unlist(lapply(as.vector(anno$Anno), function(x) strsplit(x,"[:]")[[1]][2]))
Idents(obj) = as.vector(anno$celltype)

count = GetAssayData(object = obj, slot = "counts")
anno = anno[colnames(obj), c("celltype", "group")]
anno$celltype = gsub(" ", "_", anno$celltype)
table(anno$celltype)
sum(rownames(anno) != colnames(count))

mouse_gene_sub = mouse_gene[rownames(mouse_gene) %in% rownames(count),]
count = count[rownames(mouse_gene_sub),]
rownames(count) = as.vector(mouse_gene_sub$ID)

### Generate DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = anno,
                              design = ~ celltype)

rowRanges(dds) <- ebg[rownames(count)]

### obtain normalized FPKM value, and output this dataset for further use
fpkm.size <- fpkm(dds, robust = TRUE)

mouse_gene_sub$rowSums = Matrix::rowSums(fpkm.size)
mouse_gene_sub = mouse_gene_sub[mouse_gene_sub$gene_short_name %in% rownames(exp),]

### exclude duplicated names 
tmp = mouse_gene_sub %>% group_by(gene_short_name) %>% slice_max(order=rowSums, n=1) %>% as.data.frame()
x = table(tmp$gene_short_name)
tmp[tmp$gene_short_name %in% rownames(x)[x>1],]
tmp = tmp[!tmp$gene_short_name %in% rownames(x)[x>1],]

fpkm = fpkm.size[as.vector(tmp$ID),]
rownames(fpkm) = as.vector(tmp$gene_short_name)
colnames(fpkm) = as.vector(anno$celltype)

fpkm = round(fpkm, digits = 2)
write.table(fpkm, "E60_sig.txt", quote=F, sep="\t")

sum(!rownames(fpkm) %in% rownames(exp))
exp_mix = exp[rownames(fpkm), pd$time == "E6"]
write.table(exp_mix, "E60_mix.txt", quote=F, sep="\t")





#### E6.5

time_i = "E6.5"
obj = readRDS(paste0("obj_", time_i, "_exp.rds"))
anno = readRDS(paste0("obj_", time_i, "_anno.rds"))
anno = anno[colnames(obj),]
anno$celltype = unlist(lapply(as.vector(anno$Anno), function(x) strsplit(x,"[:]")[[1]][2]))
Idents(obj) = as.vector(anno$celltype)

obj$keep = !anno$group %in% c(351, 352)

obj = subset(obj, keep)

count = GetAssayData(object = obj, slot = "counts")
anno = anno[colnames(obj), c("celltype", "group")]
anno$celltype = gsub(" ", "_", anno$celltype)
table(anno$celltype)
sum(rownames(anno) != colnames(count))

mouse_gene_sub = mouse_gene[rownames(mouse_gene) %in% rownames(count),]
count = count[rownames(mouse_gene_sub),]
rownames(count) = as.vector(mouse_gene_sub$ID)

### Generate DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = anno,
                              design = ~ celltype)

rowRanges(dds) <- ebg[rownames(count)]

### obtain normalized FPKM value, and output this dataset for further use
fpkm.size <- fpkm(dds, robust = TRUE)

mouse_gene_sub$rowSums = Matrix::rowSums(fpkm.size)
mouse_gene_sub = mouse_gene_sub[mouse_gene_sub$gene_short_name %in% rownames(exp),]

### exclude duplicated names 
tmp = mouse_gene_sub %>% group_by(gene_short_name) %>% slice_max(order=rowSums, n=1) %>% as.data.frame()
x = table(tmp$gene_short_name)
tmp[tmp$gene_short_name %in% rownames(x)[x>1],]
tmp = tmp[!tmp$gene_short_name %in% rownames(x)[x>1],]

fpkm = fpkm.size[as.vector(tmp$ID),]
rownames(fpkm) = as.vector(tmp$gene_short_name)
colnames(fpkm) = as.vector(anno$celltype)

fpkm = round(fpkm, digits = 2)
write.table(fpkm, "E65_sig.txt", quote=F, sep="\t")

sum(!rownames(fpkm) %in% rownames(exp))
exp_mix = exp[rownames(fpkm), pd$time == "E6.5"]
write.table(exp_mix, "E65_mix.txt", quote=F, sep="\t")





#### E7

time_i = "E7"
obj = readRDS(paste0("obj_", time_i, "_exp.rds"))
anno = readRDS(paste0("obj_", time_i, "_anno.rds"))
anno = anno[colnames(obj),]
anno$celltype = unlist(lapply(as.vector(anno$Anno), function(x) strsplit(x,"[:]")[[1]][2]))
Idents(obj) = as.vector(anno$celltype)


count = GetAssayData(object = obj, slot = "counts")
anno = anno[colnames(obj), c("celltype", "group")]
anno$celltype = gsub(" ", "_", anno$celltype)
table(anno$celltype)
sum(rownames(anno) != colnames(count))

mouse_gene_sub = mouse_gene[rownames(mouse_gene) %in% rownames(count),]
count = count[rownames(mouse_gene_sub),]
rownames(count) = as.vector(mouse_gene_sub$ID)

### Generate DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = anno,
                              design = ~ celltype)

rowRanges(dds) <- ebg[rownames(count)]

### obtain normalized FPKM value, and output this dataset for further use
fpkm.size <- fpkm(dds, robust = TRUE)

mouse_gene_sub$rowSums = Matrix::rowSums(fpkm.size)
mouse_gene_sub = mouse_gene_sub[mouse_gene_sub$gene_short_name %in% rownames(exp),]

### exclude duplicated names 
tmp = mouse_gene_sub %>% group_by(gene_short_name) %>% slice_max(order=rowSums, n=1) %>% as.data.frame()
x = table(tmp$gene_short_name)
tmp[tmp$gene_short_name %in% rownames(x)[x>1],]
tmp = tmp[!tmp$gene_short_name %in% rownames(x)[x>1],]

fpkm = fpkm.size[as.vector(tmp$ID),]
rownames(fpkm) = as.vector(tmp$gene_short_name)
colnames(fpkm) = as.vector(anno$celltype)

fpkm = round(fpkm, digits = 2)
write.table(fpkm, "E70_sig.txt", quote=F, sep="\t")

sum(!rownames(fpkm) %in% rownames(exp))
exp_mix = exp[rownames(fpkm), pd$time == "E7"]
write.table(exp_mix, "E70_mix.txt", quote=F, sep="\t")





#### E7.5

time_i = "E7.5"
obj = readRDS(paste0("obj_", time_i, "_exp.rds"))
anno = readRDS(paste0("obj_", time_i, "_anno.rds"))
anno = anno[colnames(obj),]
anno$celltype = unlist(lapply(as.vector(anno$Anno), function(x) strsplit(x,"[:]")[[1]][2]))
Idents(obj) = as.vector(anno$celltype)

count = GetAssayData(object = obj, slot = "counts")
anno = anno[colnames(obj), c("celltype", "group")]
anno$celltype = gsub(" ", "_", anno$celltype)
table(anno$celltype)
sum(rownames(anno) != colnames(count))

mouse_gene_sub = mouse_gene[rownames(mouse_gene) %in% rownames(count),]
count = count[rownames(mouse_gene_sub),]
rownames(count) = as.vector(mouse_gene_sub$ID)

### Generate DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = anno,
                              design = ~ celltype)

rowRanges(dds) <- ebg[rownames(count)]

### obtain normalized FPKM value, and output this dataset for further use
fpkm.size <- fpkm(dds, robust = TRUE)

mouse_gene_sub$rowSums = Matrix::rowSums(fpkm.size)
mouse_gene_sub = mouse_gene_sub[mouse_gene_sub$gene_short_name %in% rownames(exp),]

### exclude duplicated names 
tmp = mouse_gene_sub %>% group_by(gene_short_name) %>% slice_max(order=rowSums, n=1) %>% as.data.frame()
x = table(tmp$gene_short_name)
tmp[tmp$gene_short_name %in% rownames(x)[x>1],]
tmp = tmp[!tmp$gene_short_name %in% rownames(x)[x>1],]

fpkm = fpkm.size[as.vector(tmp$ID),]
rownames(fpkm) = as.vector(tmp$gene_short_name)
colnames(fpkm) = as.vector(anno$celltype)

fpkm = round(fpkm, digits = 2)
write.table(fpkm, "E75_sig.txt", quote=F, sep="\t")

sum(!rownames(fpkm) %in% rownames(exp))
exp_mix = exp[rownames(fpkm), pd$time == "E7.5"]
write.table(exp_mix, "E75_mix.txt", quote=F, sep="\t")
