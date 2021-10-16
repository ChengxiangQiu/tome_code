#############################
#### Section1_E85_embryo ####
#############################

#### The scripts used for analyzing E85 embryo data

######################################
#### Step1_merge_pipeline_outputs ####
######################################

#### to save time, I split all the PCR samples to 10 analyzing batches
#### After running Junyue Cao's sci-RNA-seq3 pipeline
#### https://github.com/JunyueC/sci-RNA-seq3_pipeline
#### merging all the 10 batches to create the data matrix

library(Matrix)
library(dplyr)
library(stringr)
library(ggplot2)
library(gridExtra) 
library(viridis)

data_path = ""
work_path = ""

df_cell_merge = NULL
gene_count_merge = NULL

for(cnt in 1:10){
  print(cnt)
  
  load(paste0(data_path, "/nobackup/output_", cnt, "/report/sci_summary.RData"))

  print(sum(rownames(gene_count) != df_gene$gene_id))
  print(sum(colnames(gene_count) != df_cell$sample))
  
  gene_count_merge = cbind(gene_count_merge, gene_count)
  df_cell_merge = rbind(df_cell_merge, df_cell)  
}

gene_count = gene_count_merge
df_cell = df_cell_merge
rm(gene_count_merge, df_cell_merge)

### read RT_barcode and bbi-sample-sheet to extract sample name (or RT_group) for each individual cell
RT_barcode = read.table(paste0(data_path, "/RT_barcode_file.txt"), as.is=T)
names(RT_barcode) = c("RT_barcode_well", "RT_barcode_sequence")
bbi_sample = read.csv(paste0(data_path, "/bbi-sci-samplesheet.csv"), header=T, as.is=T)
names(bbi_sample) = c("RT_barcode_well", "RT_group", "genome")
RT_barcode = RT_barcode %>% 
  left_join(bbi_sample[,c("RT_barcode_well", "RT_group")], by = "RT_barcode_well")

df_cell$RT_barcode_sequence = str_sub(unlist(lapply(as.vector(df_cell$sample), function(x) strsplit(x,"[.]")[[1]][2])), -10, -1)
df_cell = df_cell %>%
  left_join(RT_barcode[,c("RT_barcode_sequence", "RT_group")], by = "RT_barcode_sequence")

df_cell$UMI_count = Matrix::colSums(gene_count)
gene_count_copy = gene_count
gene_count_copy@x[gene_count_copy@x > 0] = 1
df_cell$gene_count = Matrix::colSums(gene_count_copy)

mouse_gene = read.table(paste0("/help_code/mouse.v12.geneID.txt"), sep="\t", as.is=T, header=T)
rownames(mouse_gene) = as.vector(mouse_gene$gene_ID)
mouse_gene = mouse_gene[as.vector(df_gene$gene_id),]
df_gene$chr = as.vector(mouse_gene$V1)
rownames(df_gene) = unlist(lapply(rownames(df_gene), function(x) strsplit(x,"[.]")[[1]][1]))
df_gene$gene_id = rownames(df_gene)
gene_keep = df_gene$chr %in% paste0("chr", c(1:19, "M", "X", "Y"))

keep = df_cell$UMI_count >= 200 & 
  df_cell$gene_count >= 100 & 
  df_cell$unmatched_rate < 0.4 & 
  !is.na(df_cell$RT_group)
df_cell = df_cell[keep,]
df_gene = df_gene[gene_keep,]
gene_count = gene_count[gene_keep, keep]
rownames(gene_count) = as.vector(df_gene$gene_id)
colnames(gene_count) = as.vector(df_cell$sample)
rownames(df_cell) = as.vector(df_cell$sample)

write.csv(df_cell, paste0(work_path, "/df_cell.csv"))
write.csv(df_gene, paste0(work_path, "/df_gene.csv"))
writeMM(t(gene_count), paste0(work_path, "/gene_count.mtx"))

#### run scrublet using python to detect doublets 
#### python /help_code/run_scrublet.py           

doublet_scores_observed_cells = read.csv(paste0(work_path, "/doublet_scores_observed_cells.csv"), header=F)
doublet_scores_simulated_doublets = read.csv(paste0(work_path, "/doublet_scores_simulated_doublets.csv"), header=F)

df_cell$doublet_score = as.vector(doublet_scores_observed_cells$V1)
df_cell$detected_doublets = df_cell$doublet_score > 0.2

#### checking if sub-clusters include over 15% doublet cells 

global = read.csv(paste0(work_path, "/doublet_cluster/global.csv"), header=T)
main_cluster_list = sort(as.vector(unique(global$louvain)))

res = NULL

for(i in 1:length(main_cluster_list)){
  print(paste0(i, "/", length(main_cluster_list)))
  dat = read.csv(paste0(work_path, "/doublet_cluster/adata.obs.louvain_", (i-1), ".csv"), header=T)
  dat$louvain = as.vector(paste0("cluster_", dat$louvain))
  dat = dat %>%
    left_join(df_cell[,c("sample", "detected_doublets", "doublet_score")], by = "sample")
  
  tmp2 = dat %>%
    group_by(louvain) %>%
    tally() %>%
    dplyr::rename(n_sum = n)
  
  tmp1 = dat %>%
    filter(detected_doublets == "TRUE") %>%
    group_by(louvain) %>%
    tally() %>%
    left_join(tmp2, by = "louvain") %>%
    mutate(frac = n/n_sum) %>%
    filter(frac > 0.15)
  
  dat$doublet_cluster = dat$louvain %in% as.vector(tmp1$louvain) 
  
  p1 = ggplot(dat, aes(umap_1, umap_2, color = louvain)) + geom_point() + theme(legend.position="none") 
  p2 = ggplot(dat, aes(umap_1, umap_2, color = doublet_cluster)) + geom_point()
  p3 = ggplot(dat, aes(umap_1, umap_2, color = detected_doublets)) + geom_point()
  p4 = ggplot(dat, aes(umap_1, umap_2, color = doublet_score)) + geom_point() + scale_color_viridis(option = "plasma")
  p5 = ggplot(dat, aes(umap_1, umap_2, color = UMI_count)) + geom_point() + scale_color_viridis(option = "plasma")
  p6 = ggplot(dat, aes(umap_1, umap_2, color = gene_count)) + geom_point() + scale_color_viridis(option = "plasma")
  
  pdf(paste0(work_path, "/doublet_cluster/adata.obs.louvain_", (i-1), ".pdf"), 12, 18)
  grid.arrange(p1, p2, p3, p4, p5, p6, nrow=3, ncol=2) 
  dev.off()
  
  dat[,!colnames(dat) %in% c("umap_1", "umap_2")]
  dat$main_louvain = (i-1)
  
  res = rbind(res, dat)
}

rownames(res) = as.vector(res$sample)
res = res[rownames(df_cell),]
df_cell$doublet_cluster = res$doublet_cluster

saveRDS(list(df_cell = df_cell,
             df_gene = df_gene,
             gene_count = gene_count),
        paste0(work_path, "/dat.rds"))


