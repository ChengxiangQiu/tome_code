#############################
#### Section1_E85_embryo ####
#############################

#### The scripts used for analyzing E85 embryo data

#######################################
#### Step2_detect_doublet_clusters ####
#######################################

#### Processing using monocle/alpha

options(warn=-1)
suppressMessages(library(dplyr))
suppressMessages(library(reticulate))
suppressMessages(library(monocle))
suppressMessages(library(flexclust))
suppressMessages(library(mcclust))
suppressMessages(library(VGAM))
suppressMessages(library(ggrastr))
library(reticulate)
import("louvain")
print(packageVersion('monocle'))
work_path = ''

dat = readRDS(paste0(work_path, "/dat.rds"))
count = dat[["gene_count"]]
fd = dat[["df_gene"]]
pd = dat[["df_cell"]]

print(sum(rownames(pd) != colnames(count)))
print(sum(rownames(fd) != rownames(count)))

### cells labeled as doublets (by Scrublet) or from doublet-derived subclusters were filtered out 
pd = pd[!(pd$detected_doublets | pd$doublet_cluster),]

### For each cell, we only retain protein-coding genes, lincRNA genes and pseudogenes
fd = fd[(fd$gene_type %in% c('protein_coding', 'pseudogene', 'lincRNA')) & (!fd$chr %in% c('chrX', 'chrY')),]

count = count[rownames(fd), rownames(pd)]

### Genes expressed in less than 10 cells and cells expressing less than 100 genes were further filtered out
min.features = 100
min.cells = 10

# filter genes on the number of cells expressing
if (min.cells > 0) {
    num.cells <- Matrix::rowSums(x = count > 0)
    count <- count[which(x = num.cells >= min.cells), ]
}

# filer cells on the number of genes expressing
if (min.features > 0) {
    nfeatures <- Matrix::colSums(x = count > 0)
    count <- count[, which(x = nfeatures >= min.features)]
}

sum(!colnames(count) %in% rownames(pd))
sum(!rownames(count) %in% rownames(fd))

pd = pd[colnames(count),]
fd = fd[rownames(count),]

pData = new("AnnotatedDataFrame",data=pd) 
fData = new("AnnotatedDataFrame",data=fd)

cds <- newCellDataSet(count, 
                      phenoData = pData,
                      featureData =fData,
                      expressionFamily = negbinomial.size())

DelayedArray:::set_verbose_block_processing(TRUE)
options(DelayedArray.block.size=1000e6)
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

disp_table = dispersionTable(cds)
disp_table = disp_table %>% 
    mutate(excess_disp = (dispersion_empirical - dispersion_fit) / dispersion_fit) %>%
    arrange(plyr::desc(excess_disp))
top_subset_genes = as.character(head(disp_table, 5000)$gene_id)

cds = setOrderingFilter(cds, top_subset_genes)
cds <- preprocessCDS(cds,  method = 'PCA',
                     norm_method = 'log',
                     num_dim = 50,
                     verbose = T)

cds <- reduceDimension(cds, 
                       max_components = 2,
                       reduction_method = 'UMAP',
                       metric="cosine",
                       min_dist = 0.01,
                       n_neighbors = 50,
                       verbose = T)

cds <- clusterCells(cds,
                    method = 'louvain',
                    res = 1e-6,
                    louvain_iter = 1,
                    verbose = T)

cds <- partitionCells(cds)

saveRDS(cds, paste0(work_path, '/monocle_cds.rds'))

jpeg(paste0(work_path, '/louvain_component.jpeg'), width=10, height=8, units = 'in', res=300)
plot_cell_clusters(cds,
                   color_by = 'louvain_component',
                   cell_size = 0.1,
                   show_group_id = T) +
    theme(legend.text=element_text(size=6)) + #set the size of the text
    theme(legend.position="right") #put the color legend on the right
dev.off()

saveRDS(list(count = exprs(cds),
             pd = data.frame(pData(cds)),
             fd = data.frame(fData(cds))),
        paste0(work_path, "/dat_2.rds"))


#### We computed differential expressed genes across cell clusters with the 
#### differentialGeneTest() function in Monocle/3
#### Of note, here I use Monocle/3 rather than Monocle/alpha, because the function
#### top_markers in Monocle/3 is much faster and works well

library(monocle3)
library(dplyr)
work_path = ''

dat = readRDS(paste0(work_path, "/dat_2.rds"))
exp = dat[['count']]
pd = dat[['pd']]
fd = dat[['fd']]
fd$gene_short_name = fd$gene_name
pd$myCluster = paste0('cluster_', pd$Cluster)

myCluster_table = table(pd$myCluster)
small_cluster = names(myCluster_table)[myCluster_table < 2500]

pd_sub_1 = pd %>%
    filter(!myCluster %in% small_cluster) %>%
    group_by(myCluster) %>%
    sample_n(2500) %>%
    as.data.frame()
pd_sub_2 = pd %>%
    filter(myCluster %in% small_cluster) %>%
    as.data.frame()
pd_sub = rbind(pd_sub_1, pd_sub_2)
rownames(pd_sub) = as.vector(pd_sub$sample)
exp_sub = exp[,rownames(pd_sub)]

cds <- new_cell_data_set(exp_sub,
                         cell_metadata = pd_sub,
                         gene_metadata = fd)

cds <- preprocess_cds(cds, num_dim = 50)
cds = reduce_dimension(cds)
cds = cluster_cells(cds, resolution=1e-5)

marker_test_res <- top_markers(cds, 
                               group_cells_by="myCluster", 
                               reference_cells=1000,
                               cores = 8)

markers <- marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(10, pseudo_R2)

saveRDS(markers, file=paste0(work_path, "/marker_top10.rds"))

#### Of note, after finding top 10 genes of each cell cluster
#### we use monocle/alpha (v2.99) again :)

options(warn=-1)
suppressMessages(library(dplyr))
suppressMessages(library(reticulate))
suppressMessages(library(monocle))
suppressMessages(library(flexclust))
suppressMessages(library(mcclust))
suppressMessages(library(VGAM))
suppressMessages(library(ggrastr))
library(reticulate)
import("louvain")
print(packageVersion('monocle'))
work_path = ''

cds = readRDS(paste0(work_path, '/monocle_cds.rds'))

markers = readRDS(paste0(work_path, "/marker_top10.rds"))

cds$myCluster = paste0("cluster_", cds$Cluster)
cluster_list = as.vector(unique(cds$myCluster))
fData(cds)$gene_short_name = fData(cds)$gene_name

for(i in 1:length(cluster_list)){
    
    print(paste0(i,'/',length(cluster_list)))
    
    cds_subset = cds[,pData(cds)$myCluster == cluster_list[i]]
    
    cds_subset = setOrderingFilter(cds_subset, as.vector(markers$gene_id))
    cds_subset <- preprocessCDS(cds_subset,  
                                method = 'PCA',
                                norm_method = 'log',
                                num_dim = 10,
                                verbose = T)
    
    cds_subset <- reduceDimension(cds_subset, 
                                  max_components = 2,
                                  reduction_method = 'UMAP',
                                  metric="cosine",
                                  min_dist = 0.1,
                                  n_neighbors = 50,
                                  verbose = T)
    
    #### res = 1e-04 for most clustering analysis
    cds_subset <- clusterCells(cds_subset,
                        method = 'louvain',
                        res = 1e-4,
                        louvain_iter = 1,
                        verbose = T)
    
    saveRDS(cds_subset, paste0(work_path, '/doublet_cluster_2/', cluster_list[i], '.rds'))
    
    jpeg(paste0(work_path, '/doublet_cluster_2/', cluster_list[i], '.jpeg'), width=10, height=8, units = 'in', res=300)
    print(plot_cell_clusters(cds_subset,
                             color_by = 'Cluster',
                             show_group_id = T) +
              theme(legend.text=element_text(size=10)) + #set the size of the text
              theme(legend.position="right")) #put the color legend on the right
    dev.off()
    
    for(j in 1:length(cluster_list)){
        markers_sub = markers %>% filter(cell_group == cluster_list[j])
        jpeg(paste0(work_path, '/doublet_cluster_2/', cluster_list[i], "_", cluster_list[j], '.jpeg'), width=10, height=8, units = 'in', res=300)
        print(plot_cell_clusters(cds_subset,
                                 markers = as.character(markers_sub$gene_short_name)))
        dev.off()
    }
    
}


#### Subclusters showing low expression of target cell cluster specific markers
#### and enriched expression of non-target cell cluster specific markers were 
#### annotated as doublets derived subclusters and filtered out in visualization 
#### and downstream analysis.
