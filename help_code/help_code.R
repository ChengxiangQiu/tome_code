
###############################################
### Function: doing clustering using Seurat ###
###############################################

doClusterSeurat <- function(obj, nfeatures = 2500, resolution = 1, k.filter = 200, savePath = NULL, correctCC = FALSE, n_dim = 30, min.dist = 0.75){
    
    if(length(table(obj$group))!=1){
        
        obj.list <- SplitObject(object = obj, split.by = "group")
        
        for (i in 1:length(x = obj.list)) {
            obj.list[[i]] <- NormalizeData(object = obj.list[[i]], verbose = FALSE)
            obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]], 
                                                  selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
        }
        
        reference.list <- obj.list[names(table(obj$group))]
        obj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:n_dim, k.filter = k.filter)
        
        obj.integrated <- IntegrateData(anchorset = obj.anchors, dims = 1:n_dim)
        
        # switch to integrated assay. The variable features of this assay are
        # automatically set during IntegrateData
        DefaultAssay(object = obj.integrated) <- "integrated"
        
        obj <- obj.integrated 
        
    } else {
        
        obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
        obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
        
    }
    
    if(correctCC == TRUE){
        obj <- ScaleData(object = obj, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(obj), verbose = FALSE)
    } else {
        obj <- ScaleData(object = obj, verbose = FALSE)
    }
    
    obj <- RunPCA(object = obj, npcs = n_dim, verbose = FALSE)
    obj <- FindNeighbors(object = obj, dims = 1:n_dim)
    obj <- FindClusters(object = obj, resolution = resolution)
    obj <- RunUMAP(object = obj, reduction = "pca", dims = 1:n_dim, min.dist = min.dist)
    obj <- RunTSNE(object = obj, reduction = "pca", dims = 1:30)
    
    if (!is.null(savePath)){
        pdf(paste0(savePath,"tmp.tsne.pdf"), 10, 8)
        print(DimPlot(obj, reduction = "tsne", label = TRUE, pt.size = 1) + NoLegend())
        dev.off()
        
        pdf(paste0(savePath,"tmp.umap.pdf"), 10, 8)
        print(DimPlot(obj, reduction = "umap", label = TRUE, pt.size = 1) + NoLegend())
        dev.off()
    }
    
    obj
}




#####################################################
### Function: finding ancestor node for each node ###
#####################################################

createLineage_Knn <- function(emb, pd, reduction="umap", replication_times=500, removing_cells_ratio=0.2, k_neigh = 5){
    
    print(dim(emb))
    if(!"Anno" %in% names(pd) | !"day" %in% names(pd)) {print("Error: no Anno or day in pd")}
    if(sum(rownames(pd)!=rownames(emb))!=0) {print("Error: rownames are not matched")}
    pd$state = pd$Anno
    
    res = list()
    
    rep_i = 1
    
    while(rep_i < (replication_times+1)){
        
        sampling_index = sample(1:nrow(pd),round(nrow(pd)*(1-removing_cells_ratio)))
        
        emb_sub = emb[sampling_index,]
        pd_sub = pd[sampling_index,]
        
        irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$day)=="pre",]
        irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$day)=="nex",]
        pd_sub1 <- pd_sub[pd_sub$day == "pre",]
        pd_sub2 <- pd_sub[pd_sub$day == "nex",]
        
        pre_state_min = min(table(as.vector(pd_sub1$state)))
        
        if (pre_state_min < k_neigh & pre_state_min >= 3){
            k_neigh = pre_state_min
            print(k_neigh)
        }
        
        if (pre_state_min < 3){
            next
        }
        
        neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index
        
        tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
        for(i in 1:k_neigh){
            tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
        }
        state1 <- names(table(as.vector(pd_sub1$state)))
        state2 <- names(table(as.vector(pd_sub2$state)))
        
        tmp2 <- matrix(NA,length(state2),length(state1))
        for(i in 1:length(state2)){
            x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
            for(j in 1:length(state1)){
                tmp2[i,j] <- sum(x==state1[j])
            }
        }
        tmp2 <- tmp2/apply(tmp2,1,sum)
        tmp2 <- data.frame(tmp2)
        row.names(tmp2) = state2
        names(tmp2) = state1
        
        res[[rep_i]] = tmp2
        
        rep_i = rep_i + 1
        
    }
    
    return(res)
}


