#####################################
#### Section4_spatial_inferrence ####
#####################################

#### Inference of the approximate spatial locations of cell states during mouse gastrulation

#########################
#### Step1_corn_plot ####
#########################

#### Of note, the CIBERSORTx analysis is web-based. 
#### Please email Chengxiang Qiu (cxqiu@uw.edu) if you have any problems to 
#### perform the analysis

dat = readRDS("dat.rds")
exp = dat[['exp']]
pd = dat[['pd']]

time_point = c("E55", "E60", "E65", "E70", "E75")

celltype = NULL
for(i in time_point){
  dat = read.table(paste0(i, "_result.txt"), header=T, row.names=1, as.is=T)
  dat = dat[,!colnames(dat) %in% c("P.value", "Correlation", "RMSE")]
  celltype = c(celltype, colnames(dat))
}
celltype = as.vector(unique(celltype))


df = NULL
for(i in time_point){
  dat = read.table(paste0(i, "_result.txt"), header=T, row.names=1, as.is=T)
  dat = dat[,!colnames(dat) %in% c("P.value", "Correlation", "RMSE")]
  
  celltype_tmp = celltype[!celltype %in% colnames(dat)]
  tmp = data.frame(matrix(0, nrow(dat), length(celltype_tmp)))
  colnames(tmp) = celltype_tmp; rownames(tmp) = rownames(dat)
  
  x = cbind(dat, tmp)
  x = x[,celltype]
  
  df = rbind(df, x)
}


df = t(df)
sum(colnames(df) != as.vector(pd$name))

df = round(df, digits = 6)

write.table(df, "./my_example.txt", quote=F, sep="\t")

#### The corn plot was generated using code provided in the paper
#### https://www.nature.com/articles/s41586-019-1469-8
