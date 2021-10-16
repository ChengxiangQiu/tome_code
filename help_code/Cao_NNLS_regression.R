# here I input two matrix MM1, MM2, and then generate the correlation coefficient matrix for predicting each
# cell type in MM1 with cell types in MM2 with high expressed genes in MM1
correlation_analysis_nnls_spec_gene <- function(MM1, MM2, fold.change = 1.5, top_gene_num = 300, spec_gene_num = 300) {
  cell_names = colnames(MM1)
  gene_median_expr = apply(MM1, 1, median)
  gene_max_expr = apply(MM1, 1, max)
  correlation_matrix = lapply(cell_names, function(x) {
    gene_other_max_expr = apply(MM1[, colnames(MM1) != x],1,  max)
    target_expr = MM1[, colnames(MM1) == x]
    df_tmp = data.frame(target = target_expr, ratio = (target_expr + 1) / (gene_median_expr + 1), gene_id = row.names(MM1))
    gene_markers = (df_tmp %>% filter(ratio > fold.change) %>%  arrange(desc(ratio)) %>% head(top_gene_num))$gene_id
    df_tmp = data.frame(target = target_expr, ratio = (target_expr + 1) / (gene_other_max_expr + 1), gene_id = row.names(MM1))
    top_markers = (df_tmp %>% arrange(desc(ratio)) %>% head(spec_gene_num))$gene_id
    selected_row = (row.names(MM1) %in% c(as.character(top_markers), as.character(gene_markers)))
    MM1_filtered = MM1[selected_row, ]
    MM2_filtered = MM2[selected_row, ]
    nnls_result = nnls::nnls(MM2_filtered, MM1_filtered[, x])
    df_coef = data.frame("cell_name_MM2" = colnames(MM2), "beta" = coef(nnls_result), "cell_name_MM1" = x)
    return(df_coef)
  })
  result = do.call(rbind, correlation_matrix)
  result = result %>% spread(cell_name_MM2, beta)
  result_2 = result %>% select(-cell_name_MM1)
  rownames(result_2) = result$cell_name_MM1
  return(result_2)
}
correlation_analysis_bidirection <- function(MM1, MM2, fold.change = 1.5, top_gene_num = 300, spec_gene_num = 300) {
  com_1 = correlation_analysis_nnls_spec_gene(MM1, MM2, fold.change = fold.change, top_gene_num = top_gene_num, spec_gene_num = spec_gene_num)
  result = as.data.frame(com_1)
  result$source = rownames(result)
  result_1 = result %>% gather(key = target, value = beta_1, 1:ncol(MM2))
  com_2 = correlation_analysis_nnls_spec_gene(MM2, MM1, fold.change = fold.change, top_gene_num = top_gene_num, spec_gene_num = spec_gene_num)
  result = as.data.frame(com_2)
  result$target = rownames(result)
  result_2 = result %>% gather(key = source, value = beta_2, 1:ncol(MM1))
  result = inner_join(result_1, result_2)
  return(result)
}
correlation_analysis_cell_type_selection <- function(tmp_result, marker_num = 1) {
  tmp_result$beta = (2 * (tmp_result$beta_1 + 0.01) * (tmp_result$beta_2 + 0.01))
  tmp = tmp_result %>% select(-beta_1) %>% select(-beta_2) %>% spread(key = source, value = beta)
  whole_matrix = as.matrix(tmp %>% select(-target))
  rownames(whole_matrix) = tmp$target
  # filter the matrix for selecting the top two matched cell type for each cell type
  unique_source = unique(tmp_result$source)
  # select the top two source cell type for each target
  top_target = lapply(unique_source, function(x) {
    return((tmp_result %>% filter(source == x, beta > 0) %>% arrange(desc(beta)) %>% head(marker_num))$target)
  })
  selected_target = unique(as.vector(unlist(top_target)))
  filtered_result = tmp_result %>% filter(target %in% selected_target)
  filtered_matrix = t(whole_matrix[selected_target, ])
  plot_MM = (filtered_matrix)
  plot_MM = (plot_MM) / apply(plot_MM, 1, max)
  return(list(filtered_result, filtered_matrix, plot_MM))
}
#The main function is correlation_analysis_bidirection, it accept two normalized gene expression matrix by cell types (with same gene axis) and report the correlation coefficient matrix by non-negative linear regression
#The input is normalized gene expression matrix. Use monocle cds object as an example:
#  mouse_expr = exprs(cds)
#mouse_expr = t(t(mouse_expr) / colSums(mouse_expr)) * 100000
#mouse_expr = log(mouse_expr + 1)