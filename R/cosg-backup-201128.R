
select_top_n<-function(scores,n_top){
  d <- data.frame(
    x   = data.table::copy(scores),
    indice=seq(1,length(scores)))

  data.table::setDT(d)
  data.table::setorder(d,-x)
  n_top_indice<-d$indice[1:n_top]
  return(n_top_indice)
}

#' Marker gene identification for cell groups
#'
#' This function finds marker genes for each of the cell groups in a dataset.
#'
#' @param assay Assay to use in marker gene identification
#' @param slot Slot to pull data from
#' @param alpha The penalty factor to penalize gene expression in cells not belonging to the cluster of interest
#' @param n_genes_user Number of top ranked genes returned in the result
#' @return A list containing two dataframes for ranked marker genes' names and scores, respectively
#' @export
cosg<-function(
  object,
  groups='all',
  assay='RNA',
  slot='data',
  alpha=0.1,
  n_genes_user=100
){

  ### Obtain the cellxgene data
  genexcell<-Seurat::GetAssayData(object = object[[assay]], slot = slot)

  if (groups == 'all'){
    group_info <- Seurat::Idents(object = object)
  }else{
    object <- Seurat::subset(x = object, idents = groups)
    group_info <- Seurat::Idents(object = object)
  }

  ### unique groups
  groups_order=sort(unique(group_info))
  n_cluster=length(groups_order)

  if (n_cluster == 1){
    stop('Cannot perform marker gene identification on a single cluster.')}

  n_cell=ncol(genexcell)
  n_gene=nrow(genexcell)
  gene_name=rownames(genexcell)



  cluster_mat=matrix(0,nrow =n_cluster,ncol = n_cell)
  cluster_mat_reverse=matrix(0,nrow =n_cluster,ncol = n_cell)

  order_i=1
  ### Set gene lambda and gene omega
  for (group_i in groups_order){
    #     print(group_i)
    idx_i=group_info==group_i
    #     print(idx_i[1:4])
    cluster_mat[order_i,idx_i]=1
    cluster_mat_reverse[order_i,!idx_i]=1
    order_i=order_i+1
  }


  cluster_mat_sparse=as(cluster_mat, "dgCMatrix")
  cluster_mat_reverse_sparse=as(cluster_mat_reverse, "dgCMatrix")
  ### Calculate the cosine similarity
  cluster_cosine_m=proxyC::simil(cluster_mat_sparse,genexcell)
  cluster_cosine_m_reverse=proxyC::simil(cluster_mat_reverse_sparse,genexcell)


  rank_stats_names=data.frame(matrix(matrix(), n_genes_user, length(groups_order),
                                     dimnames=list(seq(1,n_genes_user), groups_order)),
                              stringsAsFactors=F)
  rank_stats_scores=data.frame(matrix(matrix(), n_genes_user, length(groups_order),
                                      dimnames=list(seq(1,n_genes_user), groups_order)),
                               stringsAsFactors=F)

  order_i=1
  ### Set gene lambda and gene omega
  for (group_i in groups_order){
    #     print(group_i)
    idx_i=group_info==group_i
    ## Compare the most ideal case to the worst case
    scores1=cluster_cosine_m[order_i,]
    scores2=cluster_cosine_m_reverse[order_i,]
    #     scores=as.vector(scores1-alpha*scores2)
    scores=scores1-alpha*scores2
    #     print(scores)
    global_indices = select_top_n(scores, n_genes_user)
    #         print(gene_name[global_indices[1:5]])

    rank_stats_names[,order_i]=gene_name[global_indices]
    rank_stats_scores[,order_i]=scores[global_indices]

    ### save the group names
    #     factor_name<-c(factor_name,group_i)
    order_i=order_i+1
  }

  ###
  ranks_stats=list(
    names=rank_stats_names,
    scores=rank_stats_scores

  )

  ### return
  return(ranks_stats)

}
