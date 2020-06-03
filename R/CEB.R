#' CEB: a celltypes identify platform for single-cell multi-omic data
#
#' @param dataset_list a list with multi-omic data
#' @param feature.filter a boolean variable defines whether to perform gene filtering. Default is "TRUE".
#' @param feature.filter.fraction the fraction of feature filter. Default is 0.06.
#' @param datatype defines the type of data, which could be "count", "CPM", "RPKM" and "FPKM". Default is "count".
#' @param log.trans a boolean variable defines whether to perform log transform. Default is "FALSE".
#' @param SC3 a boolean variable that defines whether to cluster cells using SC3 method. Default is "TRUE".
#' @param gene_filter a boolean variable defines whether to perform gene filtering in SC3 method. Default is "FALSE".
#' @param svm_num_cells defines the mimimum number of cells above which SVM will be run. Default is 5000.
#' @param tSNE a boolean variable that defines whether to cluster cells using t-SNE method. Default is "TRUE".
#' @param dimensions sets the number of dimensions wanted to be retained in t-SNE step. Default is 3.
#' @param perplexity sets the perplexity parameter for t-SNE dimension reduction. Default is 30 when number of cells >= tsne_min_cells.
#' @param tsne_min_cells defines the number of cells in input dataset below which
#' \code{tsne_min_perplexity=10} would be employed for t-SNE step. Default is 200.
#' @param tsne_min_perplexity sets the perplexity parameter of t-SNE step for small datasets (number of cells \code{<200}).
#' @param var_genes defines the number of variable genes used by t-SNE analysis, when \code{tSNE = TRUE}. Default is "NULL".
#' @param kmeans a boolean variable that defines whether to cluster cells using k-means method. Default is "TRUE".
#' @param iter.max sets the maximum number of iterations in k-means method. Default is 100.
#' @param spectral_cluster a boolean variable that defines whether to cluster cells using spectral cluster method. Default is "TRUE".
#' @param DBSCAN a boolean variable that defines whether to cluster cells using DBSCAN method. Default is "FALSE".
#' @param eps set the parameter eps of DBSCAN. Default is 270.
#' @param MinPts set the parameter MinPts of DBSCAN. Default is 5.
#' @param hierarchical_clustering a boolean variable that defines whether to cluster cells using hierarchical cluster method. Default is "FALSE".
#' @param distance_form set the distance form used in hierarchical cluster method. Default is "euclidean".
#' @param SEED sets the seed of the random number generator. Setting the seed to a fixed value can
#' produce reproducible clustering results. Default is 1.
#'
#' @return a summary of the ensemble clustering result
#' @examples
#' data(data_RNA)
#' data(data_Methy)
#' dataset_list <- list(data_RNA, data_Methy)
#'
#' data(data_reallabel)
#'
#' CEB_result <- CEB(dataset_list, feature.filter = T, feature.filter.fraction = 0.06,
#' datatype = "count", log.trans = F,
#' SC3 = T, gene_filter = F, svm_num_cells = 5000,
#' tSNE = T, dimensions = 3, perplexity = 30, tsne_min_cells = 200,
#' tsne_min_perplexity = 10, var_genes = NULL,
#' kmeans = T, iter.max = 100,
#' spectral_cluster = T,
#' DBSCAN = F, eps = 270, MinPts = 5,
#' hierarchical_clustering = F, distance_form = "euclidean",
#' SEED = 1)
#'
#' library("mclust")
#'
#' # for data_Clark_multiomic
#' ARI_individual_clustering <- NULL
#' i = 0
#' for(i in 1:2){
#'   for(j in 1:4){
#'     ARI_individual_clustering <- c(ARI_individual_clustering,
#'                                    c(adjustedRandIndex(CEB_result[[1]][[i]][j,], reallabel)))
#'   }
#' }
#' # data_Clark_RNA_filter_log_count + data_Clark_methy_filter_log_count的各4种聚类方法共8个结果ARI值
#' ARI_individual_clustering
#' write.csv(t(ARI_individual_clustering),file="D:\\R_work\\CEB\\ARI_individual_clustering_multiomic.csv")
#'
#'
#' @export CEB_result
CEB <- function(dataset_list, feature.filter = T, feature.filter.fraction = 0.06, datatype = "count", log.trans = F,
                SC3 = T, gene_filter = F, svm_num_cells = 5000,
                tSNE = T, dimensions = 3, perplexity = 30, tsne_min_cells = 200, tsne_min_perplexity = 10, var_genes = NULL,
                kmeans = T, iter.max = 100,
                spectral_cluster = T,
                DBSCAN = F, eps = 270, MinPts = 5,
                hierarchical_clustering = F, distance_form = "euclidean",
                SEED = 1){

  # -----------------------data processing---------------------------
  source('FUNCTION_feature_filter.R')

  dataset_filter_list <- NULL
  if(feature.filter) {
    for(i in dataset_list)
      dataset_filter_list <- c(dataset_filter_list, list(feature_filter(i, fraction = feature.filter.fraction)))
  } else{
    dataset_filter_list <- dataset_list
  }

  inputTags_filter_list <- NULL
  for(i in dataset_filter_list){
    inputTags_filter_list = c(inputTags_filter_list, list(as.matrix(i)))
  }

  # ----------------------run my_individual_clustering----------------------------
  source('FUNCTION_my_individual_clustering.R')

  cluster_results_Clark_list <- NULL
  cluster_results_Clark <- NULL
  cluster_number_Clark <- NULL

  for(i in 1:2){
    word = paste("----------calculating dataset", i, "------------", sep = "")
    message(word)
    cluster_results_filter_log <- NULL
    cluster_results_filter_log <- my_individual_clustering(inputTags_filter_list[[i]], datatype = datatype, log.trans = log.trans,
                                                           SC3 = SC3, gene_filter = gene_filter, svm_num_cells = svm_num_cells,
                                                           tSNE = tSNE, dimensions = dimensions, perplexity = perplexity, tsne_min_cells = tsne_min_cells, tsne_min_perplexity = tsne_min_perplexity, var_genes = var_genes,
                                                           kmeans = kmeans, iter.max = iter.max,
                                                           spectral_cluster = spectral_cluster,
                                                           DBSCAN = DBSCAN, eps = eps, MinPts = MinPts,
                                                           hierarchical_clustering = hierarchical_clustering, distance_form = distance_form,
                                                           SEED = SEED)
    cluster_results_Clark_list <- c(cluster_results_Clark_list, list(cluster_results_filter_log[[1]]))
    cluster_results_Clark <-rbind(cluster_results_Clark, cluster_results_filter_log[[1]])
    cluster_number_Clark <-rbind(cluster_results_Clark, cluster_results_filter_log[[2]])

  }

  # ------------------------run cluster_ensemble--------------------------
  source('FUNCTION_my ensemble functions.R')

  ensemble_results <- NULL
  ensemble_results_list <- NULL
  for(i in cluster_results_Clark_list){
    ensemble_results_each <- cluster_ensemble(i, program.dir = ".", CSPA = TRUE, SEED = 1)
    ensemble_results <- rbind(ensemble_results, ensemble_results_each)
    ensemble_results_list <- rbind(ensemble_results_list, list(ensemble_results_each))
  }

  CEB_output <- list(cluster_results_Clark_list, ensemble_results_list)
  return(CEB_output)



}

