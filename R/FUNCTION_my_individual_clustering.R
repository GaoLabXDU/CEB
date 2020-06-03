#individuals cluster methods
# SC3¡¢SIMLR¡¢tSNE+Kmeans¡¢kmeans¡¢Spectral Clustering¡¢DBSCAN¡¢hierarchical clustering

my_sc3 <- function(inputTags, datatype, log.trans, sc3_gene_filter, svm_num_cells, SEED){
    
    library("SC3")
    library("SingleCellExperiment")
  
    exp_cell_exprs <- NULL
    sc3OUTPUT <- NULL

    # Data transformation
    if (datatype == "count") {
        ### For count data, it would be normalized by the total cound number and then log2 transformed
        exp_cell_exprs <- SingleCellExperiment(assays = list(counts = inputTags))
        normcounts(exp_cell_exprs) <- t(t(inputTags)/colSums(inputTags))*1000000
        if (log.trans == TRUE){
          logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)
        } else {
          logcounts(exp_cell_exprs) <- normcounts(exp_cell_exprs)
        }
        
    } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
        ### For CPM, FPKM, RPKM or TPM data, it would be log2 transformed
        exp_cell_exprs <- SingleCellExperiment(assays = list(normcounts = inputTags))
        if (log.trans == TRUE){
          logcounts(exp_cell_exprs) <- log2(normcounts(exp_cell_exprs) + 1)
        } else{
          logcounts(exp_cell_exprs) <- normcounts(exp_cell_exprs)
        }
        
    } 

    rowData(exp_cell_exprs)$feature_symbol <- rownames(exp_cell_exprs)
    exp_cell_exprs <- exp_cell_exprs[!duplicated(rowData(exp_cell_exprs)$feature_symbol), ]

    ### Estimating optimal number of clustering
    exp_cell_exprs <- sc3_estimate_k(exp_cell_exprs)
    optimal_K <- metadata(exp_cell_exprs)$sc3$k_estimation

    ### Clustering by SC3 at the optimal K
    if (ncol(inputTags) < svm_num_cells){
        #print(optimal_K)
        exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = sc3_gene_filter, n_cores = 1, rand_seed = SEED)
        
    } else if (ncol(inputTags) >= svm_num_cells){
        ### Runing SVM
        exp_cell_exprs <- sc3(exp_cell_exprs, ks = optimal_K, biology = FALSE, gene_filter = sc3_gene_filter,
                            svm_max = svm_num_cells, svm_num_cells = svm_num_cells, n_cores = 1, rand_seed = SEED)
        exp_cell_exprs <- sc3_run_svm(exp_cell_exprs, ks = optimal_K)
    }

    ### Exporting SC3 results
    p_Data <- colData(exp_cell_exprs)
    col_name <- paste("sc3_", optimal_K, "_clusters", sep = '')
    sc3OUTPUT <- p_Data[, grep(col_name, colnames(p_Data))]
    return(sc3OUTPUT)
}


my_tSNE_kmeans <- function(inputTags, datatype, log.trans, dimensions, perplexity, k.min, k.max, var_genes, SEED){
    
    library(Rtsne)
    library(ADPclust)
    
    input_lcpm <- NULL
    tsne_input <- NULL
    tsne_output <- NULL
    tsne_kmeansOUTPUT <- NULL
    adpOUTPUT <- NULL

    ### Data transformation
    inputTags_tsne <- inputTags
    if (datatype == "count") {
      ### If the input data is original count data or CPM, it would be tranformed to CPM
      if (log.trans == TRUE){
        tsne_input <- log2(t(t(inputTags_tsne)/colSums(inputTags_tsne))*1000000+1)
      } else (
        tsne_input <- t(t(inputTags_tsne)/colSums(inputTags_tsne))*1000000
      )
      
    } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
      ### If the input data is FPKM or RPKM, we use the transformed TPM data generated before as the input
      if (log.trans == TRUE){
        tsne_input <- log2(inputTags_tsne + 1)
      } else {
        tsne_input <- inputTags_tsne
      }
    }

    if (is.null(var_genes)){
        set.seed(SEED)
        tsne_output <- Rtsne(t(tsne_input), dims = dimensions, perplexity = perplexity, check_duplicates = FALSE)
    } else{
        se_genes = rep(NA, nrow(tsne_input))
        for (i in 1:nrow(tsne_input)){
            se_genes[i] = sqrt(var(tsne_input[i,])/length(tsne_input[i,]))
        }
        decreasing_rank = order(se_genes, decreasing = TRUE)

        set.seed(SEED)

        tsne_output <- Rtsne(t(tsne_input[decreasing_rank[1:var_genes],]), dims = dimensions, perplexity = perplexity)
    }

    ### Determining the optimal cluster number (k) and centroid by ADPclust
    adpOUTPUT <- adpclust(tsne_output$Y, htype = "amise",centroids="auto", nclust = k.min:k.max)

    ### Clustering the cells by kmeans
    tsne_kmeansOUTPUT <- kmeans(tsne_output$Y, tsne_output$Y[adpOUTPUT$centers,], adpOUTPUT$nclust)

    return(tsne_kmeansOUTPUT)
}

my_simple_kmeans <- function(inputTags, datatype, log.trans, iter.max){
  if (datatype == "count"){  # transformed to CPM
    if (log.trans == TRUE){
      input_transformed <- log2(t(t(inputTags)/colSums(inputTags))*1000000+1)
    } else (
      input_transformed <- t(t(inputTags)/colSums(inputTags))*1000000
    )
  } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
    if (log.trans == TRUE){
      input_transformed <- log2(inputTags + 1)
    } else (
      input_transformed <- inputTags
    )
  }
  source('FUNCTION_4 methods for K.R')
  k = best_k_ADP(t(inputTags))
  kmeans_output <- kmeans(t(input_transformed), k, iter.max = iter.max, nstart = 1)
  
}

my_spectral_clustering <- function(inputTags, datatype, log.trans, iter.max){
  library("kernlab")
  if (datatype == "count"){  # transformed to CPM
    if (log.trans == TRUE){
      input_transformed <- log2(t(t(inputTags)/colSums(inputTags))*1000000+1)
    } else (
      input_transformed <- t(t(inputTags)/colSums(inputTags))*1000000
    )
  } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
    if (log.trans == TRUE){
      input_transformed <- log2(inputTags + 1)
    } else (
      input_transformed <- inputTags
    )
  }
  source('FUNCTION_4 methods for K.R')
  k = best_k_ADP(t(inputTags))
  spectral_output <- specc(t(input_transformed), k, iterations = iter.max)
  
}

my_DBSCAN <- function(inputTags, datatype, log.trans, eps, MinPts, SEED){
  library("fpc")
  if (datatype == "count"){  # transformed to CPM
    if (log.trans == TRUE){
      input_transformed <- log2(t(t(inputTags)/colSums(inputTags))*1000000+1)
    } else (
      input_transformed <- t(t(inputTags)/colSums(inputTags))*1000000
    )
  } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
    if (log.trans == TRUE){
      input_transformed <- log2(inputTags + 1)
    } else (
      input_transformed <- inputTags
    )
  }
  dbscan_output <- dbscan(t(input_transformed), eps = eps, MinPts = MinPts, seeds = SEED)
  
}

my_hclust <- function(inputTags, datatype, log.trans, k = 3, distance_form){
  if (datatype == "count"){  # transformed to CPM
    if (log.trans == TRUE){
      input_transformed <- log2(t(t(inputTags)/colSums(inputTags))*1000000+1)
    } else (
      input_transformed <- t(t(inputTags)/colSums(inputTags))*1000000
    )
  } else if (datatype == "CPM" || datatype == "FPKM" || datatype == "RPKM" || datatype == "TPM") {
    if (log.trans == TRUE){
      input_transformed <- log2(inputTags + 1)
    } else (
      input_transformed <- inputTags
    )
  }
  hclust_dist = dist(t(input_transformed), method = distance_form)  # ¼ÆËã¾àÀë
  hclust_output <- hclust(hclust_dist, method = "complete")  # ¾ÛÀà
  hclust_res = cutree(hclust_output, k = k)
  return(hclust_res)
  
}

#---------------------------------------------------------------------
my_individual_clustering <- function(inputTags, datatype = "count", log.trans = T,
                                SC3 = TRUE, gene_filter = FALSE, svm_num_cells = 5000, 
                                tSNE = TRUE, dimensions = 3, perplexity = 30, tsne_min_cells = 200, tsne_min_perplexity = 10, var_genes = NULL,
                                kmeans = TRUE, iter.max = 100,
                                spectral_cluster = TRUE,
                                DBSCAN = TRUE, eps = 270, MinPts = 5,
                                hierarchical_clustering = TRUE, distance_form = "euclidean",
                                SEED = 1){

    cluster_number <- NULL
    cluster_results <- NULL
    inputTags = as.matrix(inputTags)
    
    
    ##### SC3
    if(SC3 == TRUE){
        message("Performing SC3 clustering...")

        sc3OUTPUT <- my_sc3(inputTags = inputTags, datatype = datatype, log.trans = log.trans, sc3_gene_filter = gene_filter,
                            svm_num_cells = svm_num_cells, SEED = SEED)
        cluster_results <- rbind(cluster_results, matrix(c(sc3OUTPUT), nrow = 1, byrow = TRUE))
        cluster_number <- c(cluster_number, max(c(sc3OUTPUT)))
    }
    
    
    ##### tSNE+kmeans
    if(tSNE == TRUE){
        message("Performing tSNE + k-means clustering...")

        ### Dimensionality reduction by Rtsne
        if(length(inputTags[1,]) < tsne_min_cells) {
            perplexity = tsne_min_perplexity
        }

        tsne_kmeansOUTPUT <- my_tSNE_kmeans(inputTags = inputTags, datatype = datatype, log.trans = log.trans, dimensions = dimensions,
                                            perplexity = perplexity, k.min = 2, k.max = max(cluster_number), var_genes = var_genes, SEED = SEED)
        cluster_results <- rbind(cluster_results, matrix(c(tsne_kmeansOUTPUT$cluster), nrow = 1, byrow = TRUE))
        cluster_number <- c(cluster_number, max(c(tsne_kmeansOUTPUT$cluster)))
    }
    
    
    ##### simples kmeans
    if(kmeans == TRUE){
      message("Performing k-means clustering...")
      
      kmeans_output <- my_simple_kmeans(inputTags, datatype = datatype, log.trans = log.trans, iter.max = iter.max)
      cluster_results <- rbind(cluster_results, matrix(c(kmeans_output$cluster), nrow = 1, byrow = TRUE))
      cluster_number <- c(cluster_number, max(c(kmeans_output$cluster)))
    }
    
    
    ##### spectral clustering
    if(spectral_cluster == TRUE){
      message("Performing Spectral clustering...")
      spectral_cluster_output <- my_spectral_clustering(inputTags, datatype = datatype, log.trans = log.trans, iter.max = iter.max)
      cluster_results <- rbind(cluster_results, matrix(c(spectral_cluster_output@.Data), nrow = 1, byrow = TRUE))
      cluster_number <- c(cluster_number, max(c(spectral_cluster_output)))
    }
    
    
    ##### DBSCAN clustering
    if(DBSCAN == TRUE){
      message("Performing DBSCAN clustering...")
      dbscan_output <- my_DBSCAN(inputTags, datatype = datatype, log.trans = log.trans, eps = eps, MinPts = MinPts, SEED = SEED)
      cluster_results <- rbind(cluster_results, matrix(c(dbscan_output$cluster), nrow = 1, byrow = TRUE))
      cluster_number <- c(cluster_number, max(c(dbscan_output$cluster)))
    }
    
    
    ##### hierarchical clustering
    if(hierarchical_clustering == TRUE){
      message("Performing hierarchical clustering...")
      hclust_output <- my_hclust(inputTags, datatype = datatype, log.trans = log.trans, distance_form = distance_form)
      cluster_results <- rbind(cluster_results, matrix(c(hclust_output), nrow = 1, byrow = TRUE))
      cluster_number <- c(cluster_number, max(c(hclust_output)))
    }
    
    individual_clustering_output <- list(cluster_results, cluster_number)
    return(individual_clustering_output)
}
