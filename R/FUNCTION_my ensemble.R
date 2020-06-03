# FUNCTION my ensemble by CSPA


#run cluster_ensemble
cluster_ensemble <- function(cluster_results,
                 program.dir = ".",
                 k_min = 2,
                 k_max = NULL,
                 CSPA = TRUE,
                 SEED = 1
){
  #set parameters
  set.seed(SEED)
  # Example for Input data
  #cluster_results = matrix(data = c(1,1,1,2,2,3,3,2,2,2,3,3,1,1,1,1,2,2,3,3,3),nrow = 3, byrow = T)
  
  # Transforming missing data into 0
  for (i in 1:nrow(cluster_results)){
    cluster_results[i,which(is.na(cluster_results[i,]))] <- 0
  }
  
  ###remapping cluster labels where needed ：如果聚类标签不是从1开始的，则都map成从1开始。
  cluster_results = labels_remapping(cluster_results)
  
  if(is.null(k_max)){
    k_max <- max(cluster_results)
  }
  
  #run cspa
  cspa_result <- cspa(object = cluster_results, k = k_max)
  
  #cluster_ensembles <- list()
  
  return(cspa_result)
}

#Example for test
#cluster_results = rbind(c(1,1,1,1,2,2,2), c(1,1,NA,1,2,2,1))
#cluster_results = matrix(data = c(1,1,1,2,2,3,3,2,2,2,3,3,1,1,1,1,2,2,3,3,3),nrow = 3, byrow = T)
#ensemble_results <- cluster_ensemble(cluster_results, program.dir = "D:/R_work/18_SAFE-clustering", CSPA = TRUE, SEED = 123)
#ensemble_results


#-----------------------------------------------------------------------


