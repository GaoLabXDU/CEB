# FUNCTION_4 methods for K


#---------------------4----PAM(Partitioning Around Medoids) 围绕中心点的分割算法---------------------------
best_k_PAM <- function(dataset){
  library(fpc)
  pamk.best <- pamk(dataset)
  result <- list(pamk.best$nc, pamk.best$pamobject$clustering)
  return(result)
}

#---------------------5----Calinsky criterion---------------------------
best_k_Calinsky <- function(dataset, k.min = 1, k.max = 10, iter = 10){
  library(vegan)
  ca_clust <- cascadeKM(dataset, inf.gr = k.min, sup.gr = k.max, iter = iter)
  #ca_clust$results
  calinski.best <- as.numeric(which.max(ca_clust$results[2,]))
  return(calinski.best)
}

#--------------------6-----Affinity propagation (AP) clustering---------------------------
best_k_AP <- function(dataset){
  library(apcluster)
  ap_clust <- apcluster(negDistMat(r=2), dataset)
  return(length(ap_clust@clusters))
}

#--------------------7-----轮廓系数Average silhouette method---------------------------
best_k_SI <- function(dataset, k.min = 2, k.max = 10){
  library(factoextra)
  #fviz_nbclust(dataset, kmeans, method = "silhouette")
  #return(K_SI)
  library(fpc)
  si_value = NULL
  for (i in k.min:k.max){
	kmeans_output <- kmeans(dataset,i)  #K聚类结果存于kmeans_output变量
	stats=cluster.stats(dist(dataset), kmeans_output$cluster)  #求出聚类评价统计量
	si_value[i]=stats$avg.silwidth  #将结果存入si_value
  }
  return(which.max(si_value)) 
  
}

#-------------------------ADPclust---------------------------
best_k_ADP <- function(dataset, k.min = 2, k.max = 10){
  library(ADPclust)
  adpOUTPUT <- adpclust(dataset, htype = "amise",centroids="auto", nclust = k.min:k.max)
  return(adpOUTPUT$nclust)
}
























