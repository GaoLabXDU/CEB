# 拆解的SNFtool::spectralClustering的函数。
# Similarity is the similarity matrix of samples
# k is the numbers of clusters

spectral_clustering <- function(Similarity, k, normal = FALSE){
  d <- rowSums(Similarity);d
  d[d == 0] <- .Machine$double.eps  #把相似性为0的地方改为一个非常小的数（机器最大精度），防止除数为0
  D <- diag(d) # 生成对角线为d的对角矩阵
  L <- D - Similarity
  
  Di <- diag(1/sqrt(d))
  NL <- Di %*% L %*% Di
  
  eig <- eigen(NL)  #计算特征值和特征向量
  res <- sort(abs(eig$values), index.return = TRUE)  # 按特征值从小到大进行排序
  U <- eig$vectors[, res$ix[1:k]]  # 取出前k个最小特征值对应的特征向量
  
  # 如果参数normal = T，则执行正则化
  normalize <- function(x) x/sqrt(sum(x^2))
  if (normal == T) {
    U <- t(apply(U, 1, normalize))  #apply()是对对象应用某个函数
  }
  
  #eigDiscrete <- .discretisation(U)
  #eigDiscrete <- eigDiscrete$discrete
  #labels <- apply(eigDiscrete, 1, which.max)
  
  label <- kmeans(U,k)$cluster
  return(label)
}





