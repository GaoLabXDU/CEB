# ����SNFtool::spectralClustering�ĺ�����
# Similarity is the similarity matrix of samples
# k is the numbers of clusters

spectral_clustering <- function(Similarity, k, normal = FALSE){
  d <- rowSums(Similarity);d
  d[d == 0] <- .Machine$double.eps  #��������Ϊ0�ĵط���Ϊһ���ǳ�С������������󾫶ȣ�����ֹ����Ϊ0
  D <- diag(d) # ���ɶԽ���Ϊd�ĶԽǾ���
  L <- D - Similarity
  
  Di <- diag(1/sqrt(d))
  NL <- Di %*% L %*% Di
  
  eig <- eigen(NL)  #��������ֵ����������
  res <- sort(abs(eig$values), index.return = TRUE)  # ������ֵ��С�����������
  U <- eig$vectors[, res$ix[1:k]]  # ȡ��ǰk����С����ֵ��Ӧ����������
  
  # �������normal = T����ִ������
  normalize <- function(x) x/sqrt(sum(x^2))
  if (normal == T) {
    U <- t(apply(U, 1, normalize))  #apply()�ǶԶ���Ӧ��ĳ������
  }
  
  #eigDiscrete <- .discretisation(U)
  #eigDiscrete <- eigDiscrete$discrete
  #labels <- apply(eigDiscrete, 1, which.max)
  
  label <- kmeans(U,k)$cluster
  return(label)
}




