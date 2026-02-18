## INPUT gene and gene net OUTPUT adj matrix

library(igraph)

rm(list=ls())

setwd('E:\\SVM\\Realdata\\RD1_SVM')
adj <- read.csv('adjmatrix_comp_UNG.csv')
adj <-as.data.frame(adj[,-1])
dim(adj)
# 拉普拉斯矩阵 ------------------------------------------------------------------

# Non-Normalized Laplacian Matrix from adjacency matrix
Non.NormalizedLaplacianMatrix = function(adj){
  diag(adj) <- 0
  deg <- apply(adj,1,sum)
  D = diag(deg)
  L = D - adj             # 最普通的 L 矩阵 
  return(L)
}

L <- Non.NormalizedLaplacianMatrix(adj)

# 归一化的拉普拉斯矩阵 --------------------------------------------------------------

# Normalized Laplacian Matrix from adjacency matrix
laplacianMatrix = function(adj){
  diag(adj) <- 0                   # 邻接矩阵对角元0
  # 度矩阵元素（对角）--邻接矩阵每行元素的绝对值之和 
  deg <- apply(abs(adj),1,sum)     # abs(adj)-矩阵各元素去绝对值、1-表示按行计算，2表示按列、sum-自定义的调用函数
  p <- ncol(adj)
  L <- matrix(0,p,p)               # p*p 的0元素的 Laplaceian 矩阵
  nonzero <- which(deg!=0)         # 哪些行 元素绝对值之和不为0
  for (i in nonzero){
    for (j in nonzero){
      L[i,j] <- -adj[i,j]/sqrt(deg[i]*deg[j])  # i j 不等时（L 为对称阵）
    }
  }
  diag(L) <- 1                                 # 对角线为1
  return(L)
}

L_norm <- laplacianMatrix(adj)
write.csv(L_norm, 'structural_G.csv')


