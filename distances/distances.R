#### Other distances betwwen graphs that one can try (for temporal evolution of networks)
source("~/Dropbox/Distances/tools.R")
require('MASS')

jaccard_based_distance<-function(A1,A2){
 ## we assume that the nodes from A1 to A2 are in correspondance (identical, no new nodes or nodes disppearing) 
  N<-nrow(A1)
  Delta<-abs(A1-A2)
  #Delta<-sapply(Delta,FUN=function(x){
  #  return(prune(x,0))
  #})
  ## Delta keeps track of all the edges that are in A1 but not in A2
  #Inter=A1-Delta ## the only non zero coefficient are the edges that are present in both A1 and A2
  
  Union_rough=A1+A2
  Union<-apply(Union_rough,1,FUN=function(x){
    xx<-sapply(x,FUN=function(x){
      return(delta_c(x,0))
    })
  })
  #return(sum(Inter)/sum(Union))
  return(sum(Delta)/sum(Union))
  
}



poly_distance<-function(A1,A2,order_max,weights=NULL,alpha=1){
  ## we assume that the nodes from A1 to A2 are in correspondance (identical, no new nodes or nodes disppearing) 
  N<-nrow(A1)
  if (is.null(weights)){
    weights<-sapply(1:order_max, FUN=function(k){
      return(1/(N-1)^(alpha*(k-1)))
    })
  }
 
  decomp1<-svd(A1)
  poly1<-sapply(1:order_max,FUN=function(k){
    weights[k]*decomp1$d^k
  })
  exp1<-decomp1$u%*%diag(apply(poly1,1,sum))%*%t(decomp1$v)
  decomp2<-svd(A2)
  poly2<-sapply(1:order_max,FUN=function(k){
    weights[k]*decomp2$d^k
  })
  exp2<-decomp2$u%*%diag(apply(poly2,1,sum))%*%t(decomp2$v)
  Delta<-exp1-exp2
  return(1/N^2*sum(sapply(as.vector(Delta),FUN=function(x){x^2})))
}

ginv<-function(X, tol = sqrt(.Machine$double.eps))
{
  ## Generalized Inverse of a Matrix
  dnx <- dimnames(X)
  if(is.null(dnx)) dnx <- vector("list", 2)
  s <- svd(X)
  nz <- s$d > tol * s$d[1]
  structure(
    if(any(nz)) s$v[, nz] %*% (t(s$u[, nz])/s$d[nz]) else X,
    dimnames = dnx[2:1])
}

poly_distance2<-function(A1,A2,order_max,weights=NULL,alpha=1){
  ## we assume that the nodes from A1 to A2 are in correspondance (identical, no new nodes or nodes disppearing) 
  N<-nrow(A1)
  if (is.null(weights)){
    weights<-sapply(1:order_max, FUN=function(k){
      return(1/(N-1)^(alpha*(k-1)))
    })
  }
  
  for (k in 1:order_max){
    if (k==1){
      I=weights[k]*diag(nrow(A1))
      D=0.5*(weights[k]*(A1 %^% (k))%*%ginv(A2 %^% (k))+weights[k]*(A2 %^% (k))%*%ginv(A1 %^% (k)))
    }
    else{
      I=I+weights[k]*diag(nrow(A1))
      D=D+0.5*(weights[k]*(A1 %^% (k))%*%ginv(A2 %^% (k))+weights[k]*(A2 %^% (k))%*%ginv(A1 %^% (k)))
    }
  }
  
  norm=sum(sapply(as.vector(D-I),FUN=function(x){x^2}))/N^2
  return(norm)
}

library(sna)

hamming_based_distance<-function(A,B){
  N=nrow(A)
  return(sum(abs(A-B))/(N*(N-1)))
}



eigen_distance<-function(A1,A2,f=function(x){return(x)},p=2){
  decomp1<-svd(A1)
  poly1<-f(sort(decomp1$d))
  decomp2<-svd(A2)
  poly2<-f(sort(decomp2$d))
  delta=poly1-poly2
  return((sum((poly1-poly2)^p))^(1/p))
}








