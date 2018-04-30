#### Implementation of a few distances
require('MASS')
library(sna)


hamming_based_distance<-function(A,B){
    ## we assume that the nodes from A to B are in correspondance (identical, no new nodes or nodes disppearing)
    N=nrow(A)
    return(sum(abs(A-B))/(N*(N-1)))
}



jaccard_based_distance<-function(A1,A2){
 ## we assume that the nodes from A1 to A2 are in correspondance (identical, no new nodes or nodes disppearing) 
  N<-nrow(A1)
  Delta<-abs(A1-A2)
  Union_rough=A1+A2
  Union<-apply(Union_rough,1,FUN=function(x){
    xx<-sapply(x,FUN=function(x){
      return(delta_c(x,0))
    })
  })
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






eigen_distance<-function(A1,A2,f=function(x){return(x)},p=2,type="laplacian"){
  N=nrow(A1)
  if (type=='laplacian'){
    L1=diag(apply(A1,1,sum))-A1
    L2=diag(apply(A2,1,sum))-A2
  }
  else{
    if(type=='norm_laplacian'){
      D1=diag(1.0/sqrt(apply(A1,1,sum)))
      D1[which(is.infinite((D1)))]=0
      L1=matrix(0,N,N)
      diag(L1)=1
      L1=L1-D1%*%A1%*%D1
      
      D2=diag(1.0/sqrt(apply(A2,1,sum)))
      D2[which(is.infinite((D2)))]=0
      L2=matrix(0,N,N)
      diag(L2)=1
      L2=L2-D2%*%A2%*%D2
    }
    else{
      L1=A1
      L2=A2
    }
  }

  decomp1<-eigen(L1, only.values = TRUE,symmetric = TRUE)
  poly1<-f(sort(decomp1$values))
  decomp2<-eigen(L2, only.values = TRUE,symmetric = TRUE)
  poly2<-f(sort(decomp2$values))
  delta=poly1-poly2
  return((sum((poly1-poly2)^p))^(1/p))
}



heat_distance<-function(A1,A2, alpha,p=2){
    f=function(x){return(exp(-alpha*x))}
        
    D1=diag(1/sqrt(apply(A1,1,sum)))
    D1[which(D1==Inf)]=0
    L1=matrix(0,nrow(A1),nrow(A1))
    diag(L1)=1
    L1=L1-D1%*%A1%*%D1
    
    
    D2=diag(1/sqrt(apply(A2,1,sum)))
    D2[which(D2==Inf)]=0
    L2=matrix(0,nrow(A2),nrow(A2))
    diag(L2)=1
    L2=L2-D2%*%A2%*%D2
    
    decomp1<-svd(L1)
    decomp2<-svd(L2)
    poly1<-diag(f(decomp1$d))
    decomp2<-svd(A2)
    poly2<-diag(f(decomp2$d))
    delta=decomp1$u%*%poly1%*%decomp1$v-decomp2$u%*%poly2%*%decomp2$v
    return((sum((delta)^p))^(1/p))
}







