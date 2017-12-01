##### Tool functions (just a few utility functions)

delta<-function(x,t) ifelse(x==t,1,0)
delta_c<-function(x,t) ifelse(x!=t,1,0)

flatten<-function(A) as.vector(t(A)) ## makes it into a row vector
prune<-function(x,thres) ifelse(x>thres,x,0)
soft_threshold<-function(x,thres) ifelse(x>thres,x-thres,ifelse(x<(-1.0*thres),x+thres,0))
hard_threshold<-function(x,thres) ifelse(abs(x)>thres,x,0)

row_normalize<-function(A){
  normA<-sqrt(apply(A^2, 1, sum))
  B=A
  for ( j in 1: nrow(A)){
    B[j,]<-A[j,]/normA[j]  
  }
  return(B)
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
