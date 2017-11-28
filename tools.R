##### Tool functions

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