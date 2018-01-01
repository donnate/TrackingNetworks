###  Function for testing the effect

source("~/Dropbox/Distances/distances/spanning_trees.R")
source("~/Dropbox/Distances/distances/distances.R")
source("~/Dropbox/Distances/tests_synthetic_data/test_functions.R")
source("~/Dropbox/Distances/tests_synthetic_data/test_ADMM.R")
### generate adjacency matrix
N=30
p=0.4


deletion_process<-function(){
  A<-generate_random_adjacency(N,p, TRUE)
  ### Delete one of the edges
  B=A
  candidate_set<-which(A!=0, arr.ind = TRUE)
  ind=sample(nrow(candidate_set),1)
  i=candidate_set[ind,1]
  j=candidate_set[ind,2]
  B[i,j]=0
  B[j,i]=0
  d_i=sum(B[i,])
  d_j=sum(B[j,])
  B2=B%*%B
  d_i2=sum(B2[i,])
  d_j2=sum(B2[j,])
  
  pred_poly=2/N^2*(1+(1+d_i+d_j)/(N-1)^alpha+1/(N-1)^(2*alpha)*(d_i2+d_j2+d_i*d_j+d_i+d_j+1))
  pred_hamming=2/(N*(N-1))
  poly_distance(A,B)
  
}


plot_difference_Hamming_poly<-function(alpha,N,d_i, d_j, d_i2,d_j2){
  x=seq(from=10, to=1000, by=20)
  diff=matrix(0,length(x),1)
  it=1
  for (N in seq(from=10, to=1000, by=20)){
    pred_poly=2/N^2*(1+(1+d_i+d_j)/(N-1)^alpha+1/(N-1)^(2*alpha)*(d_i2+d_j2+d_i*d_j+d_i+d_j+1))
    pred_hamming=2/(N*(N-1))
    diff[it]=pred_poly-pred_hamming
    it=it+1
  }
  
  
  plot(x,diff)
  
}


