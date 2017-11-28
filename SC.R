##############################################################################################################################################################
######################################### Investigating the Spectral Clustering Approach ####################################################################
##############################################################################################################################################################




#########################################Spectral Clustering Algorithm ######################################################################################
### Start by creating the affinity matrix (pruned smilarity)
make.affinity <- function(S, n.neighbors=2) {
  N <- length(S[,1])
  
  if (n.neighbors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighbors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  return(A) 
}




Spectral_Clustering<-function(A,nb_clusters,laplacian="unnormalized", plot=F){
  ## Compute  degree matrix D where each diagonal value is
  ## the degree of the respective vertex and all other positions are zero
  k<-nb_clusters
  D <- diag(apply(A, 1, sum)) # sum rows
  
  if (laplacian=="unnormalized"){
    ## Compute unnormalized graph Laplacian (U=D−A)
    U <- D - A
  }
  else{
    #U <- (D %^% (-1/2)) %*% A %*% (D %^% (-1/2)) 
    U <- diag(1,nrow(A))-(D %^% (-1)) %*% A
  }
  evL <- eigen(U, symmetric=TRUE)
  ## Assuming we want k clusters, the next step is to find the k
  ## smallest eigenvectors
  Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
  if (plot==T){
    plot(Z, pch=20)
  }
  
  return(Z)
  
}



"%^%" <- function(M, power)
  with(eigen(M), vectors %*% (values^power * solve(vectors)))

#### Find right permtaion of indices
labels_assignment<-function(output,labels){
  pur=matrix(0,length(unique(output)),length(unique(labels)))
  for (j in unique(output)){
    cache=2*delta(output,j)
    cache_ref=sapply(unique(labels),FUN=function(i){
      return(delta(labels,i))
    })
    cache_ref=t(cache_ref)
    Delta=delta(cache-cache_ref,1)
    pur[j,]=2*apply(Delta,1,sum)/sum(cache)
  }
  labels_out<-apply(pur,1,which.max)
  missing_cluster<-setdiff(unique(labels),unique(labels_out))
  if(length(missing_cluster)>0){
    print ("problem ,missing cluster")
    count=matrix(0,1,length(unique(labels_out)))
    uq=unique(labels_out)
    a=0
    for (i in 1:length(uq)){
      pb=which(labels_out==uq[i])
      count[i]=length(pb)
      if (count[i]>1){
        labels_out[pb]=labels_out[pb]+(a+seq(from=0,to=length(pb)-1,by=1))
        a=a+length(pb)
      }
      
    }
  }
  for (i in labels_out){
    output[which(output==i)]=i
  }
  return(list(purity=pur, labels_clust=labels_out,labels_out=output))
}





