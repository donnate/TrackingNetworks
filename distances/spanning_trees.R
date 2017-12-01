### Compute spanning trees in graphs 

library('Matrix')
library(igraph)


#### Computes spanning tree distance. Makes adjustments if the graph becomes
#### disconnected.

get_number_spanning_trees<-function(A,adjust_disconnection=TRUE){
  ## A is the adjacency matrix of the graph
  ### Note that this works only if the graph is completely connected. Otherwise, it seems even a little meaningless.
  #A<-as.matrix(A)
  N=nrow(A)
  B_ref<-as.matrix(A)
  diag(B_ref)<-apply(A,1,sum)
  if (length(which(diag(B_ref)==0))>0 ) print(paste("disconnected component(s): ",which(diag(B_ref)==0)))
  if (length(which(diag(B_ref)==0))>0 && N>2 && adjust_disconnection==TRUE){
    print(paste("disconnected component(s): ",which(diag(B_ref)==0)))
    index=which(diag(B_ref)==0)
    index_connected<-setdiff(1:N,index)
    A<-A[index_connected,index_connected]
    return(get_number_spanning_trees(A,adjust_disconnection=TRUE))
  }
  else{
    if(N<3){
      nb=ifelse(sum(diag(B_ref))==0,0,1)
    }
    else{
      if (length(which(diag(B_ref)==0))>0 &&  adjust_disconnection==FALSE){
        nb=sum(abs(A))/(N*(N-1)) ## outputs the sparsity instead
      }
      else{
        D=sapply(1:N, FUN=function(i){
          return(sum(A[i,]))
        })
        D=diag(D)
        lambda=eigen(D-A,only.values = T)
        ll=lambda$values[which(lambda$values>10^(-7))]
        nb=-log(N)+sum(log(ll))
        if (nb==0){
          #print("disconnected components")
          graph=graph_from_adjacency_matrix(A, mode = "undirected")
          index_cliques=igraph::components(graph)
          for (l in 1:index_cliques$no){
            selection=which(index_cliques$membership==l)
            nb<-nb+get_number_spanning_trees(A[selection,selection])
          }
          
        }
      }
      
    }
    
    return(nb)
  }
  
}



ST_distance<-function(A,A_new,norm=FALSE){
  if(norm) return(abs(get_number_spanning_trees(A_new)-get_number_spanning_trees(A))/(get_number_spanning_trees(A_new,0)+get_number_spanning_trees(A,0)))
  else return(abs(get_number_spanning_trees(A_new)-get_number_spanning_trees(A)))
  
}



### security check
testit<-function(N=100,pow=2){
  G=erdos.renyi.game(N,0.3)
  #G=sample_pa(N, power = pow,directed=F)
  A=as_adjacency_matrix(G)
  D=sapply(1:N, FUN=function(i){
    return(sum(A[i,]))
  })
  D=diag(D)
  lambda=eigen(D-A,only.values = T)
  nb_true=log(1/N)+sum(log(lambda$values[1:(N-1)]))
  nb_test=get_number_spanning_trees(A)
  print(nb_test)
  return( abs(nb_true-nb_test)<0.0001)
}





