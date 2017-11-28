### Compute spanning trees in graphs 

library('Matrix')
library(igraph)


get_span_tree<-function(A, B,epsilon=10^(-8)){
  N=nrow(N)
  D=apply(A,1, FUN=function(x){
    return(sum(x))
  })
  D<-diag(D)
  lambda=eigen(D-A,only.values = T)
  non_zer0=which(lambda$values>epsilon)
  l=sum(log(lambda$values[non_zer0]))
  D2=apply(B,1, FUN=function(x){
    return(sum(x))
  })
  D2<-diag(D2)
  lambda2=eigen(D2-B,only.values = T)
  non_zer0=which(lambda2$values>epsilon)
  l2=sum(log(lambda2$values[non_zer0]))
  return(abs(l-l2))
}

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
        #B_tilde<-B_ref[1:(N-1),1:(N-1)]
        D=sapply(1:N, FUN=function(i){
          return(sum(A[i,]))
        })
        D=diag(D)
        lambda=eigen(D-A,only.values = T)
        ll=lambda$values[which(lambda$values>10^(-7))]
        nb=-log(N)+sum(log(ll))
        #nb=det(B_tilde)
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

get_nb_ST<-function(A,adjust4disconnection="addition",epsilon=10^(-8)){
  ## A is the adjacency matrix of the graph
  ### Note that this works only if the graph is completely connected. Otherwise, it seems even a little meaningless.
  #A<-as.matrix(A)
  N=nrow(A)
  G=graph_from_adjacency_matrix(A, mode = "undirected")
  t=clusters(G)
  if (adjust4disconnection=="addition"){
    if (t$no>1){
      nb_ST=0
      for (i in 1:t$no){
        index=which(t$membership==i)
        if (length(index)>1){
          A_tilde=A[index,index]
          
          nb_ST=nb_ST+get_nb_ST(A_tilde)
        }
        else{
          if (length(index)==1){
            nb_ST=nb_ST+1
          }
        }
        
      }
      return(nb_ST)
    }
    else{
      D=apply(A,1, sum)
      D<-diag(D)
      lambda=eigen(D-A,only.values = T)
      lambda=sort(lambda$values)
      return(1/length(lambda)*exp(sum(log(lambda[2:length(lambda)]))))
    }
    
  }
  else{
    D=apply(A,1, sum)
    D<-diag(D)
    lambda=eigen(D-A,only.values = T)
    lambda=sort(lambda$values)
    index=which(lambda>epsilon)
    lambda=lambda[index]
    return(1/length(lambda)*exp(sum(log(lambda))))
  }
  
}



ST_distance<-function(A,A_new,norm=FALSE){
  if(norm) return(abs(get_number_spanning_trees(A_new)-get_number_spanning_trees(A))/(get_number_spanning_trees(A_new,0)+get_number_spanning_trees(A,0)))
  else return(abs(get_number_spanning_trees(A_new)-get_number_spanning_trees(A)))
  
  #return(abs(get_number_spanning_trees(A_new)-get_number_spanning_trees(A))/(get_number_spanning_trees(A_new)+get_number_spanning_trees(A)))
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





