#### functions for generating random graphs and evolution proceeses


library(nettools)
#library(ergm)
library(igraph)



##########################################################################################################
############################         Graph generation         ############################################
##########################################################################################################

############################################      ER graph        ########################################
### ------------------------------------------------------------------------------------------------------
generate_random_adjacency<-function(N,p,sym=TRUE, plot=FALSE){
    ##  Description
    ##  -------------
    ##  Function for generating a simple ER graph, represented by is adjacency matrix
    
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  p               :   probability of edge connection in an ER graph. (float/double<1)
    ##  sym             :   boolean. Is the graph directed (FALSE) or undirected (TRUE)
    ##  plot            :   boolean. Should the final graph be plotted?
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  A               :   adjacency matrix of an ER graph

  ## no self edges in that model
  A<-matrix(0,N,N)
  
  for (i in 1:N){
    if (sym==TRUE){
      A[i,i:N]<-rbinom((N-i+1),1,p)
    }
    else{
      A[i,]<-rbinom(N,1,p)
    }
  }
  diag(A)<-0
  if(sym==TRUE){
    A<-A+t(A)
  }
  A<-delta_c(A,0)
  if (plot==TRUE){
    heatmap.2(A,Rowv=FALSE, Colv=FALSE,dendrogram = "none")
  }
  ## test if the matrix is completely connected. OW sample an edge
  ## This is merely because we want connected components
  dest_sum<-apply(A,1,sum)
  unconnected<-which(dest_sum==0)
  for ( e in unconnected){
    if (e==N){
      edge<-sample(1:(N-1),1)
    }
    else{
      if(e==1){
        edge<-sample(2:N,1)
      }
      else{
        possibilities<-c(1:(e-1), (e+1):N)
        edge<-sample(possibilities,1)
      }
    }
    
    A[e,edge]<-1
    if (sym==TRUE){
      A[edge,e]<-1
    }
  }
  return(A)
}
### -------------------------------------------------------------------------------------------------------


############################         More ''realistic'' random graphs             #########################
### -------------------------------------------------------------------------------------------------------
generate_realistic_adjacency<-function(N,opts=1,verbose=TRUE,...){
    ##  Description
    ##  -------------
    ##  Function for generating a simple ER graph, represented by is adjacency matrix
    
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  opts            :   what type of random graph (1: PA, 2: Island, 3: Dot Product, 4: SBM). Default=1
    ##  verbose         :   boolean. Should some summaries be printed as the graph is generated (for debugging)
    ##  ...             :   named arguments for the generation of the graph (as per igraph notations)
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  A               :   adjacency matrix of an ER graph
  args<-list(power=0.9,islands.n=3,islands.size=9,islands.pin=0.3,n.inter=3,K=6,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
  if (hasArg(power) )args$power=power
  if (hasArg(islands.n)) args$islands.n=islands.n
  if (hasArg(islands.size) )args$islands.size=islands.size
  if (hasArg(islands.pin)) args$islands.pin=islands.pin
  if (hasArg(n.inter)) args$n.inter=n.inter
  if (hasArg(K)) args$K=K
  if (hasArg(block.sizes)) args$block.sizes=block.sizes
  if (hasArg(pm)) args$pm=pm
  print(opts)
  if (opts==1){
    if (verbose==TRUE) print(paste("power graph: p=",args$power))
    power=args$power
    g=sample_pa(N, power,directed=FALSE)
  }
  else{
    if (opts==2){
      islands.n=args$islands.n
      islands.size=args$islands.size
      islands.pin=args$islands.pin
      n.inter=args$n.inter
      g=sample_islands(islands.n, islands.size, islands.pin, n.inter)
      if (verbose==TRUE) print(paste("island graph: islands.n=",islands.n,"islands.size"=islands.size))
    }
    else{
      if (opts==4){
        pm <- args$pm
        print(N)
        block.sizes=c(floor(N/3),floor(N/3),N-2*floor(N/3))
        g <- sample_sbm(N, pref.matrix=pm, block.sizes=block.sizes)
        if (verbose==TRUE) print(paste("stochastic block model: block size", args$block.sizes))
      }
      else{
        K=args$K
        lpvs <- matrix(rnorm(N*K), K, N)
        lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
        g <- sample_dot_product(lpvs)
        if (verbose==TRUE) print(paste("dot product graph"))
      }
    }
  }
  if (verbose==TRUE) plot(g)
  return(g)
  #m2 = ergm(data ~ edges + mutual)
  #myNet<-network.initialize(5,bipartite=3)
}
### -------------------------------------------------------------------------------------------------------




##########################################################################################################
############################         Graph evolution          ############################################
##########################################################################################################



### ------------------------------------------------------------------------------------------------------
sample_new_dest<-function(node,N){
    ##  Description
    ##  -------------
    ##  helper function. Samples a new destinaton for the edge starting at node ''node'', in a graph with N nodes
    
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  node            :   node id of the source of the edge
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  choice          :   new destination for an edge starting at ''node''
  if (node==1){
    candidates=1:(N-1)
  }
  else{ 
    if(node==N){
      candidates=2:N
    }
    else{
      candidates=c(1:(node-1),(node+1):N)
    }
  }
  choice=sample(candidates,1)
  return(choice)
}
### ------------------------------------------------------------------------------------------------------


############################################      ER graph        ########################################
### ------------------------------------------------------------------------------------------------------
random_alteration_adjacency<-function(A,prop,p){
    ##  Description
    ##  -------------
    ## A is the adjacency matrix that we wish to transform
    ## A is passed as binary matrix
    ## In this setup, we select prop % of the edges of A, and we randomly reassign them with probability p
    ##  INPUT:
    ##  =============================================================
    ##  A               :   initial (binary) adjacency matrix
    ##  prop            :   proportion of the edges that are potentially modified
    ##  p               :   probability that a candidate ''modifiable'' edge is reassigned somewhhere else
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  A               :   new modified adjacency matrix
  nb_edges_to_be_modified=ceiling(sum(A)/2*prop)
  #print(paste("initial number of edges", sum(A)/2))
  N=nrow(A)
  A_temp=A
  A_temp[upper.tri(A)]=0
  edges<-which(A_temp!=0, arr.ind=TRUE)
  index_edges_modified=sample(nrow(edges),nb_edges_to_be_modified,replace = FALSE)
  #print(length(index_edges_modified))
  for ( i in 1:length(index_edges_modified)){
    e=index_edges_modified[i]
    src=edges[e,1]
    dest=edges[e,2]
    #print(A[src,dest])
    flip=rbinom(1,1,p)
    #print(sum(flip))
    if (flip>0){
      new_dest=sample(N-sum(A[src,])-1,1)
      sel=c(src,which(A[src,]==1) )
      sel_c=setdiff(1:N,sel)
      new_dest=sel_c[new_dest]
      A[src,new_dest]=1
      A[new_dest,src]=1
    }
    A[src,dest]=0
    A[dest,src]=0
    
  }
  A<-delta_c(A,0)
  return(A)
}
### ------------------------------------------------------------------------------------------------------



################################      graph alteration process      ######################################
### ------------------------------------------------------------------------------------------------------
graph_alteration<-function(g,m,p=0.1,p_disp=0,p_creation=0.01,m_disp=0,m_creation=0){
    ##  Description
    ##  -------------
    ### This function alters the graph given as input according to the follwoing process:
    ### m edges are selected at random, and plugged elsewhere with probability  p
    ### independently, m_disp edges are deleted with probability p_disp
    ### independently, m_created edges are created with probability p_created
    ##  INPUT:
    ##  =============================================================
    ##  g               :   initial graph (igraph object)
    ##  m               :   number of edges that are potentially modified
    ##  p               :   probability that a candidate ''modifiable'' edge is reassigned somewhere else
    ##  m_disp          :   number of edges that are potentially deleted
    ##  p_disp          :   probability that a candidate ''removable'' edge is deleted
    ##  m_disp          :   number of edges that are potentially created
    ##  p_creation      :   probability that a candidate ''creatable'' edge is created
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  g_prime         :   new modified graph

  if (m<nrow(get.edgelist(g))){
    to_change<-sample(1:nrow(get.edgelist(g)),m)
  }else{
    to_change<-sample(1:nrow(get.edgelist(g)),floor(nrow(get.edgelist(g))/2))
    m=floor(nrow(get.edgelist(g))/2)
  }

  Adj=get.adjacency(g)
  N<-nrow(Adj)
  list_to_change=get.edgelist(g)[to_change,]
  #print(list_to_change)
  Adj_prime=Adj ### identical copy to the original adjacency matrix, but with 0 where there can potentially be a modification
  for (i in 1:m){
    Adj_prime[list_to_change[i,1],list_to_change[i,2]]=0
    Adj_prime[list_to_change[i,2],list_to_change[i,1]]=0
  }
  change=rbinom(m,1,p)
  for(i in 1:m){
    src=list_to_change[i,1]
    if(change[i]==1){
      no_good_dest=c(src,which(Adj[src,]==1),list_to_change[i,2] ) ## cannot choose new destination from already existing destinations
      new_dest=sample(setdiff(1:N,no_good_dest),1)
    }
    else{
        new_dest=list_to_change[i,2]
    }
    Adj_prime[src,new_dest]=1
    Adj_prime[new_dest,src]=1
  }
  
  if (m_creation>0){
      C=Adj_prime
      diag(C)<-1
      new_candidate_edges=which(C==0,arr.ind=TRUE)
      nb_created=rbinom(m_creation,1,p_creation)
      created=sample(nrow(new_candidate_edges),sum(nb_created))
      for (ne in created){
          Adj_prime[new_candidate_edges[ne,1],new_candidate_edges[ne,2]]=1
          Adj_prime[new_candidate_edges[ne,2],new_candidate_edges[ne,1]]=1
      }
  }
  if (m_disp>0){
      del_candidate_edges=which(Adj_prime==1,arr.ind=TRUE)
      nb_deleted=rbinom(m_disp,1,p_disp)
      deleted=sample(nrow(del_candidate_edges),sum(nb_deleted))
      for (ne in deleted){
          Adj_prime[del_candidate_edges[ne,1],del_candidate_edges[ne,2]]=0
          Adj_prime[del_candidate_edges[ne,2],del_candidate_edges[ne,1]]=0
      }
  }
  diag(Adj_prime)<-0
  print(paste("old:", sum(Adj)))
  print(paste("new:", sum(Adj_prime)))
  g_prime=graph_from_adjacency_matrix(Adj_prime, mode ="undirected")
  return(g_prime)
}
### ------------------------------------------------------------------------------------------------------


##############################    test_change_point detection in   ER graph  #############################
### ------------------------------------------------------------------------------------------------------
test_change_point_detection<-function(N, T,time_change_point=floor(T/2),p0=0.3,p=0.3,p_new=0.6,prop=0.05,prop_new=0.05, plot=TRUE){
    ##  Description
    ##  -------------
    ##  Compares how different distances detect a dynamical regime change
    ##  INPUT:
    ##  =============================================================
    ##  N                :   number of nodes
    ##  T                :   length of the process horizon
    ##  time_change_point:   time at which the change in dynamics occurs
    ##  p0               :   probability of edge creation in ER graph (initial graph)
    ##  prop             :   proportion of the edges that are potentially modified (first time regime)
    ##  p                :   probability that a candidate ''modifiable'' edge is reassigned somewhhere else (2nd time regime)
    ##  prop_new         :   proportion of the edges that are potentially modified (first time regime)
    ##  p_new            :   probability that a candidate ''modifiable'' edge is reassigned somewhhere else (2nd time regime)
    ##
    ##  OUTPUT
    ##  =============================================================
    ## distance         :    matrix of pairwise distances between consecutive graphs
  
  A<-generate_random_adjacency(N,p,sym=TRUE, plot=FALSE)
  distances<-matrix(0,4,(T-1))
  for ( t in 1:time_change_point){
    A_new<-random_alteration_adjacency(A,prop,p)
    distances[1,t]<-abs(get_number_spanning_trees2(A_new)-get_number_spanning_trees2(A))/(get_number_spanning_trees2(A_new)+get_number_spanning_trees2(A))
    distances[5:7,t]<-netdist(A, A_new,d = "HIM")
    #print(distances[,t])
    A<-A_new
  }
  for ( t in (time_change_point+1):(T-1)){
    A_new<-random_alteration_adjacency(A,prop_new,p_new)
    distances[1,t]<-get_number_spanning_trees2(A, A_new)
    distances[5:7,t]<-netdist(A, A_new,d = "HIM")
    A<-A_new
  }

  if (plot==TRUE){
    pic_name=paste("plot_distance_evolution.jpg")
    jpeg(pic_name)
    plot(1:(T-1),distances[5,],ylim=c(0,1), type="l",col=5,xlab="time", ylab="distance between consecutive graphs", main="Evolution of the distance between consecutive graphs")
    points(1:(T-1),distances[6,], type="l",col=6,xlab="time", ylab="distance between consecutive graphs", main="Evolution of the distance between consecutive graphs")
    points(1:(T-1),distances[7,], type="l",col=7,xlab="time", ylab="distance between consecutive graphs", main="Evolution of the distance between consecutive graphs")
    points(1:(T-1),distances[1,], type="l",col=1,xlab="time", ylab="distance between consecutive graphs", main="Evolution of the distance between consecutive graphs")
    points(1:(T-1),distances[2,], type="l",col=2,xlab="time", ylab="distance between consecutive graphs", main="Evolution of the distance between consecutive graphs")
    points(1:(T-1),distances[3,], type="l",col=3,xlab="time", ylab="distance between consecutive graphs", main="Evolution of the distance between consecutive graphs")
    points(1:(T-1),distances[4,], type="l",col=4,xlab="time", ylab="distance between consecutive graphs", main="Evolution of the distance between consecutive graphs")
    legend(1,1,c( "Ratio difference in ST","Spanning tree of change matrix", "Ratio ST change matrix/ ST original graph", "Sparsity",  "H","IM","HIM"),col=1:7,lwd=2,cex=0.4)
    dev.off()
  }
  return(distances)
  
}
### ------------------------------------------------------------------------------------------------------






############################################      distributions        ########################################
### ------------------------------------------------------------------------------------------------------
get_distribution<-function(N,p_v,prop_v,B){
    ##  Description
    ##  -------------
    ##  Compares the dstributions of consecutie distances between modified graphs for spanning tree distances
    ##  INPUT:
    ##  =============================================================
    ##  N                :   number of nodes
    ##  prop_v           :   proportion of the edges that are potentially modified (vector)
    ##  p_v              :   propbability that t amodifiable edge is modified  (vector)
    ##  B                :   number of bootstrap samples
    ##
    ##  OUTPUT
    ##  =============================================================
    ## test_dist        :    matrix of pairwise distances between consecutive graphs
  test_dist<-array(0, dim=c(length(p_v), length(prop_v), B,4))
  for (i in 1:length(p_v)){
    for (j in 1:length(prop_v)){
      test_dist[i,j,,]<-sapply(1:B,FUN=function(x){
        A<-generate_random_adjacency(N,0.4,sym=TRUE, plot=FALSE)
        A_new<-random_alteration_adjacency(A,prop_v[j],p_v[i])
        stree_new=get_number_spanning_trees2(A_new)
        stree=get_number_spanning_trees2(A)
        dist=matrix(0,4,1)
        dist[1]<-abs(stree_new-stree)/(stree_new+stree)
        dist[2:4]<-netdist(A, A_new,d = "HIM")
        return(dist)
      })
    }
  }
  return(test_dist)
}
### ------------------------------------------------------------------------------------------------------





############################################      test functions        ########################################
### ------------------------------------------------------------------------------------------------------
test_get_dist<-function(){
    ##  Description
    ##  -------------
    ##  Compares hhistogramns of distances 
  N=30
  p_v=c(0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7)
  prop_v=c(0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  B=1000
  get_res=get_distribution(N,p_v,prop_v,B)
  save(get_res,p_v,prop_v,N,file="~/Dropbox/Food_network/spanning_trees/get_res.RData")
  m=matrix(0,length(p_v),length(prop_v))
  sd=matrix(0,length(p_v),length(prop_v))
  for (i in 1:length(p_v)){
    for (j in 1:length(prop_v)){
      m[i,j]=mean(get_res[i,j,])
      sd[i,j]=sd(get_res[i,j,])
    }
  }
  par(mfrow=c(4,4))
  for (i in 1:4){
    for (j in 1:4){
      hist(get_res[i,j,], main=paste("p= ",p_v[i], "prop=",prop_v[j]),breaks=10)
    }
  }
  return(get_res)
}
### ------------------------------------------------------------------------------------------------------
