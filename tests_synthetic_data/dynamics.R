#### functions for generating random graphs and evolution proceeses


library(nettools)
library(igraph)
source("./tools.R")


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
    heatmap(A,Colv = NA, Rowv = NA,main="heatmap of the adjacency matrix ")
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
generate_realistic_adjacency<-function(N,opts=1,args_l=list(),verbose=TRUE,...){
    ##  Description
    ##  -------------
    ##  Function for generating a simple ER graph, represented by is adjacency matrix
    
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  opts            :   what type of random graph (0: ER, 1: PA, 2: Island, 3: Dot Product, 4: SBM). Default=1
    ##  args_l            :   arguments for the graph
    ##  verbose         :   boolean. Should some summaries be printed as the graph is generated (for debugging) and the output
    ##                      graph plotted?
    ##  ...             :   named arguments for the generation of the graph (as per igraph notations)
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  A               :   adjacency matrix of an ER graph
  if (length(args_l)==0){
        args_l<-list(p=0.1,power=0.9,islands.n=3,islands.size=9,islands.pin=0.3,n.inter=3,K=6,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
  }
  input_list <- as.list(substitute(list(...)))
  print(input_list)
  #print(do.call("pow",$pow))
  ### Change argument list according to what is given as argument to the function
  if (hasArg(p) ){ args_l$p=input_list$p}
  if (hasArg(pow) ){args_l$power<-input_list$pow}
  if (hasArg(islands.n)) args_l$islands.n=input_list$islands.n
  if (hasArg(islands.size) )args_l$islands.size=input_list$islands.size
  if (hasArg(islands.pin)) args_l$islands.pin=input_list$islands.pin
  if (hasArg(n.inter)) args_l$n.inter=input_list$n.inter
  if (hasArg(K)) args_l$K=input_list$K
  if (hasArg(block.sizes)) args_l$block.sizes=eval(input_list$block.sizes)
  if (hasArg(pm)) args_l$pm<-eval(input_list$pm)
  name_type_graph=c('ER', 'Power Law','Island', 'Dot Product','SBM')
  print(paste("Type of graph generated: ",name_type_graph[opts+1]))
  if (opts==0){
      g=erdos.renyi.game(N, args_l$p, type = "gnp", directed = FALSE,loops = FALSE)
  }
  else{
      if (opts==1){
        if (verbose==TRUE) print(paste("power graph: p=",args_l$power))
        power=args_l$power
        g=sample_pa(N, power,directed=FALSE)
      }
      else{
        if (opts==2){
          islands.n=args_l$islands.n
          islands.size=args_l$islands.size
          islands.pin=args_l$islands.pin
          n.inter=args_l$n.inter
          g=sample_islands(islands.n, islands.size, islands.pin, n.inter)
          if (verbose==TRUE) print(paste("island graph: islands.n=",islands.n,"islands.size"=islands.size))
        }
        else{
          if (opts==4){
            pm <- args_l$pm
            print(N)
            block.sizes=c(floor(N/3),floor(N/3),N-2*floor(N/3))
            g <- sample_sbm(N, pref.matrix=pm, block.sizes=block.sizes)
            if (verbose==TRUE) print(paste("stochastic block model: block size", args_l$block.sizes))
          }
          else{
            K=args_l$K
            lpvs <- matrix(rnorm(N*K), K, N)
            lpvs <- apply(lpvs, 2, function(x) { return (abs(x)/sqrt(sum(x^2))) })
            g <- sample_dot_product(lpvs)
            if (verbose==TRUE) print(paste("dot product graph"))
          }
        }
      }
  }
  if (verbose==TRUE){ plot(g)}
  return(g)

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


### ------------------------------------------------------------------------------------------------------
###############################      basic graph alteration function        ##############################
### ------------------------------------------------------------------------------------------------------
random_alteration_adjacency<-function(A,prop,p, p_disp=0, p_creation=0,verbose=FALSE){
    ##  Description
    ##  -------------
    ## A is the adjacency matrix that we wish to transform
    ## A is passed as binary matrix
    ## In this setup, we select prop % of the edges of A, and we randomly reassign them with probability p or delete them
    ## with probability p_disp. Some edges are also created with probability p_creation
    ##  INPUT:
    ##  =============================================================
    ##  A               :   initial (binary) adjacency matrix
    ##  prop            :   proportion of the edges that are potentially modified
    ##  p               :   probability that a candidate ''modifiable'' edge is reassigned somewhere else
    ##  p_disp          :   probability that a candidate ''modifiable'' edge is deleted (p+p_disp<1)
    ##  p_creation      :   probability that a "non-existing" edge is added to the new graph
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  A               :   new modified adjacency matrix
  nb_edges_to_be_modified=ceiling(sum(A)/2*prop)
  nb_edges_to_be_created=rbinom(1,ceiling((sum(1-A)-nrow(A))/2),p_creation) ## have to substract the diagonal
  prob=c(1-p-p_disp,p,p_disp)
  if (verbose){ print(paste("initial number of edges", sum(A)/2))}
  N=nrow(A)
  A_temp=A
  A_temp[upper.tri(A)]=0
  edges<-which(A_temp!=0, arr.ind=TRUE)
  non_edges<-which(A_temp==0, arr.ind=TRUE)
  
  #####---------------------------------
  ##### Edge creation
  #####---------------------------------
  if (nb_edges_to_be_created>0){
      index_edges_created=sample(nrow(non_edges),nb_edges_to_be_created,replace= FALSE)
      for ( i in 1:length(index_edges_modified)){
          e=index_edges_created[i]
          src=non_edges[e,1]
          dest=non_edges[e,2]
          A[src,dest]=1
          A[dest,src]=1
      }
  }
  
  #####---------------------------------
  ##### Edge modification
  #####---------------------------------
  index_edges_modified=sample(nrow(edges),nb_edges_to_be_modified,replace = FALSE)
  for ( i in 1:length(index_edges_modified)){
    e=index_edges_modified[i]
    src=edges[e,1]
    dest=edges[e,2]
    flip=rmultinom(1,1,prob)
    if (which(flip>0)>1){
        if (which(flip>0)==2){
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
    
  }
  A<-delta_c(A,0)
  return(A)
}
### ------------------------------------------------------------------------------------------------------


### ------------------------------------------------------------------------------------------------------
################################      alternative graph alteration process      ##########################
### ------------------------------------------------------------------------------------------------------
graph_alteration<-function(g,m,p=0.1,p_disp=0,p_creation=0.01,m_disp=0,m_creation=0, verbose=FALSE){
    ##  Description
    ##  -------------
    ### This is an alternative to the previous function.
    ### This function alters the graph given as input according to the follwoing process:
    ### m edges are selected at random, and plugged elsewhere with probability  p
    ### independently, m_created edges are created with probability p_created
    ### independently, m_disp edges are deleted with probability p_disp. Here this sparsification process occurs after edges  have
    ### been rewired or created.
    ##  INPUT:
    ##  =============================================================
    ##  g               :   initial graph (igraph object)
    ##  m               :   number of edges that are potentially modified
    ##  p               :   probability that a candidate ''modifiable'' edge is reassigned somewhere else
    ##  m_disp          :   number of edges that are potentially deleted
    ##  p_disp          :   probability that a candidate ''removable'' edge is deleted
    ##  m_disp          :   number of edges that are potentially created
    ##  p_creation      :   probability that a candidate ''creatable'' edge is created
    ##  verbose         :   boolean. Should intermediary messages be printed? (useful for debugging)
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  g_prime         :   new modified graph

  if (m<nrow(get.edgelist(g))){
      to_change<-sample(1:nrow(get.edgelist(g)),m)   ## sample m edges to potentially rewire
  }else{
      to_change<-sample(1:nrow(get.edgelist(g)),floor(nrow(get.edgelist(g))/2))  ## if the number of edges to change is in fact superior to the number of existing edges, change only half.
    m=floor(nrow(get.edgelist(g))/2)
  }

  Adj=get.adjacency(g)
  N<-nrow(Adj)
  list_to_change=get.edgelist(g)[to_change,]
  Adj_prime=Adj ### identical copy to the original adjacency matrix, but with 0 where there can potentially be a modification
  for (i in 1:m){
      Adj_prime[list_to_change[i,1],list_to_change[i,2]]=0  ### set the adjacency to 0 in that case
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
        new_dest=list_to_change[i,2]  ## keeps the original destination
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
  if (verbose){
      print(paste("old:", sum(Adj)))
      print(paste("new:", sum(Adj_prime)))
  }
  g_prime=graph_from_adjacency_matrix(Adj_prime, mode ="undirected")
  return(g_prime)
}
### ------------------------------------------------------------------------------------------------------






















