path2dir='~/Dropbox/TrackingNetworkChanges/'    ##### Needs to be changed accordingly
setwd(path2dir)

source("./distances/spanning_trees.R")
source("./distances/distances.R")
source("./tests_synthetic_data/test_functions.R")
library(nettools)


compare_graphs<-function(graph_seq,distance_type="poly",args=list(order_max=3,alpha=0.9)){
    ##  Description
    ##  -------------
    ##  the input is a sequence of graphs, represented by an array in which, column-wise,
    ##  each first entry indexes the time and the other two represent N x N graph adjacency matrix
    ##
    ##  INPUT:
    ##  =============================================================
    ##  graph_seq       :   (N^2+1) x T marix, where each column is a flattened adjacency matrix
    ##                      (1st entry is the time)
    ##  distance_type   :   what distance is used to compare the graphs (string with choice in  {"poly", "ST","HIM"   }. Default :"poly")
    ##  args            :   additional arguments for the distance that is being used (list with named arguments)
    ##
    ##  OUTPUT
    ##  =============================================================
    ## distances        :   matrix with pairwise distances between graphs
    
  s=dim(graph_seq)
  if (is.null(s)){
    s=c(length(graph_seq))
  }
  print(s)
  distances<-matrix(0, s[1], s[1])
  
  for (i in 1:(s[1]-1)){
    print(paste("Progress:", i/s[1]*100))
    if (distance_type=="poly"){
      if (length(s)>2){
        distances[i,(i+1):s[1]]<-sapply((i+1):s[1], FUN=function(j){
          return(poly_distance(graph_seq[i,,],graph_seq[j,,],order_max=args$order_max,alpha=args$alpha))
        })
      }
      else{
        distances[i,(i+1):s[1]]<-sapply((i+1):s[1], FUN=function(j){
          return(poly_distance(graph_seq[[i]],graph_seq[[j]],order_max=args$order_max,alpha=args$alpha))
        })
      }

    }
    else{
      if (distance_type=="ST"){
        if (length(s)>2){
          distances[i,(i+1):s[1]]<-sapply((i+1):s[1], FUN=function(j){
          return(ST_distance(graph_seq[i,,],graph_seq[j,,]))
        })
        }
        else{
          distances[i,(i+1):s[1]]<-sapply((i+1):s[1], FUN=function(j){
            return(ST_distance(graph_seq[[i]],graph_seq[[j]]))
          })
        }
      }
      else{
        if (length(s)>2){
          distances[i,(i+1):s[1]]<-sapply((i+1):s[1], FUN=function(j){
            dist=netdist(graph_seq[i,,],graph_seq[j,,],d="HIM")
            return(dist[3])
          })
        }
        else{
          distances[i,(i+1):s[1]]<-sapply((i+1):s[1], FUN=function(j){
            dist=netdist(graph_seq[[i]],graph_seq[[j]],d="HIM")
            return(dist[3])
          })
        }
      }
    }
    
  }

  return(distances+t(distances))
}







