
library(igraph)
source("~/Dropbox/Distances/distances/ADMM.R")
source("~/Dropbox/Distances/test_synthtetic_data/test_ADMM.R")
source("~/Dropbox/Distances/code/tools.R")


### Load the data from the files:
nb_graphs=49
filename='~/Dropbox/Distances/data/keys_graphs.csv'
names_graphs<-data.frame(read.csv(filename,header=FALSE))
colnames(names_graphs)<-c("name","id")
names_graphs$id<-as.numeric(names_graphs$id)
names_graphs<-names_graphs[sort(names_graphs$id,decreasing = FALSE,index.return=TRUE)$ix,]

### On va ensuite pouvoir comparer les graphes en terme de proximite tout ca
nodes<-read.csv('~/Dropbox/Distances/data/food_map.csv',header=TRUE)
nodes<-data.frame(nodes)
colnames(nodes)<-c("full_name","name","label")
get_graph_i<-function(graphs,nb_graphs,nodes){
  graph_i<-list()
  for (i in 1:nb_graphs){
    print(i)
    edges=which(graphs[[i]]!=0, arr.ind=TRUE)
    edges_id<-sapply(1:nrow(edges), FUN=function(x){
      return(c(nodes$name[edges[x,1]],nodes$name[edges[x,2]]))
    })
    edges_id<-t(edges_id)
    edges_id<-data.frame(edges_id)
    colnames(edges_id)<-c("from","to")
    graph_i[[i]]<-graph_from_data_frame(edges_id,directed=FALSE)
    
    name=paste("~/Dropbox/Food_network/Graphs/Gephi/graph_for_gephi",i,".csv",sep="")
    write.csv(edges_id,name)
  }
  return(graph_i)
}



library(gplots)
### We first need to get a sense of the density of each graphs
assess_density<-function(graph_seq, adjacency=TRUE,N=NULL){
  ## input: graph seq of adjacency matrices(preferred) or  edges list (N being the total number of edges, no self loop)
  if (adjacency==FALSE && is.null(N)){
    N=max(max(graph_seq[[1]][,1]),max(graph_seq[[1]][,2]))
  }
  else{
    if (adjacency==TRUE && is.null(N)){
      N=nrow(graph_seq[[1]])
    }
  }
  density_seq<-matrix(0,1,length(graph_seq))
  for (i in 1:length(graph_seq)){
    if (adjacency==TRUE){
      density_seq[i]=length(which(graph_seq[[i]]!=0))/(N*(N-1))
    }
    else{
      density_seq[i]=nrow(graph_seq[[i]])/(N*(N-1))
    }
    
  }
  return(density_seq)
}


density_test<-function(graphs){
  densities<-assess_density(graphs)
  t1=rep(19,8)
  t2=rep(17,8)
  t3=rep(15,8)
  t4=rep(4,2)
  t=c(t1,t2,t3,t4)
  plot(1:nb_graphs,densities,pch=t,col=1:nb_graphs,xlim=c(0,35),ylim=c(0,0.018))
  legend(27,0.018,names_graphs$name,col=1:nb_graphs,pch=t,cex=0.4)
}

### Note that this kind of structure might be enriched by looking at more data, since there
### seems to be a high correltaion between the sparsity observed 

#### Regress each graph as a function of the others:





graph_expressivity<-function(graph_seq,method="Lasso"){
  ### Input: a sequence of graphs, each containing the same nb of nodes (identified)
  ## Output: a similarity matrix based on  the expressivity
  nb_graphs<-length(graph_seq)
  N<-ncol(graph_seq[[1]])
  Sim<-matrix(0,nb_graphs,nb_graphs)
  for (i in 1:nb_graphs){
    
    A=graph_seq[[i]]
    index_reg<-setdiff(1:nb_graphs,i)
    Av <- as.vector(t(A)) ##
    #B<-as.vector(graph_seq[[index_reg[1]]])
    B=matrix(0,length(Av),nb_graphs-1 )
    it=0
    for (j in index_reg){
      #print(j)
      it<-it+1
      B[,it]<-as.vector(t(graph_seq[[j]]))
    }
    ### Apply ADMM to get expression of A in terms of the B_k
    if (method=="ADMM_L1"){
     reg<-ADMM(Av,B,rho=1,eta=0.1,Z_alpha=100,Z_penalty=1,Time=60,dir=1,plot=TRUE)
    }
    else{
      if (method=="ADMM_L2"){
        reg<-ADMM(Av,B,rho=1,eta=0.1,Z_alpha=100,Z_penalty=2,Time=60,dir=1,plot=TRUE)
      }
      else{
        reg<-cv.glmnet(B,Av,alpha=1)
      }
    }
   if (method!="Lasso"){
       opt_index=reg$suggested.index
       Sim[i,index_reg]=reg$lambda[,opt_index]
   }
   else{
       best.lambda <- which(reg$lambda==reg$lambda.min)
       Sim[i,index_reg]=reg$glmnet.fit$beta[,best.lambda]
   }
     
   
    print(paste("Finished", i/nb_graphs,"%"))
  }
  return(Sim)
}

row_normalize<-function(X){
  if (is.null(nrow(X))){
    normX=sqrt(sum(X^2))
    X<-X/normX
  }
  else{
    normX<-sqrt(apply(X^2,1,FUN=sum))
    for ( i in 1:nrow(X)){
      X[i,]<-X[i,]/normX[i]
    }
  }
  return(X)
  
}





apply_analysis<-function(graphs){
  test=graph_expressivity(graphs)
  hist(test); ## plot the histograms
  print(c(min(test),mean(test),max(test)))
  ### Normalize the rows, to get percentage of how much is explained....
  Sim<-row_normalize(test)
  Sim2<-prune(Sim,0.01)
  ### Check to see if graph is connected
  D<-apply(Sim, 1,sum)
}

##### Plot topology
library(plotly)
grad_innovation<-function(x,y) ifelse(x==y,0,ifelse(x>y,1,-1))
grad_innovation_weighted<-function(x,y) x-y
plot_innovation<-function(A,B,filter=TRUE){
  ## A and B are graphs with similar nodes
  ### Start by filtering the nodes (just in case), so that the problem remains tractable
  if ( filter==TRUE){
    set_nodes=filter_nodes(A,B)
    A<-A[set_nodes,set_nodes]
    B<-B[set_nodes,set_nodes]
  }
  Z<-grad_innovation(A,B)
  Zw<-as.matrix(A-B)
  plot_ly(z = ~Zw) %>% add_surface()
}

#### get only rtest_Lassoelevant nodes
filter_nodes<-function(A,B){
  deg1<-apply(A,1,sum)
  index1=which(deg1==0)
  deg2<-apply(B,1,sum)
  index2=which(deg2==0)
  absent=intersect(index1,index2)
  set=setdiff(1:nrow(A),absent)
  return(set)
}





#### Try to cluster by origin
#### There are several problems that have to be dealt with:
## For the weighted graph, we actually have to take into account the fact that we might not have as many recipes as we ought to for each graph
### it might create some unbalance in the shift that we see
### Wealso have to find a way to say something interesting about this inferred graph structure
### What can se say????
#### Use DP for this? Random Surfers

### Test implementation distances with these graphs:



test_lasso_distance<-function(graphs){
  ### Convert the graph sequence into an array
  ### Dammit we are going to have to use sparse matrix
  K=length(graphs)
  N=nrow(graphs[[1]])
  data<-Matrix(0,K, N^2, sparse=TRUE)
  data_t<-c()
  for (k in 1:K){
    print(k)
    #t=as(as.matrix(graphs[[k]]),"sparseVector")
    #print("done flattening")
    
    data[k,]<-as(as.matrix(graphs[[k]]),"sparseVector")
    print("done flattening")
    #c=sum(delta_c(data[k,],0))
    #print(paste("non zero",c))
    data[k,]<-data[k,]
    ### normalize them by the Lasso distance
  }
  save(data,file="flatten_data.RData")
  Sim=get_Similarity(data)
  save(Sim,file="Sim_data.RData")
  SimH=array(0, c(nrow(Sim),ncol(Sim),3))
  for (i in 1:nrow(Sim)){
    SimH[i,,]<-sapply(1:nrow(Sim), FUN=function(j){
      return(netdist(graphs[[i]], graphs[[j]],d = "HIM"))
    })
  }
  save(data, Sim,SimH,file="Similarities.RData")
  sparsity<-apply(Sim,1, FUN=function(x){
    return(sum(delta_c(x,0)))})
  heatmap.2(Sim, Rowv=FALSE, Colv=FALSE, trace="none")
  ### threshold to get "significant" edges
#   nb_nbrs=3 ### keep 3 closest neighbors
#   for (i in 1:nrow(Sim)){
#     closest=which
#   }
#   ##
#   for (k in 2:15){
#     pk<-post_processing(Sim,k,true_labels = 1:nb_graphs,go_for_plotting = FALSE)
#     ### Run spectral clustering 
#     SC<-pk$Z
#     km <- pk$clusters
#     if(go_for_plotting==TRUE){
#       data_F<-data.frame(SC,row.names=names_graphs$name)
#       names_proj<-sapply(1:k,FUN=function(x){
#         return(paste("Z",x,sep=""))
#       })
#       colnames(data_F)<-names_proj
#       plot(data_F, col=pk$clusters, xlim=c(-0.5,1),pch=20,main="Projections colored by kmeans label")
#       plot(data_F$Z2~data_F$Z1, col=pk$clusters,xlim=c(-.5,0.7), pch=20,main="Projections colored by kmeans label")
#       with(data_F, text(Z2~Z1, labels = row.names(data_F), pos = 4,cex=0.5))
#     }
#     
#   }
#   ### Study who is the closest neighbor in the graph
#   ### i.e who has the maximum lasso coefficient with a positive sign
#   closest_neighbor<-matrix(0,nrow(Sim),2)
#   for (j in 1:nrow(Sim)){
#     closest_neighbor[j,1]<-which.max(Sim[j,])
#     closest_neighbor[j,2]<-max(Sim[j,])
#   }
  ### Or just take the MINIMUM spanning tree to infer a structure
  gg<-graph_from_adjacency_matrix(Sim,weighted=TRUE,mode="upper")
  plot(mst(gg,weights = -edge_attr(gg, "weight") ), vertex.label=row.names(data_F),vertex.size=5)
  
  
  graph_closest<-data.frame(closest_neighbor,row.names=1:nb_graphs)
  colnames(graph_closest)<-c("to","weight")
  graph_closest["from"]<-1:nb_graphs
  graph_relationship<-graph_from_edgelist(as.matrix(graph_closest[,c("from","to")]),directed=TRUE)
  
  plot( graph_relationship, layout=layout_with_fr, vertex.size=graph_closest$weight, vertex.label=row.names(data_F))
  ### Study sparsity of the solution
  return(list(Sim=Sim))
}

#heatmap.2(dist$Sim1,Rowv=FALSE,Colv=FALSE,dendrogram ="none",trace="none")
#### Problem: it's a huge graph that takes  a long time to be put in the right matrix form.
#### Maybe we would  be more luck with the correlation matrix from the financial series 
#### Compare with the distances:
comp_dist<-function(){
  nb_clst=3
  load("~/Dropbox/Food_network/code/RData/HIM_dist_food.RData")
  par(mfrow=c(1,3))
  heatmap.2(dist$Sim1,Rowv=F,Colv=F, dendrogram="none", trace="none", main="Hamming Distances between food graphs")
  heatmap.2(dist$Sim2,Rowv=F,Colv=F, dendrogram="none", trace="none", main="IM Distances between food graphs")
  heatmap.2(dist$Sim3,Rowv=F,Colv=F, dendrogram="none", trace="none", main="HIM Distances between food graphs")
  
  
  h2<-hist(dist$Sim2, xlim=c(0,0.25),col=rgb(0,1,0,1),breaks = 15,xlab="distance",main="distances between food graphs")
  h3<-hist(dist$Sim3,col=rgb(0,0,1,1), add=T,breaks = 15)
  h1<-hist(dist$Sim1, col=rgb(1,0,0,1),add=T)
  legend(0.15,300,legend=c("H","IM","HIM"),col=c(rgb(1,0,0,0.5),rgb(0,1,0,1),rgb(0,0,1,1)),cex=0.7,pch=19)
  
  ### Try to apply Spectral Clustering
  alpha_1<-quantile(dist$Sim1,.5)
  p1<-post_processing(dist$Sim1,thres=alpha_1,nb_clusters = nb_clst,true_labels = 1:nrow(dist$Sim1))
  data_F<-data.frame(p1$Z,row.names=names_graphs$name)
  names_proj<-sapply(1:nb_clst,FUN=function(x){
    return(paste("Z",x,sep=""))
  })
  colnames(data_F)<-names_proj
  gg<-graph_from_adjacency_matrix(dist$Sim1,weighted=TRUE,mode="upper")
  tree_int=mst(gg,weights = -edge_attr(gg, "weight") )
  plot(mst(gg,weights = -edge_attr(gg, "weight") ), vertex.label=row.names(data_F),vertex.size=5)
  plot(data_F$Z2~data_F$Z1, col=p1$clusters,xlim=c(-.5,0.7), pch=20,main="Projections colored by kmeans label")
  with(data_F, text(Z2~Z1, labels = row.names(data_F), pos = 4,cex=0.5))
}
tree_max=mst(gg,weights = -edge_attr(gg, "weight"))
set.seed(1492)
l <- layout.fruchterman.reingold(tree_max, niter=5000, area=vcount(gg)^4*10)


 
plot(tree_max, layout=l, 
     vertex.label=row.names(data_F),
     edge.arrow.size=0.5, 
     vertex.label.cex=0.75, 
     vertex.label.family="Helvetica",
     vertex.label.font=1,
     vertex.shape="circle", 
     vertex.size=5, 
     edge.width=0.5)


### Count pure edges
gg<-graph_from_adjacency_matrix(test_Lasso,weighted=TRUE,mode="upper")
tree_stock=mst(gg,weights = -edge_attr(gg, "weight") )
stock_adj=as(as_adj(tree_stock,"both"),"matrix")
cache=matrix(FALSE, nrow(stock_adj),ncol(stock_adj))
for (i in 3:(nrow(stock_adj)-2)){
  cache[i,(i-2):(i+2)]=TRUE
}
cache[1,1:3]=TRUE
cache[1,1:4]=TRUE
cache[nrow(stock_adj),(nrow(stock_adj)-2):nrow(stock_adj)]=TRUE
cache[nrow(stock_adj)-1,(nrow(stock_adj)-3):nrow(stock_adj)]=TRUE
stock_adj_cache=matrix(0, nrow(stock_adj),ncol(stock_adj))
stock_adj_cache[cache]=stock_adj[cache]

###permuting the labels:
permute_matrix<-function(M){
  M_permuted=apply(M,1,FUN=function(x){
    return(sample(x, length(x),replace = FALSE))
  })
}
B=1000
test_dist<-sapply(1:B, FUN=function(b){
  samp=permute_matrix(stock_adj)
  return(sum(samp[cache]))
})

hist(test_dist)
