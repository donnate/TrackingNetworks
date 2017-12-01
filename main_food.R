

path2dir='~/Dropbox/TrackingNetworkChanges/'    ##### Needs to be changed accordingly
library(igraph)
setwd(path2dir)
source("./tests_synthetic_data/test_functions.R")
source("./distances/spanning_trees.R")
source("./distances/distances.R")

### Load the data from the files:
nb_graphs=49
filename='./data/recipes/Graphs/keys_graphs.csv'
names_graphs<-data.frame(read.csv(filename,header=FALSE))
colnames(names_graphs)<-c("name","id")
names_graphs$id<-as.numeric(names_graphs$id)
names_graphs<-names_graphs[sort(names_graphs$id,decreasing = FALSE,index.return=TRUE)$ix,]


graph_i<-list()
for (i in 1:49){
  filename=paste('./heat_distances/tests_real/Food/graphs/graph',toString(i-1),'.csv',sep='')
  graph_i[[i]]=read.csv(filename,header=TRUE)
  row.names(graph_i[[i]])=graph_i[[i]][,1]
  graph_i[[i]]=graph_i[[i]][,2:ncol(graph_i[[i]])]
}



DistancesJaccard=matrix(0,49,49)
DistancesHamming=matrix(0,49,49)
Distances=matrix(0,49,49)
DistancesIM=matrix(0,49,49)
DistancesHIM=matrix(0,49,49)
DistancesPoly=matrix(0,49,49)
DistancesPoly2=matrix(0,49,49)
DistancesEigen=matrix(0,49,49)


for (i in 2:49){
  print(i)
  A=graph_i[[i]]
  for (j in 1:(i-1)){
    print(j)
    A_new=graph_i[[j]]
    DistancesJaccard[i,j]=jaccard_based_distance(A,A_new)
    dist<-netdist(A, A_new,d = "HIM")
    DistancesHamming[i,j]=dist[1]
    DistancesIM[i,j]=dist[2]
    DistancesHIM[i,j]=dist[3]
    DistancesPoly[i,j]=poly_distance(A, A_new,order_max=3,alpha=0.5)
    DistancesPoly2[i,j]=poly_distance(A, A_new,order_max=5,alpha=0.9)
    alpha=0.9
    DistancesEigen[i,j]=eigen_distance(A, A_new,function(x){return(-alpha*x)},p=2)
  }
}



DistancesEigen_tot2=DistancesEigen_tot+t(DistancesEigen_tot)
DistancesPoly_tot2=DistancesPoly_tot+t(DistancesPoly_tot)
DistancesPoly2_tot2=DistancesPoly2_tot+t(DistancesPoly2_tot)
DistancesIM_tot2=DistancesIM_tot+t(DistancesIM_tot)
DistancesHIM_tot2=DistancesHIM_tot+t(DistancesHIM_tot)
DistancesHamming_tot2=DistancesHamming_tot+t(DistancesHamming_tot)
DistancesJaccard_tot2=DistancesJaccard_tot+t(DistancesJaccard_tot)



