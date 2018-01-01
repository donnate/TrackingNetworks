##### Functions for comparing different dustances' ability to pick up change points

dir="~/TrackingNetworkChanges"
setwd(dir)
source("./tests_synthetic_data/test_functions.R")
source("./spanning_trees.R")
source("./distances.R")


##########################################################################################################
##########################################################################################################
############################         Change point detection          #####################################
##########################################################################################################
##########################################################################################################




##########################################################################################################
##########        Change point detection in independent  graphs  (graph type recognition task)    ########
##########################################################################################################











### ------------------------------------------------------------------------------------------------------
##############################      random ER graph comparison with change point       ###################
### ------------------------------------------------------------------------------------------------------
random_comp_tests<-function(N=30,p=c(0.05,0.2),T=21,verbose=TRUE,save_graph_seq=T,name_file_ext=""){
    ##  Description
    ##  -------------
    ##  Function for generating T+1 random ER graphs, and computing distance between consecutive graphs
    ##  (using different metrics: spanning trees, Hamming, IM, etc).
    ##  At time T, the type of ER graphs generated changes (that is, the parameter p changes from p[2]
    ##  to p[3]), and we wish to study the effect on the registered distance. In this case, there is no
    ##  evolution (in the sense that each graph is a perturbed version of the previous: here every graph
    ##  is independent from the previous): instead, we wish to assess whether the distances are able to
    ##  recognize graphs belonging to the same family.
    
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  p               :   probability of edge connection in an ER graph. (2-d vector of float/double<1), the first
    ##                      (resp. the second) entry corresponding to the ER parameter for the first (resp. second) phase
    ##  T               :   number of different graphs that are generated. (int>1)
    ##  verbose         :   boolean. Should intermediary summary messages be printed? (useful for
    ##                      debugging. Default=TRUE
    ##  save_graph_seq  :   boolean.Should the graphs be saved? Default=TRUE
    ##  name_file_ext   :   string (name of the file were the graphs should be saved if save_graph_seq
    ##                      is TRUE. Default=""
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  list            :   list with named arguments dist_random (distance between consecutive  graphs)
    ##                      and data (a matrix where each column corresponds to a flattened graph)
    
  print("further tests: comparison of completely random graphs, with generating regime change")
  A<-vector("list",T)
  A_new<-vector("list",T-1)
  A[[1]]<-generate_random_adjacency(N,p[1], TRUE)
  graph_seq<-matrix(0,T,N^2)
  graph_seq[1,]<-as.vector(t(A[[1]]))
  
  dist_random<-matrix(0, 16,T-1)
  for (t in 1:(T-1)){
    if (verbose==TRUE) print(paste("time: ",t))
    if(t<10){
      A_new[[t]]<-generate_random_adjacency(N,p[1], TRUE)
    }
    else{
      A_new[[t]]<-generate_random_adjacency(N,p[2], TRUE)
    }
    
    stree_new=get_number_spanning_trees2(A_new[[t]])
    stree=get_number_spanning_trees2(A[[t]])
    ### surface Plot for Z
    Z=A_new[[t]]-A[[t]]
    if (t==1){
      df <- data.frame(x = rep(seq_len(ncol(Z)), each = nrow(Z)),
                       y = rep(seq_len(nrow(Z)), times = ncol(Z)),
                       z = c(Z))
    }
    else{
      df <- df+data.frame(x = rep(seq_len(ncol(Z)), each = nrow(Z)),
                       y = rep(seq_len(nrow(Z)), times = ncol(Z)),
                       z = c(Z))
    }
    
    dist_random[1,t]<-1/N*abs(log(stree_new)-log(stree))
    dist_random[2:4,t]<-netdist(A[[t]], A_new[[t]],d = "HIM")
    dist_random[5,t]<-poly_distance(A[[t]], A_new[[t]],order_max=3,alpha=0.5)
    dist_random[6,t]<-poly_distance(A[[t]], A_new[[t]],order_max=5,alpha=0.7)
    dist_random[7,t]<-poly_distance(A[[t]], A_new[[t]],order_max=3,alpha=0.9)
    dist_random[8,t]=abs(sp_Anew-sp_A)/(sp_A+sp_Anew)
    dist_random[9,t]<-poly_distance2(A[[t]], A_new[[t]],order_max=5,alpha=0.5)
    dist_random[10,t]<-poly_distance2(A[[t]], A_new[[t]],order_max=3,alpha=0.9)
    it_alpha=0
    
    for (alpha in alphas){
      if (it_alpha==0){
        cols=paste("f(l)=exp(-",alpha,"l)",sep="")
      }
      else{
        cols=c(cols,paste("f(l)=exp(-",alpha,"l)",sep=""))
      }
      dist_random[11+it_alpha,t]<-eigen_distance(A[[t]], A_new[[t]],function(x){return(-alpha*x)},p=2)
      it_alpha=it_alpha+1
      
    }
    dist_random[11+it_alpha,t]<-eigen_distance(A[[t]], A_new[[t]],function(x){ifelse(x<2,x,0)})
    it_alpha=it_alpha+1
    cols=c(cols,"f(l)= l. 1(l<2)")
    
    A[[t+1]]<-A_new[[t]]
    graph_seq[t+1,]<-as.vector(t(A[[t+1]]))
  }
  dist_random=data.frame(dist_random,row.names = c("ST","Hamming","IM","HIM","alpha=0.5,order=5","alpha=0.5,order=3","alpha=0.9,order=3","ST norm","poly2 alpha=0.5,order=5","poly2 alpha=0.9,order=3",cols))
  
  plot(1:(T-1),dist_random[2,],ylim=c(0,1), type="l",col=2,xlab="time", ylab="distance", main="Distance between consecutive graphs")
  points(1:(T-1),dist_random[1,], type="l",col=1)
  points(1:(T-1),dist_random[3,], type="l",col=3,xlab="time", ylab="distance between consecutive graphs")
  points(1:(T-1),dist_random[4,], type="l",col=4,xlab="time", ylab="distance between consecutive graphs")
  points(1:(T-1),dist_random[5,]/max(dist_random[5,]), type="l",col=5,xlab="time", ylab="distance between consecutive graphs")
  
  legend(1,1,legend=c( "ST","H","IM","HIM",paste("Poly x",max(dist_random[5,]))),col=1:5,lwd=2,cex=0.7)
  
  
  if (save_graph_seq){
    #dist=data.frame(dist,row.names = c("Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=5"))
    write.table(dist_random,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_matrices_",N,"_",p,"_",name_file_ext,"distances.txt",sep=""),row.names = TRUE,col.names = FALSE)
    write.table(graph_seq,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_matrices_",N,"_",p,"_",name_file_ext,".txt",sep=""),row.names = FALSE,col.names = FALSE)
  }
  return(list(dist_random=dist_random,data=graph_seq))
}
### ------------------------------------------------------------------------------------------------------




### ------------------------------------------------------------------------------------------------------
#######################      random ``realistic'' graph comparison with change point #####################
### ------------------------------------------------------------------------------------------------------
random_comp_tests_realistic<-function(N=30,opts=1,args, args2,T=21,verbose=FALSE,save_graph_seq=T,name_file_ext=""){
    
    ##  Description
    ##  -------------
    ##  Function for generating T+1 random graphs, and computing distance between consecutive graphs
    ##  (using different metrics: spanning trees, Hamming, IM, etc).
    ##  At time T, the type of graphs generated changes (that is, the graph parameters change from args
    ##  to args2), and we wish to study the effect on the registered distance. In this case, there is no
    ##  evolution (in the sense that each graph is a perturbed version of the previous: here every graph
    ##  is independent from the previous): instead, we wish to assess whether the distances are able to
    ##  recognize graphs belonging to the same family.
    
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  opts            :   what type of random graph (1: PA, 2: Island, 3: Dot Product, 4: SBM). Default=1
    ##  args            :   arguments for the graph in the first dynamic period (named list
    ##                      as per igraph documentation)
    ##  args2           :   arguments for the graph in the second dynamic period (named list
    ##                      as per igraph documentation)
    ##  T               :   number of different graphs that are generated. (int>1)
    ##  verbose         :   boolean. Should intermediary summary messages be printed? (useful for
    ##                      debugging. Default=TRUE
    ##  save_graph_seq  :   boolean.Should the graphs be saved? Default=TRUE
    ##  name_file_ext   :   string (name of the file were the graphs should be saved if save_graph_seq
    ##                      is TRUE. Default=""
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  list            :   list with named arguments dist_random (distance between consecutive  graphs)
    ##                      and data (a matrix where each column corresponds to a flattened graph)
    
    
  print("further tests: comparison of completely random graphs, with generating regime change")
  #### Compare against totally random networks
  A<-vector("list",T)
  A_new<-vector("list",T-1)
  Ag<-generate_realistic_adjacency(N,args,opts, verbose)
  plot(Ag, main="initial graph")
  A[[1]]<-as(get.adjacency(Ag),"matrix")
  graph_seq<-matrix(0,T,N^2)
  graph_seq[1,]<-as.vector(t(A[[1]]))
  
  dist_random<-matrix(0, 16,T-1)
  for (t in 1:(T-1)){
    if (verbose==TRUE) print(paste("time: ",t))
    if(t<10){
      #Ag_new<-graph_alteration(Ag,m,p_disp,p_creation,loc=loc)
      Ag_new<-generate_realistic_adjacency(N,args,opts, verbose)
      A_new[[t]]<-get.adjacency(Ag_new)
    }
    else{
      #Ag_new<-graph_alteration(Ag,m2,p_disp2,p_creation,loc=loc)
      Ag_new<-generate_realistic_adjacency(N,args2,opts, verbose)
      A_new[[t]]<-as(get.adjacency(Ag_new),"matrix")
    }
    stree_new=get_nb_ST(A_new[[t]])
    stree=get_nb_ST(A[[t]])
    dist[1,t]<-1/N*abs(log(stree_new)-log(stree))
    dist[2:4,t]<-netdist(A[[t]], A_new[[t]],d = "HIM")
    dist[5,t]<-poly_distance(A[[t]], A_new[[t]],order_max=3,alpha=0.5)
    dist[6,t]<-poly_distance(A[[t]], A_new[[t]],order_max=5,alpha=0.7)
    dist[7,t]<-poly_distance(A[[t]], A_new[[t]],order_max=3,alpha=0.9)
    dist[8,t]=abs(sp_Anew-sp_A)/(sp_A+sp_Anew)
    dist[9,t]<-poly_distance2(A[[t]], A_new[[t]],order_max=5,alpha=0.5)
    dist[10,t]<-poly_distance2(A[[t]], A_new[[t]],order_max=3,alpha=0.9)
    it_alpha=0
    
    for (alpha in alphas){
      if (it_alpha==0){
        cols=paste("f(l)=exp(-",alpha,"l)",sep="")
      }
      else{
        cols=c(cols,paste("f(l)=exp(-",alpha,"l)",sep=""))
      }
      dist[11+it_alpha,t]<-eigen_distance(A[[t]], A_new[[t]],function(x){return(-alpha*x)},p=2)
      it_alpha=it_alpha+1
      
    }
    dist[11+it_alpha,t]<-eigen_distance(A[[t]], A_new[[t]],function(x){ifelse(x<2,x,0)})
    it_alpha=it_alpha+1
    cols=c(cols,"f(l)= l. 1(l<2)")
    
    A[[t+1]]<-A_new[[t]]
    graph_seq[t+1,]<-as.vector(t(A[[t+1]]))
  }
  dist=data.frame(dist,row.names = c("ST","Hamming","IM","HIM","alpha=0.5,order=5","alpha=0.5,order=3","alpha=0.9,order=3","ST norm","poly2 alpha=0.5,order=5","poly2 alpha=0.9,order=3",cols))
  
  #par(mfrow=c(2,2))
  #plot(1:(T-1),dist_random[1,],ylim=c(0,1), type="l",col=2,xlab="time", ylab="distance ST", main="Distance between consecutive graphs-ST")
  plot(1:(T-1),dist_random[2,],ylim=c(0,1), type="l",col=2,xlab="time", ylab="distance", main="Distance between consecutive graphs")
  #points(1:(T-1),dist_random[6,], type="l",col=6,xlab="time", ylab="distance between consecutive graphs")
  #points(1:(T-1),dist_random[7,], type="l",col=7,xlab="time", ylab="distance between consecutive graphs")
  points(1:(T-1),dist_random[1,], type="l",col=1)
  #points(1:(T-1),dist_random[2,], type="l",col=2,xlab="time", ylab="distance between consecutive graphs")
  points(1:(T-1),dist_random[3,], type="l",col=3,xlab="time", ylab="distance between consecutive graphs")
  points(1:(T-1),dist_random[4,], type="l",col=4,xlab="time", ylab="distance between consecutive graphs")
  points(1:(T-1),dist_random[5,]/max(dist_random[5,]), type="l",col=5,xlab="time", ylab="distance between consecutive graphs")
  
  legend(1,1,legend=c( "ST","H","IM","HIM",paste("Poly x",max(dist_random[5,]))),col=1:5,lwd=2,cex=0.7)
  
  
  if (save_graph_seq){
    dist=data.frame(dist,row.names = c("Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=5"))
    write.table(dist,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_comp_real_matrices_",N,"_",m,"_",name_file_ext,"distances.txt",sep=""),row.names = TRUE,col.names = FALSE)
    write.table(graph_seq,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_comp_real_matrices_",N,"_",m,"_",name_file_ext,".txt",sep=""),row.names = FALSE,col.names = FALSE)
  }
  return(list(dist_random=dist_random,data=graph_seq))
}
### ------------------------------------------------------------------------------------------------------















