##### Design other tests to assess the pertinence of the metrics proposed
##### Note that we would have to prove in each case that these are indeed metrics.
dir="/Users/cdonnat/Dropbox/TrackingNetworkChanges"
setwd(dir)
source("./tests_synthetic_data/test_functions.R")
source("./spanning_trees.R")
source("./distances.R")



random_comp_tests<-function(N=30,p=0.05,prop=0.05,T=21,verbose=TRUE,compare_LASSO=FALSE,save_graph_seq=T,name_file_ext=""){
  print("further tests: comparison of completely random graphs, with generating regime change")
  #### Compare against totally random networks
  A<-vector("list",T)
  A_new<-vector("list",T-1)
  A[[1]]<-generate_random_adjacency(N,p, TRUE)
  graph_seq<-matrix(0,T,N^2)
  graph_seq[1,]<-as.vector(t(A[[1]]))
  
  dist_random<-matrix(0, 16,T-1)
  for (t in 1:(T-1)){
    if (verbose==TRUE) print(paste("time: ",t))
    if(t<10){
      A_new[[t]]<-generate_random_adjacency(N,p, TRUE)
    }
    else{
      A_new[[t]]<-generate_random_adjacency(N,0.6, TRUE)
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
  
  if(compare_LASSO==TRUE){
    Sim_Lasso=get_Similarity(graph_seq)
    heatmap.2(Sim_Lasso+t(Sim_Lasso),Rowv=F, Colv=F, dendrogram="none",main="heatmap of the Lasso Similarity")
  }
  else{
    Sim_Lasso=FALSE
  }
  if (save_graph_seq){
    #dist=data.frame(dist,row.names = c("Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=5"))
    write.table(dist_random,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_matrices_",N,"_",p,"_",name_file_ext,"distances.txt",sep=""),row.names = TRUE,col.names = FALSE)
    write.table(graph_seq,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_matrices_",N,"_",p,"_",name_file_ext,".txt",sep=""),row.names = FALSE,col.names = FALSE)
  }
  return(list(dist_random=dist_random,data=graph_seq,Sim_Lasso=Sim_Lasso))
}


random_comp_tests_realistic<-function(N=30,m,m2,p_disp.p_disp2,p_creation=0.01,T=21,loc=FALSE,verbose=FALSE,compare_LASSO=FALSE,save_graph_seq=T,name_file_ext=""){
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
      Ag_new<-generate_realistic_adjacency(N,args,opts, verbose)
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
  
  if(compare_LASSO==TRUE){
    Sim_Lasso=get_Similarity(graph_seq)
    heatmap.2(Sim_Lasso+t(Sim_Lasso),Rowv=F, Colv=F, dendrogram="none",main="heatmap of the Lasso Similarity")
  }
  else{
    Sim_Lasso=FALSE
  }
  if (save_graph_seq){
    dist=data.frame(dist,row.names = c("Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=5"))
    write.table(dist,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_comp_real_matrices_",N,"_",m,"_",name_file_ext,"distances.txt",sep=""),row.names = TRUE,col.names = FALSE)
    write.table(graph_seq,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_comp_real_matrices_",N,"_",m,"_",name_file_ext,".txt",sep=""),row.names = FALSE,col.names = FALSE)
  }
  return(list(dist_random=dist_random,data=graph_seq,Sim_Lasso=Sim_Lasso))
}


random_comp_chnge_point_tests<-function(N,p,prop=0.05,T=21,verbose=TRUE,compare_LASSO=FALSE,save_graph_seq=T,name_file_ext=""){
  print("further tests: comparison for random innovation, with regime change at t=10 and 16:")
  print(paste("At t=10, the probability of a new link increases from ", p , "to 0.6, (prop=5% of edge changes)"))
  print(paste("At t=16, the probability of a new link goes back to ", p , "but 20% of the edges are unplugged and replugged"))
  #### Compare against totally random networks
  A<-vector("list",T)
  A_new<-vector("list",T-1)
  A[[1]]<-generate_random_adjacency(N,p, TRUE)
  plot(graph_from_adjacency_matrix(A[[1]],mode = "undirected"), main="initial graph")
  graph_seq<-matrix(0,T,N^2)
  graph_seq[1,]<-as.vector(t(A[[1]]))
  dist_changepnt<-matrix(0, 16,T-1)
  prop=0.05
  T=21
  for (t in 1:(T-1)){
    if (verbose==TRUE) print(paste("time: ",t))
    if(t<10){
      A_new[[t]]<-random_alteration_adjacency(A[[1]],0.05,p)
    }
    else{
      if(t<16){
        A_new[[t]]<-random_alteration_adjacency(A[[t]],0.05,0.6)
      }
      else{
        A_new[[t]]<-random_alteration_adjacency(A[[t]],0.2,p)
      }
    }
    
    graph_seq[t+1,]<-as.vector(t(A[[t+1]]))
    stree_new=get_nb_ST(A_new[[t]])
    stree=get_nb_ST(A[[t]])
    dist_changepnt[1,t]<-1/N*abs(log(stree_new)-log(stree))
    dist_changepnt[2:4,t]<-netdist(A[[t]], A_new[[t]],d = "HIM")
    dist_changepnt[5,t]<-poly_distance(A[[t]], A_new[[t]],order_max=3,alpha=0.5)
    dist_changepnt[6,t]<-poly_distance(A[[t]], A_new[[t]],order_max=5,alpha=0.7)
    dist_changepnt[7,t]<-poly_distance(A[[t]], A_new[[t]],order_max=3,alpha=0.9)
    dist_changepnt[8,t]=abs(sp_Anew-sp_A)/(sp_A+sp_Anew)
    dist_changepnt[9,t]<-poly_distance2(A[[t]], A_new[[t]],order_max=5,alpha=0.5)
    dist_changepnt[10,t]<-poly_distance2(A[[t]], A_new[[t]],order_max=3,alpha=0.9)
    it_alpha=0
    
    for (alpha in alphas){
      if (it_alpha==0){
        cols=paste("f(l)=exp(-",alpha,"l)",sep="")
      }
      else{
        cols=c(cols,paste("f(l)=exp(-",alpha,"l)",sep=""))
      }
      dist_changepnt[11+it_alpha,t]<-eigen_distance(A[[t]], A_new[[t]],function(x){return(-alpha*x)},p=2)
      it_alpha=it_alpha+1
      
    }
    dist_changepnt[11+it_alpha,t]<-eigen_distance(A[[t]], A_new[[t]],function(x){ifelse(x<2,x,0)})
    it_alpha=it_alpha+1
    cols=c(cols,"f(l)= l. 1(l<2)")
    A[[t+1]]<-A_new[[t]]
  }
  dist_changepnt=data.frame(dist_changepnt,row.names = c("ST","Hamming","IM","HIM","alpha=0.5,order=5","alpha=0.5,order=3","alpha=0.9,order=3","ST norm","poly2 alpha=0.5,order=5","poly2 alpha=0.9,order=3",cols))
  
  #dist_changepnt=data.frame(dist_changepnt,row.names = c("Hamming1","Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=5","alpha=0.5,order=3","alpha=0.9,order=3","ST norm","poly2 alpha=0.5,order=5","poly2 alpha=0.9,order=3",cols))
  
  
  
  plot(1:(T-1),dist_changepnt[2,],ylim=c(0,1), type="l",col=2,xlab="time", ylab="distance between consecutive graphs", main="Evolution of the distance between consecutive graphs (with change point)")
  points(1:(T-1),dist_changepnt[1,], type="l",col=1,xlab="time", ylab="distance between consecutive graphs")
  points(1:(T-1),dist_changepnt[3,], type="l",col=3,xlab="time", ylab="distance between consecutive graphs")
  points(1:(T-1),dist_changepnt[4,], type="l",col=4,xlab="time", ylab="distance between consecutive graphs")
  points(1:(T-1),dist_changepnt[5,]/max(dist_changepnt[5,]), type="l",col=5,xlab="time", ylab="distance between consecutive graphs")
  legend(1,1,c( "ST-based dist",  "H","IM","HIM",paste("Poly x",max(dist_changepnt[5,]))),col=1:5,lwd=2,cex=0.7)
  if(compare_LASSO==TRUE){
    Sim_Lasso=get_Similarity(graph_seq)
    heatmap.2(Sim_Lasso+t(Sim_Lasso),Rowv=F, Colv=F, dendrogram="none",main="heatmap of the Lasso Similarity")
  }
  else{
    Sim_Lasso=FALSE
  }
  if (save_graph_seq){
    #dist=data.frame(dist_changepnt,row.names = c("ST-based dist",  "H","IM","HIM","alpha=0.5,order=5"))
    write.table(dist_changepnt,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_comp_changepnt_",N,"_",m,"_",name_file_ext,"distances.txt",sep=""),row.names = TRUE,col.names = FALSE)
    write.table(graph_seq,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_comp_changepnt_",N,"_",p,"_",name_file_ext,".txt",sep=""),row.names = FALSE,col.names = FALSE)
  }
  return(list(dist_changepnt=dist_changepnt,data=graph_seq,Sim_Lasso=Sim_Lasso))
}
