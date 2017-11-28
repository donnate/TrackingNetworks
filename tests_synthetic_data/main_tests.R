dir="/Users/cdonnat/Dropbox/Distances/distances"
setwd(dir)
source("~/Dropbox/Distances/tests_synthetic_data/test_functions.R")
source("spanning_trees.R")
source("distances.R")
source("../SC.R")
### generate adjacency matrix
N=30
p=0.4


#test_tree=get_spanning_trees(A, 1)


## function for visualizing the results
test_basic<-function(N,p,save_graph_seq=T, name_file_ext=""){
  
  B<-500
  print("Investigating the stability of the number of spanning trees:")
  print(paste("For B=",B,"trials, get number of spanning trees in a random Matrix with edge probability p=",p," and N=",N ,"nodes "))
  nb_spanning_trees<-matrix(0,B,1)
  for (i in 1:B){
    A<-generate_random_adjacency(N,p, TRUE)
    nb_spanning_trees[i]<-get_number_spanning_trees2(A)
    ## Note that all edges are weigthed equally in that case
  }
  summary(nb_spanning_trees)
  par(mfrow=c(1,2))
  plot(nb_spanning_trees)
  hist(nb_spanning_trees)
  if (save_graph_seq){
    write.table(graph_seq,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_basic_",N,"_",p,"_",name_file_ext,".txt",sep=""))
  }
  return(nb_spanning_trees)
}




#### Test evolution distances

test_smooth_RD_changes<-function(N,p,prop,T=11,verbose=TRUE, go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, name_file_ext="",path_to_graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/',name_graph=''){
  if (initial_message==TRUE){
    print("Investigating the distance for evolution of a graph over time with random changes:")
    print(paste("At each time iteration, ",100*prop, "% of nodes are randomly reassigned (ER-model)"))
  #print(paste("For B=",B,"trials, get number of spanning trees in a random Matrix with edge probability p=",p," and N=",N ,"nodes "))
  }
  alphas=c(0.1,0.4,0.9,1.2,3)
  A<-generate_random_adjacency(N,p, TRUE)
  dist<-matrix(0, 12+length(alphas)+1,T-1)
  graph_seq<-matrix(0,T,N^2)
  graph_seq[1,]<-as.vector(t(A))
  library(nettools)
  for (t in 1:(T-1)){
    if (verbose==TRUE) print(paste("time: ",t))
    A_new<-random_alteration_adjacency(A,prop,p)
    graph_seq[t+1,]<-as.vector(t(A_new))
    hamming<-hamming_based_distance(A, A_new)
    dist[1,t]<-hamming
    dist[2,t]<-jaccard_based_distance(A,A_new)
    sp_Anew=get_nb_ST(A_new)
    sp_A=get_nb_ST(A)
    #dist[3,t]<-abs(sp_Anew-sp_A)/(sp_A+sp_Anew)
    dist[3,t]=abs(log(sp_A)-log(sp_Anew))
    dist[10,t]<-2*abs(sp_Anew-sp_A)/(sp_A+sp_Anew)
    dist[4:6,t]<-netdist(A, A_new,d = "HIM")
    dist[7,t]<-poly_distance(A, A_new,order_max=5,alpha=0.5)
    dist[8,t]<-poly_distance(A, A_new,order_max=3,alpha=0.5)
    dist[9,t]<-poly_distance(A, A_new,order_max=3,alpha=0.9)
    dist[10,t]=ST_distance(A,A_new,norm=TRUE)
    dist[11,t]<-poly_distance2(A, A_new,order_max=5,alpha=0.5)
    dist[12,t]<-poly_distance2(A, A_new,order_max=3,alpha=0.9)
    it_alpha=0
    
    
    for (alpha in alphas){
      if (it_alpha==0){
        cols=paste("f(l)=exp(-",alpha,"l)",sep="")
      }
      else{
        cols=c(cols,paste("f(l)=exp(-",alpha,"l)",sep=""))
      }
      dist[13+it_alpha,t]<-eigen_distance(A, A_new,function(x){return(-alpha*x)},p=2)
      it_alpha=it_alpha+1
      
    }
    dist[13+it_alpha,t]<-eigen_distance(A, A_new,function(x){ifelse(x<2,x,0)})
    it_alpha=it_alpha+1
    cols=c(cols,"f(l)= l. 1(l<2)")
  }
  
  if(go_plot==TRUE){
    plot(1:(T-1), 1-dist[2,], type='l',col=2, xlab="time", ylab="distance", ylim=c(0,1),main="Distances for  small consecutive random changes",cex=0.7)
    for (j in 3:6){
      points(1:(T-1), dist[j,], type='l',col=j)
    }
    legend(1,1,c("Jaccard","Spanning Trees","Hamming","IM","HIM"),col=2:6,lty=1,lwd = 2,cex=0.7)
   
    
     plot(1:(T-1), dist[7,], type='l',col=7, xlab="time", ylab="distance", ylim=c(0,max(dist[7:9,])),main="Polynomial Distances for  small consecutive random changes",cex=0.7)
    for (j in 8:9){
      points(1:(T-1), dist[j,], type='l',col=j)
    }
    legend(1,max(dist[7:9,]),c("alpha=0.5,order=5","alpha=0.5,order=3","alpha=0.9,order=3",paste("poly 2 alpha=0.5,order=5 (x", signif(max(dist[11,]),4),")"),paste("poly2 alpha=0.9,order=3 (x", signif(max(dist[12,]),4),")")),col=7:9,lty=1,lwd = 2,cex=0.7)
    ### Global plot
    pdf('/Users/cdonnat/Dropbox/Distances/write_up/plot/smooth_rdm_changes_R.pdf',width=6,height=4,paper='special')
    plot(1:(T-1), 1-dist[2,], type='l',col=2, xlab="time", ylab="distance", ylim=c(0,1),main="Distances for  small consecutive random changes",cex=0.7)
    for (j in 3:12){
      if (j>7){
        points(1:(T-1), dist[j,]/max(dist[j,]), type='l',col=j)
      }
      else{
        points(1:(T-1), dist[j,], type='l',col=j)
      }
    }
    #legend("topright", inset=c(-0.2,0), legend=c("A","B"), pch=c(1,3), title="Group")
    par(xpd=TRUE)
    legend("topright", inset=c(-0.5,0), legend=c("Jaccard","Spanning Trees","Hamming","IM","HIM",paste("alpha=0.5,order=5 (x", signif(max(dist[7,]),4),")"),paste("alpha=0.5,order=3 (x", signif(max(dist[8,]),4),")"),paste("alpha=0.9,order=3 (x", signif(max(dist[9,]),4),")"),paste("ST norm (x", signif(max(dist[10,]),4),")"),paste("poly 2 alpha=0.5,order=5 (x", signif(max(dist[11,]),4),")"),paste("poly2 alpha=0.9,order=3 (x", signif(max(dist[12,]),4),")")),col=2:12,lty=1,lwd = 2,cex=0.7)
    #legend(1,max(dist[7:9,]),c("alpha=0.5,order=5","alpha=0.5,order=3","alpha=0.9,order=3"),col=7:9,lty=1,lwd = 2,cex=0.7)  }
    dev.off()
  }
  dist=data.frame(dist,row.names = c("Hamming1","Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=5","alpha=0.5,order=3","alpha=0.9,order=3","ST norm","poly2 alpha=0.5,order=5","poly2 alpha=0.9,order=3",cols))
  if (save_graph_seq){
    #write.csv(graph_seq,sep=" ",file = paste(path_to_graph,N,"_",p,"_",name_file_ext,".txt",sep=""),row.names = FALSE,col.names = FALSE)
    write.table(dist,sep=",",file = paste(path_to_graph,N,"_",p,"_",name_graph,"distances.txt",sep=""),row.names = TRUE,col.names = FALSE)
    write.table(graph_seq,sep=",",file = paste(path_to_graph,N,"_",p,"_",name_graph,".txt",sep=""),row.names = FALSE,col.names = FALSE)
  }
  
  return(list(dist=dist, data=graph_seq))
}



test_smooth_Realistic_changes<-function(N,m, p_disp=0.1,p_creation=0.01,T=11,opts=1,loc=FALSE,verbose=FALSE,go_plot=TRUE,initial_message=TRUE,save_graph_seq=T,path_to_graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/', name_file_ext="",very_verbose=F){
  Ag<-generate_realistic_adjacency(N,opts=opts, verbose=TRUE)
  A<-as(get.adjacency(Ag),"matrix")
  graph_seq<-matrix(0,T,nrow(A)^2)
  graph_seq[1,]<-as.vector(t(A))
  dist<-matrix(0, 18,T-1)
  
  if (initial_message==TRUE){
    print("Investigating the distance for evolution of a graph over time with random changes:")
    print(paste("A graph with N=", nrow(A), "nodes is considered as per creation model ",opts))
    print(paste("At each time iteration, ",m, "edges are either randomly reassigned or disppear with proba ",p_disp))
    print(paste("Probablity of creation a a new edge for each vertex= ",p_creation))
    #print(paste("For B=",B,"trials, get number of spanning trees in a random Matrix with edge probability p=",p," and N=",N ,"nodes "))
  }

  if (is.null(args)){
    args<-list(power=0.9,islands.n=3,islands.size=10,islands.pin=0.3,n.inter=3,K=6,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
  }
  
  
  for (t in 1:(T-1)){
    if (verbose==TRUE) print(paste("time: ",t))
    Ag_new<-graph_alteration(Ag,m,p_disp,p_creation)
    A_new=as(get.adjacency(Ag_new),"matrix")
    if (very_verbose) plot(Ag_new)
    graph_seq[t+1,]<-as.vector(t(A_new))
    hamming<-hamming_based_distance(A, A_new)
    dist[1,t]<-sum(abs(A-A_new))/(N*(N-1))
    dist[2,t]<-jaccard_based_distance(A,A_new)
    sp_Anew=get_nb_ST(A_new)
    sp_A=get_nb_ST(A)
    #dist[3,t]<-abs(sp_Anew-sp_A)/(sp_A+sp_Anew)
    dist[3,t]=abs(log(sp_A)-log(sp_Anew))
    #dist[3,t]<-abs(sp_Anew-sp_A)/(sp_A+sp_Anew)
    dist[4:6,t]<-netdist(A, A_new,d = "HIM")
    dist[7,t]<-poly_distance(A, A_new,order_max=3,alpha=0.5)
    dist[8,t]<-poly_distance(A, A_new,order_max=5,alpha=0.7)
    dist[9,t]<-poly_distance(A, A_new,order_max=3,alpha=0.9)
    dist[10,t]=abs(sp_A-sp_Anew)/(sp_A+sp_Anew)
    dist[11,t]<-poly_distance2(A, A_new,order_max=5,alpha=0.5)
    dist[12,t]<-poly_distance2(A, A_new,order_max=3,alpha=0.9)
    it_alpha=0
    alphas=c(0.1,0.4,0.9,1.2,3)
    for (alpha in alphas){
      if (it_alpha==0){
        cols=paste("f(l)=exp(-",alpha,"l)",sep="")
      }
      else{
        cols=c(cols,paste("f(l)=exp(-",alpha,"l)",sep=""))
      }
      dist[13+it_alpha,t]<-eigen_distance(A, A_new,function(x){return(-alpha*x)},p=2)
      it_alpha=it_alpha+1
      
    }
    dist[13+it_alpha,t]<-eigen_distance(A, A_new,function(x){ifelse(x<2,x,0)})
    it_alpha=it_alpha+1
    cols=c(cols,"f(l)= l. 1(l<2)")
    Ag<-Ag_new
  }
  dist=data.frame(dist,row.names = c("Hamming1","Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=3","alpha=0.7,order=3","alpha=0.9,order=3","ST norm","poly2 alpha=0.5,order=5","poly2 alpha=0.9,order=3",cols))
  max_poly=max(dist[7,])
  if (go_plot==TRUE){
    plot(1:(T-1), 1-dist[2,], type='l',col=2, xlab="time", ylab="distance", ylim=c(0,1),main="distances with a change of regime at t=10")
    for (j in 3:6){
      points(1:(T-1), dist[j,], type='l',col=j)
    }
    for (j in 7:12){
      points(1:(T-1), dist[j,]/max(dist[j,]), type='l',col=j)
    }
    par(xpd=TRUE)
    legend("topright", inset=c(-0.5,0), legend=c("Jaccard","Spanning Trees","Hamming","IM","HIM",paste("alpha=0.5,order=5 (x", signif(max(dist[7,]),4),")"),paste("alpha=0.5,order=3 (x", signif(max(dist[8,]),4),")"),paste("alpha=0.9,order=3 (x", signif(max(dist[9,]),4),")"),paste("ST norm (x", signif(max(dist[10,]),4),")"),paste("poly 2 alpha=0.5,order=5 (x", signif(max(dist[11,]),4),")"),paste("poly2 alpha=0.9,order=3 (x", signif(max(dist[12,]),4),")")),col=2:12,lty=1,lwd = 2,cex=0.7)
  }
  if (save_graph_seq){
    #dist=data.frame(dist,row.names = c("Hamming1","Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=3","alpha=0.7,order=5"))
    write.table(dist,sep=" ",file = paste(path_to_graph,"test_smooth_realistic_changes_",N,"_",opts,"_",name_file_ext,"distances.txt",sep=""),row.names = TRUE,col.names = FALSE)
    write.table(graph_seq,sep=" ",file = paste(path_to_graph,"test_smooth_realistic_changes_",N,"_",opts,"_",name_file_ext,".txt",sep=""),row.names = FALSE,col.names = FALSE)
  }
  return(list(dist=dist, data=graph_seq))
}




### Test evolution distance with drastic change point:

test_change_point<-function(N=30,p=0.4,prop=0.05,prop2=0.2,p2=0.4, T=21,verbose=FALSE,go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, path_to_graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/',name_file_ext=""){
  if (initial_message==TRUE){
    print("Investigating the ability of the distance to detect dynamic regime changes")
    print(paste("An ER graph with N=", N, "nodes is considered, with edge proba=",p))
    print(paste("At time t=10, the proportion of nodes randomly reassigned changes from", prop, " to ",prop2, "with new connection proba= ", p2))
  }
  alphas=c(0.1,0.4,0.9,1.2,3)
  A<-generate_random_adjacency(N,p)
  graph_seq<-matrix(0,T,N^2)
  graph_seq[1,]<-as.vector(t(A))
  dist<-matrix(0, 18,T-1)
  for (t in 1:(T-1)){
    if (verbose==TRUE) print(paste("time: ",t))
    if(t<10){
      A_new<-random_alteration_adjacency(A,prop,p)
    }
    else{
      A_new<-random_alteration_adjacency(A,prop2,p2)
    }
    graph_seq[t+1,]<-as.vector(t(A_new))
    hamming<-hamming_based_distance(A, A_new)
    dist[1,t]<-1/(N*(N-1))*sum(abs(A-A_new))
    dist[2,t]<-jaccard_based_distance(A,A_new)
    sp_Anew=get_nb_ST(A_new)
    sp_A=get_nb_ST(A)
    #dist[3,t]<-abs(sp_Anew-sp_A)/(sp_A+sp_Anew)
    dist[3,t]=abs(log(sp_A)-log(sp_Anew))
    dist[4:6,t]<-netdist(A, A_new,d = "HIM")
    dist[7,t]<-poly_distance(A, A_new,order_max=3,alpha=0.5)
    dist[8,t]<-poly_distance(A, A_new,order_max=5,alpha=0.7)
    dist[9,t]<-poly_distance(A, A_new,order_max=3,alpha=0.9)
    dist[10,t]=abs(sp_Anew-sp_A)/(sp_A+sp_Anew)
    dist[11,t]<-poly_distance2(A, A_new,order_max=5,alpha=0.5)
    dist[12,t]<-poly_distance2(A, A_new,order_max=3,alpha=0.9)
    it_alpha=0
    
    for (alpha in alphas){
      if (it_alpha==0){
        cols=paste("f(l)=exp(-",alpha,"l)",sep="")
      }
      else{
        cols=c(cols,paste("f(l)=exp(-",alpha,"l)",sep=""))
      }
      dist[13+it_alpha,t]<-eigen_distance(A, A_new,function(x){return(-alpha*x)},p=2)
      it_alpha=it_alpha+1
      
    }
    dist[13+it_alpha,t]<-eigen_distance(A, A_new,function(x){ifelse(x<2,x,0)})
    it_alpha=it_alpha+1
    cols=c(cols,"f(l)= l. 1(l<2)")
    A<-A_new
  }
  dist=data.frame(dist,row.names = c("Hamming1","Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=5","alpha=0.5,order=3","alpha=0.9,order=3","ST norm","poly2 alpha=0.5,order=5","poly2 alpha=0.9,order=3",cols))
  
  
  max_poly=max(dist[7,])
  if (go_plot==TRUE){
    plot(1:(T-1), 1-dist[2,], type='l',col=2, xlab="time", ylab="distance", ylim=c(0,1),main="distances with a change of regime at t=10")
    for (j in 3:6){
      points(1:(T-1), dist[j,], type='l',col=j)
    }
    for (j in 7:12){
      points(1:(T-1), dist[j,]/max(dist[j,]), type='l',col=j)
    }
    par(xpd=TRUE)
    legend("topright", inset=c(-0.5,0), legend=c("Jaccard","Spanning Trees","Hamming","IM","HIM",paste("alpha=0.5,order=5 (x", signif(max(dist[7,]),4),")"),paste("alpha=0.5,order=3 (x", signif(max(dist[8,]),4),")"),paste("alpha=0.9,order=3 (x", signif(max(dist[9,]),4),")"),paste("ST norm (x", signif(max(dist[10,]),4),")"),paste("poly 2 alpha=0.5,order=5 (x", signif(max(dist[11,]),4),")"),paste("poly2 alpha=0.9,order=3 (x", signif(max(dist[12,]),4),")")),col=2:12,lty=1,lwd = 2,cex=0.7)
  }
  if (save_graph_seq){
    #dist=data.frame(dist,row.names = c("Hamming1","Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=3","alpha=0.7,order=5"))
    write.table(dist,sep=" ",file = paste(path_to_graph,"test_change_points_",N,"_",p,"_",name_file_ext,"distances.txt",sep=""),row.names = TRUE,col.names = FALSE)
    write.table(graph_seq,sep=" ",file = paste(path_to_graph,"test_change_points_",N,"_",p,"_",name_file_ext,".txt",sep=""),row.names = FALSE,col.names = FALSE)
  }
  return(list(dist=dist, data=graph_seq))
}

test_change_point_realistic<-function(N=30,m=10,m2=15,p=0.1,p2=0.1,p_disp=0.0,p_disp2=0.0,p_creation=0.01,opts=1, T=21,verbose=TRUE,m_disp=0,m_disp2=0,m_creation=5,compare_LASSO=FALSE,go_plot=TRUE,initial_message=TRUE,very_verbose=FALSE,save_graph_seq=FALSE,path_to_graph='/Users/cdonnat/Dropbox/Distances/tests_synthetic_data/generated_graphs/',name_file_ext=""){
  if (initial_message==TRUE){
    print("Investigating the ability of the distance to detect dynamic regime changes")
    
    print(paste("A graph with N=", N, "nodes is considered as per creation model ",opts))
    print(paste("At time t=10, the number of randomly reassigned edges changes from", m, " to ",m2, "with new change proba= ", p2, " vs ", p, "and deletion",p_disp ))
  }
  alphas=c(0.1,0.4,0.9,1.2,3)
  Ag<-generate_realistic_adjacency(N,opts=opts, verbose=verbose)
  A<-as(get.adjacency(Ag),"matrix")
  graph_seq<-matrix(0,T,nrow(A)^2)
  graph_seq[1,]<-as.vector(t(A))
  dist<-matrix(0, 18,T-1)
  for (t in 1:(T-1)){
    if (verbose==TRUE) print(paste("time: ",t))
    if(t<7){
      Ag_new<-graph_alteration(Ag,m=m,p=p,p_disp=p_disp,p_creation=p_creation,m_disp=m_disp,m_creation=m_creation)
    }
    else{
      if (t<14) Ag_new<-graph_alteration(Ag,m=m2,p=p2,p_disp=p_disp2,p_creation=p_creation,m_disp=m_disp2,m_creation=m_creation)
      else{
        Ag_new<-graph_alteration(Ag,m=m,p=p,p_disp=p_disp,p_creation=p_creation,m_disp=m_disp,m_creation=m_creation)
        print("changed back")
      }
    }
    if (very_verbose) plot(Ag_new)
    A_new=as(get.adjacency(Ag_new),"matrix")
    graph_seq[t+1,]<-as.vector(t(A_new))
    hamming<-hamming_based_distance(A, A_new)
    dist[1,t]<-1/(N*(N-1))*sum(abs(A-A_new))
    dist[2,t]<-jaccard_based_distance(A,A_new)
    sp_Anew=get_nb_ST(A_new)
    sp_A=get_nb_ST(A)
    #dist[3,t]<-abs(sp_Anew-sp_A)/(sp_A+sp_Anew)
    dist[3,t]=abs(log(sp_A)-log(sp_Anew))
    dist[4:6,t]<-netdist(A, A_new,d = "HIM")
    dist[7,t]<-poly_distance(A, A_new,order_max=3,alpha=0.5)
    dist[8,t]<-poly_distance(A, A_new,order_max=5,alpha=0.7)
    dist[8,t]<-poly_distance(A, A_new,order_max=5,alpha=0.7)
    dist[9,t]<-poly_distance(A, A_new,order_max=3,alpha=0.9)
    dist[10,t]=abs(sp_Anew-sp_A)/(sp_A+sp_Anew)
    dist[11,t]<-poly_distance2(A, A_new,order_max=5,alpha=0.5)
    dist[12,t]<-poly_distance2(A, A_new,order_max=3,alpha=0.9)
    it_alpha=0
    
    for (alpha in alphas){
      if (it_alpha==0){
        cols=paste("f(l)=exp(-",alpha,"l)",sep="")
      }
      else{
        cols=c(cols,paste("f(l)=exp(-",alpha,"l)",sep=""))
      }
      dist[13+it_alpha,t]<-eigen_distance(A, A_new,function(x){return(-alpha*x)},p=2)
      it_alpha=it_alpha+1
      
    }
    dist[13+it_alpha,t]<-eigen_distance(A, A_new,function(x){ifelse(x<2,x,0)})
    it_alpha=it_alpha+1
    cols=c(cols,"f(l)= l. 1(l<2)")
    A<-A_new
  }
  dist=data.frame(dist,row.names = c("Hamming1","Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=5","alpha=0.5,order=3","alpha=0.9,order=3","ST norm","poly2 alpha=0.5,order=5","poly2 alpha=0.9,order=3",cols))
  
  max_poly=max(dist[7,])
  if(compare_LASSO==TRUE){
    Sim_Lasso=get_Similarity(graph_seq)
    Sim_Lasso=Sim_Lasso+t(Sim_Lasso)
    if (go_plot==TRUE){
      heatmap.2(Sim_Lasso,Rowv=FALSE,Colv=FALSE,dendrogram = "none")
    }
  }
  else{
    Sim_Lasso=FALSE
  }
  if (go_plot==TRUE){
    plot(1:(T-1), 1-dist[2,], type='l',col=2, xlab="time", ylab="distance", ylim=c(0,1),main="distances with a change of regime at t=10")
    for (j in 3:6){
      points(1:(T-1), dist[j,], type='l',col=j)
    }
    for (j in 7:12){
      points(1:(T-1), dist[j,]/max(dist[j,]), type='l',col=j)
    }
    par(xpd=TRUE)
    legend("topright", inset=c(-0.5,0), legend=c("Jaccard","Spanning Trees","Hamming","IM","HIM",paste("alpha=0.5,order=5 (x", signif(max(dist[7,]),4),")"),paste("alpha=0.5,order=3 (x", signif(max(dist[8,]),4),")"),paste("alpha=0.9,order=3 (x", signif(max(dist[9,]),4),")"),paste("ST norm (x", signif(max(dist[10,]),4),")"),paste("poly 2 alpha=0.5,order=5 (x", signif(max(dist[11,]),4),")"),paste("poly2 alpha=0.9,order=3 (x", signif(max(dist[12,]),4),")")),col=2:12,lty=1,lwd = 2,cex=0.7)
  }
  #dist=data.frame(dist,row.names = c("Hamming1","Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=3","alpha=0.7,order=5","alpha=0.9,order=3","norm ST","poly2 alpha=0.5,order=5","poly2 alpha=0.9,order=3"))
  
  if (save_graph_seq){
    #dist=data.frame(dist,row.names = c("Hamming1","Jaccard","Spanning Trees","Hamming","IM","HIM","alpha=0.5,order=3","alpha=0.7,order=5","alpha=0.9,order=3","norm ST","poly2 alpha=0.5,order=5","poly2 alpha=0.9,order=3"))
    lab=c("pa","island","dot","SBM")
    write.table(dist,sep=" ",file = paste(path_to_graph,"test_change_point_realistic_",N,"_",lab[opts],"_",name_file_ext,"distances.txt",sep=""),row.names = TRUE,col.names = FALSE)
    write.table(graph_seq,sep=" ",file = paste(path_to_graph,"test_change_point_realistic_",N,"_",lab[opts],"_",name_file_ext,".txt",sep=""),row.names = FALSE,col.names = FALSE)
  }
  return(list(dist=dist, data=graph_seq,Sim_Lasso=Sim_Lasso))
}


### What seems to come out of this study: the Spanning tree distance definitey needs some improvement
### it seems very unstable as it is. Maybe it's all wrong though, who knows

tree_test<-function(){
  print("test a simple tree-like structure to see if the nb of spanning trees is correctly determined")
  tree<-matrix(0,15,15)
  tree[1,2:3]<-1
  tree[2,4:5]<-1
  tree[3,6:7]<-1
  tree[4,8:9]<-1
  tree[5,10:11]<-1
  tree[6,12:13]<-1
  tree[7,14:15]<-1
  tree<-tree+t(tree)
  print(paste("method 2: ", sp_tree=get_number_spanning_trees2(tree), "spanning trees"))
  print(paste("method 1: ", sp_tree=get_number_spanning_trees(tree), "spanning trees"))
  ### One added edge
  tree2<-tree
  tree2[14,15]<-1
  tree2[15,14]<-1
  print(paste("method 2: ", sp_tree=get_number_spanning_trees2(tree), "spanning trees"))
  print(paste("method 1: ", sp_tree=get_number_spanning_trees(tree), "spanning trees"))
  print(paste("custom distances: ",get_number_spanning_trees3(tree, tree2)))
  print(paste("HIM distances: ",netdist(tree, tree2)))
  tree<-matrix(0,16,16)
  tree[1,2:3]<-1
  tree[2,4:5]<-1
  tree[3,6:7]<-1
  tree[4,8:9]<-1
  tree[5,10:11]<-1
  tree[6,12:13]<-1
  tree[7,14:15]<-1
  
}





##### Design other tests to assess the pertinence of the metrics proposed
##### Note that we would have to prove in each case that these are indeed metrics.


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




##### Tests on real data: the Stock Options dataset"
stock_option_test<-function(){
  load("/Users/cdonnat/Dropbox/Sto/Code/seq_graph_10days.RData")
  K<-length(graph_seq)
  
  
  
  g<-list()
  thres=0.55
  for (k in 1:K){
    g[[k]]<-sapply(corr_full[[k]],FUN=function(x){
      return(prune(x,thres))
    })
  }
  
  for (k in 1:K){
    #graph_seq[[k]]<-matrix(graph_seq[[k]],nrow(corr_full[[k]]),ncol(corr_full[[k]]))
    if (length(g[[k]])<ncol(stock_data)^2){
      temp<-data.frame(matrix(0,ncol(stock_data),ncol(stock_data)),row.names=colnames(stock_data))
      colnames(temp)=colnames(stock_data)
      temp[colnames(corr_full[[k]]),colnames(corr_full[[k]])]<-sapply(g[[k]],FUN=function(x){
        return(delta_c(x,0))
      })
      g[[k]]<-temp
    }
    else{
      g[[k]]<-data.frame(matrix(sapply(g[[k]],FUN=function(x){ return(delta_c(x,0))}),ncol(stock_data),ncol(stock_data)),row.names=row.names(corr_full[[k]]))
      colnames(g[[k]])<-colnames(corr_full[[k]])
    }
    
  }
  return(apply_distances(g,try_lasso=TRUE, file_dir="/Users/cdonnat/Dropbox/Food_network/code/RData/save_Sim_stocks.RData"))
}

apply_distances<-function(graph_seq,try_lasso=FALSE,file_dir){
  K=length(graph_seq)
  Sim1<-matrix(0,K,K)
  Sim2<-matrix(0,K,K)
  Sim3<-matrix(0,K,K)
  if(try_lasso==TRUE){
    Sim4<-matrix(0,K,K)
    Sim5<-matrix(0,K,K)
    Sim6<-matrix(0,K,K)
  }
  for (k in 1:(K-1)){
    #active_nodes_k=which(apply(g1,1,sum)>0)
    print(k)
    
    for( kk in (k+1):K){
      #active_nodes_kk=which(apply(g2,1,sum)>0)
      #candidate_nodes=union(active_nodes_k,active_nodes_kk)
      #g1_temp<-g1[candidate_nodes,candidate_nodes]
      #g2_temp<-g2[candidate_nodes,candidate_nodes]
      dist=netdist(graph_seq[[k]],graph_seq[[kk]],d = "HIM")
      print(dist)
      Sim1[k,kk]<-dist[1]
      Sim2[k,kk]<-dist[2]
      Sim3[k,kk]<-dist[3]
      
    }
    save(Sim1,Sim2,Sim3,file=file_dir)
  }
  res=list(Sim1=Sim1+t(Sim1),Sim2=Sim2+t(Sim2),Sim3=Sim3+t(Sim3) )
  if(try_lass0==TRUE){
    test_LASSO=test_lasso_distance(graph_seq,method="LASSO")
    test_ADMM1=test_lasso_distance(graph_seq,method="ADMM1")
    test_ADMM2=test_lasso_distance(graph_seq,method="ADMM2")
    res=list(Sim1=Sim1+t(Sim1),Sim2=Sim2+t(Sim2),Sim3=Sim3+t(Sim3),test_LASSO=test_LASSO,test_ADMM1=test_ADMM1,test_ADMM2=test_ADMM2 )
  }
  return(res)
}



test_bootstrap<-function(N,args,opts=1,B){
  boot_samples<-sapply(1:B, FUN=function(b){
    if (opts==0){
      p=args$p
      A<-generate_random_adjacency(N,p, TRUE)
      A2<-generate_random_adjacency(N,p, TRUE)
      
    }
    else{
      A<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
      A<-get.adjacency(A)
      A2<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
      A2<-get.adjacency(A2)
    }
    
    dist<-matrix(0,6,1)
    dist[1]=poly_distance(A, A2,order_max=5,alpha=0.5)
    dist[2]=poly_distance(A, A2,order_max=5,alpha=1.2)
    dist[3]=poly_distance(A, A2,order_max=3,alpha=0.5)
    dist[4]=poly_distance(A, A2,order_max=3,alpha=0.9)
    dist[5]=poly_distance(A, A2,order_max=3,alpha=1)
    dist[6]=poly_distance(A, A2,order_max=1,alpha=1)
    return(dist)
  })
}

test_bootstrap_shiny<-function(N,alpha,order_max,args,opts=1,B){
  boot_samples<-sapply(1:B, FUN=function(b){
    if (opts==0){
      p=args$p
      A<-generate_random_adjacency(N,p, TRUE)
      A2<-generate_random_adjacency(N,p, TRUE)
      
    }
    else{
      A<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
      A<-as.matrix(get.adjacency(A))
      A2<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
      A2<-as.matrix(get.adjacency(A2))
    }
    dist=poly_distance(A, A2,order_max,alpha=alpha)
    return(dist)
  })
}

test_bootstrap_shiny2<-function(N,alpha,order_max,args,opts=1,opts2=2,B){
  boot_samples<-sapply(1:B, FUN=function(b){
    if (opts==0){
      p=args$p
      A<-generate_random_adjacency(N,p, TRUE)
    }
    else{
      A<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
      A<-as.matrix(get.adjacency(A))
    }
    if (opts==0){
      p=args$p
      A2<-generate_random_adjacency(N,p, TRUE)
      
    }
    else{
      A2<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
      A2<-as.matrix(get.adjacency(A2))
    }
    dist=poly_distance(A, A2,order_max,alpha=alpha)
    return(dist)
  })
}


test_bootstrap_shiny_ST<-function(N,alpha,order_max,args,opts=1,B){
  boot_samples<-sapply(1:B, FUN=function(b){
    if (opts==0){
      p=args$p
      A<-generate_random_adjacency(N,p, TRUE)
      A2<-generate_random_adjacency(N,p, TRUE)
      
    }
    else{
      A<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
      A<-as.matrix(get.adjacency(A))
      A2<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
      A2<-as.matrix(get.adjacency(A2))
      
    }
    stree<-get_number_spanning_trees2(A)
    stree2<-get_number_spanning_trees2(A2)
    dist=abs(stree-stree2)/(stree+stree2)
    return(dist)
  })
}


test_bootstrap_shiny_HIM<-function(N,alpha,order_max,args,opts=1,B){
  boot_samples<-sapply(1:B, FUN=function(b){
    if (opts==0){
      p=args$p
      A<-generate_random_adjacency(N,p, TRUE)
      A2<-generate_random_adjacency(N,p, TRUE)
      
    }
    else{
      A<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
      A<-as.matrix(get.adjacency(A))
      A2<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
      A2<-as.matrix(get.adjacency(A2))
      
    }
    dist=netdist(A, A2,d = "HIM")
    return(dist)
  })
}

tb<-function(){
  args<-list(power=0.9,islands.n=3,islands.size=10,islands.pin=0.3,n.inter=3,K=3,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
  N=30
  opts=1
  B=1000
  test_bootstrap1<-test_bootstrap(N,args,opts,B)
}

compare_spanning_trees<-function(graph_seq,args=list(order_max=3, alpha=1)){
    s=dim(graph_seq)
    N=sqrt(s[2])
    dist1<-matrix(0, s[1], s[1])
    dist2<-matrix(0, s[1], s[1])
    dist3<-matrix(0, s[1], s[1])
    dist4<-matrix(0, s[1], s[1])
    dist5<-matrix(0, s[1], s[1])
    for (i in 1:(s[1]-1)){
        dist5[i,(i+1):s[1]]<-sapply((i+1):s[1], FUN=function(j){
            return(poly_distance(matrix(graph_seq[i,], N,N),matrix(graph_seq[j,], N,N),order_max=args$order_max,alpha=args$alpha))
        })
        dist1[i,(i+1):s[1]]<-sapply((i+1):s[1], FUN=function(j){
                return(ST_distance(matrix(graph_seq[i,], N,N),matrix(graph_seq[j,], N,N)))
            })
        temp<-sapply((i+1):s[1], FUN=function(j){
                dist=netdist(matrix(graph_seq[i,], N,N),matrix(graph_seq[j,], N,N),d="HIM")
                return(dist)
            })
        dist2[i,(i+1):s[1]]=temp[1,]
        dist3[i,(i+1):s[1]]=temp[2,]
        dist4[i,(i+1):s[1]]=temp[3,]
    }
    
  g1<-graph_from_adjacency_matrix(dist1,weighted=TRUE,mode="upper")
  mst1<-mst(g1,weights = edge_attr(g1, "weight")) 
  g2<-graph_from_adjacency_matrix(dist2,weighted=TRUE,mode="upper")
  mst2<-mst(g2,weights = edge_attr(g2, "weight") )
  g3<-graph_from_adjacency_matrix(dist3,weighted=TRUE,mode="upper")
  mst3<-mst(g3,weights = edge_attr(g3, "weight") )
  g4<-graph_from_adjacency_matrix(dist4,weighted=TRUE,mode="upper")
  mst4<-mst(g4,weights = edge_attr(g4, "weight") )
  g5<-graph_from_adjacency_matrix(dist5,weighted=TRUE,mode="upper")
  mst5<-mst(g5,weights = edge_attr(g5, "weight") )
  Sim_Lasso=get_Similarity(graph_seq)
  Sim_Lasso=Sim_Lasso+t(Sim_Lasso)
  g6<-graph_from_adjacency_matrix(Sim_Lasso,weighted=TRUE,mode="undirected")
  go_plot=TRUE
  if (go_plot==TRUE){
    plot(mst(g1,weights = edge_attr(g1, "weight") ), vertex.size=5,main="Spanning Tree from Poly")
    plot(mst(g2,weights = edge_attr(g2, "weight") ), vertex.size=5,main="Spanning Tree from Hamming")
    plot(mst(g3,weights = edge_attr(g3, "weight") ), vertex.size=5,main="Spanning Tree from IM")
    plot(mst(g4,weights = edge_attr(g4, "weight") ), vertex.size=5,main="Spanning Tree from HIM")
    plot(mst(g5,weights = edge_attr(g5, "weight") ), vertex.size=5,main="Spanning Tree from ST Sim")
    plot(mst(g6,weights = -edge_attr(g6, "weight") ), vertex.size=5,main="Spanning Tree from Lasso Sim")
  }
  
  
  return(list(Sim_Lasso=Sim_Lasso, dist1=dist1,dist2=dist2,dist3=dist3,dist4=dist4,dist5=dist5,mst1=mst1,mst2=mst2,mst3=mst3,mst4=mst4,mst5=mst5))
}
