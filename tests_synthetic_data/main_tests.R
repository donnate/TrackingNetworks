######## In this module ar e

dir="/Users/cdonnat/Dropbox/TrackingNetworkChanges"
setwd(dir)
source("./tests_synthetic_data/test_functions.R")
source("./spanning_trees.R")
source("./distances.R")

########### Default parameters for random graph generation  #####################
N=30
p=0.4
#################################################################################















##########################################################################################################
##########################################################################################################
########### Functions for computing distances  for a sequence of evolving graphs #########################
##########################################################################################################
##########################################################################################################









##########################################################################################################
########### First set of functions: simple "smooth" evolution ############################################
##########################################################################################################



### -------------------------------------------------------------------------------------------------------
test_smooth_RD_changes<-function(N,p,prop,T=11,verbose=TRUE, go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, name_file_ext="",path_to_graph='./tests_synthetic_data/generated_graphs/',name_graph='',path2plot='./plots/'){

    ##  Description
    ##  -------------
    ##  Function generating an ER graph, and letting it evolve through time. Distances between consecutive graphs are computed
    ##  for a variety of distances (Jaccard, Hammong, IM, HIM, polynomial, ST)
    ##
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  p               :   probability of edge connection in an ER graph. (float/double<1)
    ##  prop            :   proportion of edges that are reassigned
    ##  T               :   length of the horizon (evolution period) (int>1)
    ##  verbose         :   bolean. should the method be verbose?
    ##  go_plot         :   should the algorithm print the dstance curves (bolean, Default: TRUE)
    ##  initial_message :   should a summary first message be printed? (bolean, Default: TRUE)
    ##  save_graph_seq  :   should the whole graph sequence be saved in a matrix? (bolean, Default: TRUE)
    ##  name_file_ext   :   string for the name of the graph seq. Default=""
    ##  path_to_graph   :   string indicating where the graph seq should be saved. Default='./tests_synthetic_data/generated_graphs/'
    ##  name_graph      :   name of the graph sequence. Default=''
    ##  path2plot       :   string indicating where the plots should be saved.  Default='./plots/'
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  results         : named list with arguments distances (vector of T distances between consecutive graphs) and graph_seq: an array of all the flattened       adjacency matrices of the graphs




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
    pdf(paste(path2plots,'smooth_rdm_changes_R.pdf',sep=''),width=6,height=4,paper='special')
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
    write.table(dist,sep=",",file = paste(path_to_graph,N,"_",p,"_",name_graph,"distances.txt",sep=""),row.names = TRUE,col.names = FALSE)
    write.table(graph_seq,sep=",",file = paste(path_to_graph,N,"_",p,"_",name_graph,".txt",sep=""),row.names = FALSE,col.names = FALSE)
  }
  
  return(list(dist=dist, data=graph_seq))
}
### -------------------------------------------------------------------------------------------------------









### -------------------------------------------------------------------------------------------------------
test_smooth_Realistic_changes<-function(N,m, p_disp=0.1,p_creation=0.01,T=11,opts=1,verbose=FALSE,go_plot=TRUE,initial_message=TRUE,save_graph_seq=T,path_to_graph='./tests_synthetic_data/generated_graphs/', name_file_ext="",very_verbose=F,path2plot='./plots/',args=list()){
    
    
    
    ##  Description
    ##  -------------
    ##  Function generating a more "realistic" random graph (PA, Island, SBM,WS), and letting it evolve through time. Different distances between consecutive graphs are computed
    ##
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  m               :   number of edges transformed between time t and t+1 (int)
    ##  p_disp          :   probability that one of the edges chosen as "transformed" disappears-- it is reassigned elswhere withprob 1-p_disp. (float/double<1)
    ##  T               :   length of the horizon (evolution period) (int>1)
    ##  opts            :   what type of random graph (1: PA, 2: Island, 3: Dot Product, 4: SBM). Default=1
    ##  verbose         :   bolean. should the method be verbose?
    ##  go_plot         :   should the algorithm print the dstance curves (bolean, Default: TRUE)
    ##  initial_message :   should a summary first message be printed? (bolean, Default: TRUE)
    ##  save_graph_seq  :   should the whole graph sequence be saved in a matrix? (bolean, Default: TRUE)
    ##  name_file_ext   :   string for the name of the graph seq. Default=""
    ##  path_to_graph   :   string indicating where the graph seq should be saved. Default='./tests_synthetic_data/generated_graphs/'
    ##  name_graph      :   name of the graph sequence. Default=''
    ##  path2plot       :   string indicating where the plots should be saved.  Default='./plots/'
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  results         : named list with arguments distances (vector of T distances between consecutive graphs) and graph_seq: an array of all the flattened       adjacency matrices of the graphs)
    
  
  ##### 1> Generate initial random matrix
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

  if (is.null(args) || length(args)==0){
    args<-list(power=0.9,islands.n=3,islands.size=10,islands.pin=0.3,n.inter=3,K=6,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
  }
  
  
  
  ##### 2> Make graph evolve and compute adjacent matrices
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















##########################################################################################################################################################
########### Second set of functions: Test evolution distance with drastic change point in evolution mechanism ############################################
##########################################################################################################################################################

### -------------------------------------------------------------------------------------------------------
### first function for the simple ER case :

test_change_point<-function(N=30,p0=0.4,p=0.4,prop=0.05,prop2=0.2,p2=0.4, T=21,verbose=FALSE,go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, path_to_graph='./tests_synthetic_data/generated_graphs/',name_file_ext="",path2plot='./plots/'){
    
    
    ##  Description
    ##  -------------
    ##  Function generating a more "realistic" random graph (PA, Island, SBM,WS), and letting it evolve through time. Different distances between consecutive graphs are computed
    ##
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  p0               :   probability of edge formation in initial ER graph (float/double<1)
    ##  prop            :   prop of edges that  are reassigned  in first dynamic regime (float/double<1)
    ##  prop2           :   prop of edges that  are reassigned  in second dynamic regime (float/double<1)
    ##  p               :   probability that one of the edges chosen as "transformed" disappears in the first dynamic regime -- it is reassigned elswhere withprob 1-p (float/double<1)
    ##  p2               :   probability that one of the edges chosen as "transformed" disappears in the first dynamic regime -- it is reassigned elswhere withprob 1-p2 (float/double<1)
    ##  T               :   length of the horizon (evolution period) (int>1)
    ##  verbose         :   bolean. should the method be verbose?
    ##  go_plot         :   should the algorithm print the dstance curves (bolean, Default: TRUE)
    ##  initial_message :   should a summary first message be printed? (bolean, Default: TRUE)
    ##  save_graph_seq  :   should the whole graph sequence be saved in a matrix? (bolean, Default: TRUE)
    ##  name_file_ext   :   string for the name of the graph seq. Default=""
    ##  path_to_graph   :   string indicating where the graph seq should be saved. Default='./tests_synthetic_data/generated_graphs/'
    ##  name_graph      :   name of the graph sequence. Default=''
    ##  path2plot       :   string indicating where the plots should be saved.  Default='./plots/'
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  results         : named list with arguments distances (vector of T distances between consecutive graphs) and graph_seq: an array of all the flattened       adjacency matrices of the graphs)
    
    
  if (initial_message==TRUE){
    print("Investigating the ability of the distance to detect dynamic regime changes")
    print(paste("An ER graph with N=", N, "nodes is considered, with edge proba=",p))
    print(paste("At time t=10, the proportion of nodes randomly reassigned changes from", prop, " to ",prop2, "with new connection proba= ", p2))
  }
  alphas=c(0.1,0.4,0.9,1.2,3)
  A<-generate_random_adjacency(N,p0)
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
### -------------------------------------------------------------------------------------------------------



### -------------------------------------------------------------------------------------------------------
### second function for more ``realistic'' random structures :

test_change_point_realistic<-function(N=30,m=10,m2=15,p=0.1,p2=0.1,p_disp=0.0,p_disp2=0.0,p_creation=0.01,opts=1, T=21,verbose=TRUE,m_disp=0,m_disp2=0,m_creation=5,compare_LASSO=FALSE,go_plot=TRUE,initial_message=TRUE,very_verbose=FALSE,save_graph_seq=FALSE,path_to_graph='./plots/',name_file_ext=""){
  if (initial_message==TRUE){
      
      
      
      ##  Description
      ##  -------------
      ##  Function generating a more "realistic" random graph (PA, Island, SBM,WS), and letting it evolve through time. Different distances between consecutive graphs are computed
      ##
      ##  INPUT:
      ##  =============================================================
      ##  N               :   Number of nodes in the graphs (int)
      ##  m               :   probability of edge formation in initial ER graph (float/double<1)
      ##  m               :   number of edges that  are reassigned  in first dynamic regime (float/double<1)
      ##  m2              :   number of edges that  are reassigned  in second dynamic regime (float/double<1)
      ##  p_disp          :   probability that one of the edges chosen as "transformed" disappears in the first dynamic regime -- it is reassigned elswhere withprob 1-p_disp  (float/double<1)
      ##  p_disp2         :   probability that one of the edges chosen as "transformed" disappears in the first dynamic regime -- it is reassigned elswhere withprob 1-p_disp2 (float/double<1)
      ##  p_creation      :   probability that an edge is created (for all time regimes)
      ##  T               :   length of the horizon (evolution period) (int>1)
      ##  opts            :   what type of random graph (1: PA, 2: Island, 3: Dot Product, 4: SBM). Default=1
      ##  verbose,very_verbose:   bolean. should the method be verbose?
      ##  m_disp          :   number of edges that disappear (with prob 1) in the first time regime (int. Default: 0)
      ##  m_disp2         :   number of edges that disappear (with prob 1) in the second time regime (int. Default: 0)
      ##  go_plot         :   should the algorithm print the dstance curves (bolean, Default: TRUE)
      ##  initial_message :   should a summary first message be printed? (bolean, Default: TRUE)
      ##  save_graph_seq  :   should the whole graph sequence be saved in a matrix? (bolean, Default: TRUE)
      ##  name_file_ext   :   string for the name of the graph seq. Default=""
      ##  path_to_graph   :   string indicating where the graph seq should be saved. Default='./tests_synthetic_data/generated_graphs/'
      ##  name_graph      :   name of the graph sequence. Default=''
      ##  path2plot       :   string indicating where the plots should be saved.  Default='./plots/'
      ##
      ##  OUTPUT
      ##  =============================================================
      ##  results         : named list with arguments distances (vector of T distances between consecutive graphs) and graph_seq: an array of all the flattened       adjacency matrices of the graphs)
    
    
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













