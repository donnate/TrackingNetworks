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














########################################################################################################
########################################################################################################
########################################################################################################
########### Functions for computing distances  for a sequence of evolving graphs #######################
########################################################################################################
########################################################################################################
########################################################################################################







########################################################################################################
########### First set of functions: simple "smooth" evolution ##########################################
########################################################################################################



### -------------------------------------------------------------------------------------------------------
test_smooth_RD_changes<-function(N,p0,p,p_disp,p_creation,prop,T=11,alphas=c(0.1,0.4,0.9,1.2,3),verbose=TRUE, go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, name_file_ext="",path_to_graph='./tests_synthetic_data/generated_graphs/',name_graph='',path2plot='./plots/'){

    ##  Description
    ##  -------------
    ##  Function generating an ER graph, and letting it evolve through time. Distances between consecutive graphs are computed
    ##  for a variety of distances (Jaccard, Hammong, IM, HIM, polynomial, ST)
    ##
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  p0              :   parameter for initial ER graph: probability of edge connection in an ER graph. (float/double<1)
    ##  prop            :   proportion of the edges that are potentially modified
    ##  p               :   probability that a candidate ''modifiable'' edge is reassigned somewhere else
    ##  p_disp          :   probability that a candidate ''modifiable'' edge is deleted (p+p_disp<1)
    ##  p_creation      :   probability that a "non-existing" edge is added to the new graph
    ##  T               :   length of the horizon (evolution period) (int>1)
    ##  alphas          :   list of parameters to try for the exponential eigenvalue distances (f(lambda)=exp(-alpha*lambda))
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
  }
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
test_smooth_Realistic_changes<-function(N,m,m_disp=0,m_creation=0, p_modified, p_disp=0.1,p_creation=0.01,T=11,opts=1,alphas=c(0.1,0.4,0.9,1.2,3),verbose=FALSE,go_plot=TRUE,initial_message=TRUE,save_graph_seq=T,path_to_graph='./tests_synthetic_data/generated_graphs/', name_file_ext="",very_verbose=F,path2plot='./plots/',...){
    
    
    
    ##  Description
    ##  -------------
    ##  Function generating a more "realistic" random graph (PA, Island, SBM,WS), and letting it evolve through time.
    ##  The evolution mechanism is the following: we chose m edges: these selected edges are
    ##  potentially modifiable edges. Each candidate edge is then rewired with probability p, or disappears with probability p_disp, and stays in place
    ##  with probability 1-p_disp-p. We also allow for the possibility of edges to be created with probability p_creation
    ##  Different distances between consecutive graphs are computed
    ##
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  m               :   number of edges that can be potentially rewired between time t and t+1 (int)
    ##  m_disp          :   number of edges that can be potentially deleted between time t and t+1 (int)
    ##  m_created       :   number of edges that can be potentially created between time t and t+1 (int)
    ##  p_modified      :   probability that one of the m modifiable edges is rewired (float/double<1)
    ##  p_disp          :   probability that one of the edges chosen as "transformed" disappears-- it is reassigned elswhere withprob 1-p_disp. (float/double<1)
    ##  p_creation      :   probability that one "inexisting" edge appears.
    ##  alphas          :   list of parameters to try for the exponential eigenvalue distances (f(lambda)=exp(-alpha*lambda))
    ##  T               :   length of the horizon (evolution period) (int>1)
    ##  opts            :   what type of random graph (0: ER, 1: PA, 2: Island, 3: Dot Product, 4: SBM). Default=1
    ##  verbose         :   bolean. should the method be verbose?
    ##  go_plot         :   should the algorithm print the dstance curves (bolean, Default: TRUE)
    ##  initial_message :   should a summary first message be printed? (bolean, Default: TRUE)
    ##  save_graph_seq  :   should the whole graph sequence be saved in a matrix? (bolean, Default: TRUE)
    ##  name_file_ext   :   string for the name of the graph seq. Default=""
    ##  path_to_graph   :   string indicating where the graph seq should be saved. Default='./tests_synthetic_data/generated_graphs/'
    ##  name_graph      :   name of the graph sequence. Default=''
    ##  path2plot       :   string indicating where the plots should be saved.  Default='./plots/'
    ##  ...             :   additional named parameters for initial graph generation (as per igraph nomenclature)
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  results         : named list with arguments distances (vector of T distances between consecutive graphs) and graph_seq: an array of all the flattened       adjacency matrices of the graphs)
    
  
  ##### 1> Generate initial random matrix
  args<-list(p=0.1,power=0.9,islands.n=3,islands.size=9,islands.pin=0.3,n.inter=3,K=6,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
  ### Change argument list according to what is given as argument to the function
  if (hasArg(p) )args$p=p
  if (hasArg(power) )args$power=power
  if (hasArg(islands.n)) args$islands.n=islands.n
  if (hasArg(islands.size) )args$islands.size=islands.size
  if (hasArg(islands.pin)) args$islands.pin=islands.pin
  if (hasArg(n.inter)) args$n.inter=n.inter
  if (hasArg(K)) args$K=K
  if (hasArg(block.sizes)) args$block.sizes=block.sizes
  if (hasArg(pm)) args$pm=pm
  
  Ag<-generate_realistic_adjacency(N,opts=opts, args=args, verbose=TRUE)
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


  
  
  ##### 2> Make graph evolve and compute adjacent matrices
  for (t in 1:(T-1)){
    if (verbose==TRUE) print(paste("time: ",t))
    Ag_new<-graph_alteration(Ag,m,p_modified, p_disp,p_creation,m_disp=m_disp,m_creation=m_creation)
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
##########################################################################################################################################################
########### Second set of functions: Test evolution distance with drastic change point in evolution mechanism ############################################
##########################################################################################################################################################
##########################################################################################################################################################






### -------------------------------------------------------------------------------------------------------
### first function for the simple ER case :

test_change_point<-function(N=30,p0=0.4,p=0.4,prop=0.05,prop2=0.2,p2=0.4, T=21,verbose=FALSE,go_plot=TRUE,initial_message=TRUE,save_graph_seq=T, path_to_graph='./tests_synthetic_data/generated_graphs/',name_file_ext="",path2plot='./plots/'){
    
    
    ##  Description
    ##  -------------
    ##  Function generating a more "realistic" random graph (PA, Island, SBM,WS), and letting it evolve through time.
    ##  Again, the evolution mechanism is the following: we chose prop % of the edges: these selected edges are
    ##  potentially modifiable edges. Each candidate edge is then rewired with probability p, and stays in place
    ##  with probability 1-p.
    ##  Different distances between consecutive graphs are computed
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

test_change_point_realistic<-function(N=30,m=c(10,15),p_mod=c(0.1,0.1),p_disp=c(0.0,0.0),p_creation=c(0.0,0.0),opts=1, T=21,verbose=TRUE,m_disp=0,m_disp2=0,m_creation=5,go_plot=TRUE,initial_message=TRUE,very_verbose=FALSE,save_graph_seq=FALSE,path_to_graph='./plots/',name_file_ext="",...){
  if (initial_message==TRUE){
      
      
      
      ##  Description
      ##  -------------
      ##  Function generating a more "realistic" random graph (PA, Island, SBM,WS), and letting it evolve through time.
      ##  The evolution mechanism is the following: we chose m edges: these selected edges are
      ##  potentially modifiable edges. Each candidate edge is then rewired with probability p, or disappears with probability p_disp, and stays in place
      ##  with probability 1-p_disp-p. Here, we add a change point in the evolution mechanism. At time t=10, the number of modifiable edges changes to m2,
      ##  the rewiring probabiity changes from p to p2, and the disappearance probability changes from p_disp to p_disp2.
      ##  Different distances between consecutive graphs are computed.
      ##
      ##  INPUT:
      ##  =============================================================
      ##  N               :   Number of nodes in the graphs (int)
      ##  m               :   2-d vector for the number of edges that  potentially are rewired (int)
      ##  m_disp          :   2-d vector for the number of edges that  potentially are deleted in firs/second  dynamic regime (int)
      ##  m_creation      :   2-d vector for the number of edges that  are created  in first/second  dynamic regime (int)
      ##  p_mod           :   2-d vector  for the probability that one of the edges chosen as "potentially rewired" is rewired in the first dynamic regime
      ##                      (float/double<1)
      ##  p_disp          :   2-d vector  for the probability that one of the edges chosen as "potentially deleted" disappears in the first dynamic regime
      ##                      (float/double<1)
      ##  p_creation      :   2-d vector  for the probability that an edge is created (float/double<1)
      ##  T               :   length of the horizon (evolution period) (int>1)
      ##  opts            :   what type of random graph (0: ER, 1: PA, 2: Island, 3: Dot Product, 4: SBM). Default=1
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
      ##  ...             :   additional named parameters for the initial graph generation (as per igraph nomenclature)
      ##
      ##  OUTPUT
      ##  =============================================================
      ##  results         : named list with arguments distances (vector of T distances between consecutive graphs) and graph_seq: an array of all the flattened       adjacency matrices of the graphs)
    
    
    print("Investigating the ability of the distance to detect dynamic regime changes")
    
    print(paste("A graph with N=", N, "nodes is considered as per creation model ",opts))
    print(paste("At time t=10, the number of randomly reassigned edges changes from", m, " to ",m2, "with new change proba= ", p2, " vs ", p, "and deletion",p_disp ))
  }
  alphas=c(0.1,0.4,0.9,1.2,3)
  ##### 1> Generate initial random matrix
  args<-list(p=0.1,power=0.9,islands.n=3,islands.size=9,islands.pin=0.3,n.inter=3,K=6,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
  ### Change argument list according to what is given as argument to the function
  if (hasArg(p) )args$p=p
  if (hasArg(power) )args$power=power
  if (hasArg(islands.n)) args$islands.n=islands.n
  if (hasArg(islands.size) )args$islands.size=islands.size
  if (hasArg(islands.pin)) args$islands.pin=islands.pin
  if (hasArg(n.inter)) args$n.inter=n.inter
  if (hasArg(K)) args$K=K
  if (hasArg(block.sizes)) args$block.sizes=block.sizes
  if (hasArg(pm)) args$pm=pm
  Ag<-generate_realistic_adjacency(N,opts=opts, verbose=verbose)
  A<-as(get.adjacency(Ag),"matrix")
  graph_seq<-matrix(0,T,nrow(A)^2)
  graph_seq[1,]<-as.vector(t(A))
  dist<-matrix(0, 18,T-1)
  for (t in 1:(T-1)){
    if (verbose==TRUE) print(paste("time: ",t))
    if(t<7){
      Ag_new<-graph_alteration(Ag,m=m[1],p=p_mod[1],p_disp=p_disp[1],p_creation=p_creation[1],m_disp=m_disp[1],m_creation=m_creation[1])
    }
    else{
      if (t<14) Ag_new<-graph_alteration(Ag,m=m[2],p=p_mod[2],p_disp=p_disp[2],p_creation=p_creation[2],m_disp=m_disp[2],m_creation=m_creation[2])
      else{
        Ag_new<-graph_alteration(Ag,m=m[1],p=p_mod[1],p_disp=p_disp[1],p_creation=p_creation[1],m_disp=m_disp[1],m_creation=m_creation[1])
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
  return(list(dist=dist, data=graph_seq))
}




##########################################################################################################
############################         Change point detection in ER evolution          #####################
##########################################################################################################

### ------------------------------------------------------------------------------------------------------
############################    2-phase dynamic regime  in   evolving ER graph  ##########################
### ------------------------------------------------------------------------------------------------------
test_change_point_detection_ER_two_phases<-function(N, T,time_change_point=floor(T/2),p0=0.3,params_evolution=c(0.3,0,0),params_evolution_new=c(0.6,0,0),prop=0.05,prop_new=0.05, plot=TRUE){
    ##  Description
    ##  -------------
    ##  Compares how different distances detect a dynamical regime change. At time t=time_change_point, the evolution
    ##  process changes.
    ##  INPUT:
    ##  =============================================================
    ##  N                   :   number of nodes
    ##  T                   :   length of the process horizon
    ##  time_change_point   :   time at which the change in dynamics occurs
    ##  p0                  :   probability of edge creation in ER graph (initial graph)
    ##  prop                :   proportion of the edges that are potentially modified (first time regime)
    ##  params_evolution    :   parameters for the first evolution phase:
    ##                          params_evolution[1]=p=probability that a candidate ''modifiable'' edge is reassigned somewhere else,
    ##                          params_evolution[2]=p_disp,
    ##                          params_evolution[3]=p_creation
    ##  params_evolution_new:   parameters for the second evolution phase
    ##  prop_new            :   proportion of the edges that are potentially modified (first time regime)
    ##  p_new               :   probability that a candidate ''modifiable'' edge is reassigned somewhhere else (2nd time regime)
    ##
    ##  OUTPUT
    ##  =============================================================
    ## distance             :    matrix of pairwise distances between consecutive graphs
    
    A<-generate_random_adjacency(N,p0,sym=TRUE, plot=FALSE)
    distances<-matrix(0,4,(T-1))
    for ( t in 1:time_change_point){
        A_new<-random_alteration_adjacency(A,prop,params_evolution[1],params_evolution[2],params_evolution[3])
        distances[1,t]<-abs(get_number_spanning_trees2(A_new)-get_number_spanning_trees2(A))/(get_number_spanning_trees2(A_new)+get_number_spanning_trees2(A))
        distances[5:7,t]<-netdist(A, A_new,d = "HIM")
        #print(distances[,t])
        A<-A_new
    }
    for ( t in (time_change_point+1):(T-1)){
        A_new<-random_alteration_adjacency(A,prop_new,params_evolution_new[1],params_evolution_new[2],params_evolution_new[3])
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

### ------------------------------------------------------------------------------------------------------
###############      ER graph: evolution with three phases        ########################################
### ------------------------------------------------------------------------------------------------------

test_change_point_detection_ER_three_phases<-function(N,p,prop=c(0.05,0.2),p_disp=c(0,0,0),p_creation=c(0,0,0),T=21,verbose=TRUE,save_graph_seq=T,name_file_ext=""){
    ##  Description
    ##  -------------
    ##  Function for generating an ER random graph, and computing distance between consecutive graphs
    ##  as this graph evolves from time point to time point (using different metrics: ST, Hamming, IM..).
    ##  The evolution mechanism is again the same: we chose prop % of the edges: these selected edges are
    ##  potentially modifiable edges. Each candidate edge is then rewired with probability p[1], and stays in place
    ##  with probability 1-p[1].
    ##  At time t=10, the dynamics change (that is, the probability that a candidate modifiable edge
    ##  is rewired changes from p[1] to p[2]. At time t=16, the probability goes back to p[1], but the proportion
    ##  of modifiable edges oncreases from prop[1] to prop[2];
    ##
    ##  INPUT:
    ##  =============================================================
    ##  N               :   Number of nodes in the graphs (int)
    ##  opts            :   what type of random graph (1: PA, 2: Island, 3: Dot Product, 4: SBM). Default=1
    ##  p               :   vector of length 3 or float. If vector, the first entry corresponds to the original
    ##                      ER graph parameter, and the two next entries to the two distinct rewiring probabilities.
    ##                      If float, then a 3-d vector is generated as [p,p,0.6]
    ##  prop            :   vector of length 2 or float. If vector, the two entries to the two distinct proportions
    ##                      of modifiable edges.If float, then a 2-d vector is generated as [prop,0.2]
    ##  p_disp          :   3-d vector for the probabilities of disparition in the three phases
    ##  p_creation      :   3-d vector for the probabilities of creation in the three phases
    ##  T               :   number of different graphs that are generated. (int>1)
    ##  verbose         :   boolean. Should intermediary summary messages be printed? (useful for
    ##                      debugging. Default=TRUE
    ##  save_graph_seq  :   boolean.Should the graphs be saved? Default=TRUE
    ##  name_file_ext   :   string (name of the file were the graphs should be saved if save_graph_seq
    ##                      is TRUE. Default=""
    ##
    ##  OUTPUT
    ##  =============================================================
    ##  list            :   list with named arguments dist_changepnt (distance between consecutive  graphs)
    ##                      and data (a matrix where each column corresponds to a flattened graph)
    if (length(p)==1){
        prob=p
        p=matrix(0,3)
        p[1]=prob
        p[2]=prob
        p[3]=0.6
    }
    if (length(prop)==1){
        prop0=prop
        prop=matrix(0,2)
        prop[1]=prop0
        prop[2]=0.2
    }
    
    print("further tests: comparison for random innovation, with regime change at t=10 and 16:")
    print(paste("At t=10, the probability of a new link increases from ", p[2] , "to ", p[3], "(with identical changed proportion prop=",100*prop[1],"% of edges changed)"))
    print(paste("At t=16, the probability of a new link goes back to ", p[2]  , "but",100*prop[2],"% of the edges are unplugged and replugged"))
    #### Compare against totally random networks
    A<-vector("list",T)
    A_new<-vector("list",T-1)
    A[[1]]<-generate_random_adjacency(N,p[1], TRUE)
    plot(graph_from_adjacency_matrix(A[[1]],mode = "undirected"), main="initial graph")
    graph_seq<-matrix(0,T,N^2)
    graph_seq[1,]<-as.vector(t(A[[1]]))
    dist_changepnt<-matrix(0, 16,T-1)
    prop=0.05
    T=21
    for (t in 1:(T-1)){
        if (verbose==TRUE) print(paste("time: ",t))
        if(t<10){
            A_new[[t]]<-random_alteration_adjacency(A[[1]],prop[1],p[2], p_disp[1], p_creation[1])
        }
        else{
            if(t<16){
                A_new[[t]]<-random_alteration_adjacency(A[[t]],prop[1],p[3], p_disp[2], p_creation[2])
            }
            else{
                A_new[[t]]<-random_alteration_adjacency(A[[t]],prop[2],p[2], p_disp[3], p_creation[3])
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
    
    
    
    
    plot(1:(T-1),dist_changepnt[2,],ylim=c(0,1), type="l",col=2,xlab="time", ylab="distance between consecutive graphs", main="Evolution of the distance between consecutive graphs (with change point)")
    points(1:(T-1),dist_changepnt[1,], type="l",col=1,xlab="time", ylab="distance between consecutive graphs")
    points(1:(T-1),dist_changepnt[3,], type="l",col=3,xlab="time", ylab="distance between consecutive graphs")
    points(1:(T-1),dist_changepnt[4,], type="l",col=4,xlab="time", ylab="distance between consecutive graphs")
    points(1:(T-1),dist_changepnt[5,]/max(dist_changepnt[5,]), type="l",col=5,xlab="time", ylab="distance between consecutive graphs")
    legend(1,1,c( "ST-based dist",  "H","IM","HIM",paste("Poly x",max(dist_changepnt[5,]))),col=1:5,lwd=2,cex=0.7)
    
    if (save_graph_seq){
        #dist=data.frame(dist_changepnt,row.names = c("ST-based dist",  "H","IM","HIM","alpha=0.5,order=5"))
        write.table(dist_changepnt,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_comp_changepnt_",N,"_",m,"_",name_file_ext,"distances.txt",sep=""),row.names = TRUE,col.names = FALSE)
        write.table(graph_seq,sep=" ",file = paste("/Users/cdonnat/Dropbox/Distances/write_up/plot/plot_analysis/test_random_comp_changepnt_",N,"_",p,"_",name_file_ext,".txt",sep=""),row.names = FALSE,col.names = FALSE)
    }
    return(list(dist_changepnt=dist_changepnt,data=graph_seq))
}
### ------------------------------------------------------------------------------------------------------












