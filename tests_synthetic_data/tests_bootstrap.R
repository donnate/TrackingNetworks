source("./tests_synthetic_data/dynamics.R")
source("./distances/spanning_trees.R")
source("./distances/distances.R")

#### bootstrap tests for understanding the distance between two graphs of the same type
#### These functions were designed so as to be then used as part of shiny applictions
#### which might explain some of their layout



### --------------------------------------------------------------------
###########        First simple bootstrap test       ###################
### --------------------------------------------------------------------
test_bootstrap<-function(N,args,opts=1,B,param2test=data.frame(order_max=c(5,5,3,3,3,1),alpha=c(0.5,1.2,0.5,0.9,1,1))){
  boot_samples<-sapply(1:B, FUN=function(b){
    if (opts==0){
      p=args$p
      A<-generate_random_adjacency(N,p, TRUE)
      A2<-generate_random_adjacency(N,p, TRUE)
      
    }
    else{
      A<-generate_realistic_adjacency(N,args_l=args,opts=opts, verbose=FALSE)
      A<-get.adjacency(A)
      A2<-generate_realistic_adjacency(N,args_l=args,opts=opts, verbose=FALSE)
      A2<-get.adjacency(A2)
    }
    
    dist<-matrix(0,nrow(param2test),1)
    for ( i in 1:nrow(param2test)){
      dist[i]=poly_distance(A, A2,order_max=param2test$order_max[i],alpha=param2test$alpha[i])
    }

    return(dist)
  })
}
### --------------------------------------------------------------------









###############################################################################################
###########   Functions for bootstrap generation ( one for each different distance) ###########
###############################################################################################

### --------------------------------------------------------------------
###########  bootstrap test for polynomial distances    ################
### --------------------------------------------------------------------
test_bootstrap_shiny<-function(N,alpha,order_max,args,opts=1,B){
  boot_samples<-sapply(1:B, FUN=function(b){
    if (opts==0){
      p=args$p
      A<-generate_random_adjacency(N,p, TRUE)
      A2<-generate_random_adjacency(N,p, TRUE)
      
    }
    else{
      A<-generate_realistic_adjacency(N,args_l=args,opts=opts,verbose=FALSE)
      A<-as.matrix(get.adjacency(A))
      A2<-generate_realistic_adjacency(N,args_l=args,opts=opts, verbose=FALSE)
      A2<-as.matrix(get.adjacency(A2))
    }
    dist=poly_distance(A, A2,order_max,alpha=alpha)
    return(dist)
  })
}
### --------------------------------------------------------------------


### --------------------------------------------------------------------
####### bootstrap test for polynomial distances  (2nd version)  ########
### --------------------------------------------------------------------
test_bootstrap_shiny2<-function(N,alpha,order_max,args,opts=1,opts2=2,B){
  boot_samples<-sapply(1:B, FUN=function(b){
    if (opts==0){
      p=args$p
      A<-generate_random_adjacency(N,p, TRUE)
    }
    else{
      A<-generate_realistic_adjacency(N,args_l=args,opts=opts, verbose=FALSE)
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
### --------------------------------------------------------------------



### --------------------------------------------------------------------
###########  bootstrap test for Spanning Tree distances    #############
### --------------------------------------------------------------------
test_bootstrap_shiny_ST<-function(N,args,opts=1,B,normalize_ST=FALSE){
  boot_samples<-sapply(1:B, FUN=function(b){
    print(b)
    if (opts==0){
      p=args$p
      A<-generate_random_adjacency(N,p, TRUE)
      A2<-generate_random_adjacency(N,p, TRUE)
      
    }
    else{
      A<-generate_realistic_adjacency(N,args_l=args,opts=opts, verbose=FALSE)
      A<-as.matrix(get.adjacency(A))
      A2<-generate_realistic_adjacency(N,args_l=args,opts=opts, verbose=FALSE)
      A2<-as.matrix(get.adjacency(A2))
      
    }
    stree<-get_number_spanning_trees(A)
    stree2<-get_number_spanning_trees(A2)
  
    dist=ifelse(normalize_ST,abs(stree-stree2)/abs(0.5*(stree+stree2)),abs(stree-stree2))
    return(dist)
  })
}
### --------------------------------------------------------------------


### --------------------------------------------------------------------
###########  bootstrap test for HIM distances    ################
### --------------------------------------------------------------------
test_bootstrap_shiny_HIM<-function(N,args,opts=1,B){
  boot_samples<-sapply(1:B, FUN=function(b){
    if (opts==0){
      p=args$p
      A<-generate_random_adjacency(N,p, TRUE)
      A2<-generate_random_adjacency(N,p, TRUE)
      
    }
    else{
      A<-generate_realistic_adjacency(N,args_l=args,opts=opts, verbose=FALSE)
      A<-as.matrix(get.adjacency(A))
      A2<-generate_realistic_adjacency(N,args_l=args,opts=opts, verbose=FALSE)
      A2<-as.matrix(get.adjacency(A2))
      
    }
    dist=netdist(A, A2,d = "HIM")
    return(dist)
  })
}
### --------------------------------------------------------------------






###############################################################################################
#########################         Test function         #######################################
###############################################################################################

### --------------------------------------------------------------------
tb<-function(){
  args<-list(power=0.9,islands.n=3,islands.size=10,islands.pin=0.3,n.inter=3,K=3,block.sizes=c(10,10,10),pm=cbind( c(.4,0.1, .001), c(.1,0.2, .01),c(.001,0.01, .5)))
  N=30
  opts=1
  B=1000
  test_bootstrap1<-test_bootstrap(N,args,opts,B)
}
### --------------------------------------------------------------------







###############################################################################################
#################            Comparison of distances function         #########################
###############################################################################################




### --------------------------------------------------------------------
compare_several_distances<-function(graph_seq,args=list(order_max=3, alpha=1)){
    ##  Description
    ##  -------------
    ##  Compares the different distances by compting the pairwise distance matrix between the graphs in
    ##  graph_seq, and plotting the Minimum Spanning Trees (for each type of distance: Hamming, IM, HIM,ST,
    ##  polynomial-- the latter with parameters as specified by the users) obtained from this
    ##  pairwise distance matrix.Ideally, we expect that if graph_seq represents a sequence of
    ##  consecutive graphs, the MST will link consecutive graphs.
    ##
    ##  INPUT:
    ##  =============================================================
    ##  graph_seq        :   a sequence of graphs (list of size K x N x N where K is the number of graphs
    ##                       and N is the number of nodes.)
    ##  args             :   list (with names arguments) for the polynomial distances
    ##
    ##  OUTPUT
    ##  =============================================================
    ## test_dist        :    matrix of pairwise distances between consecutive graphs
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
  go_plot=TRUE
  if (go_plot==TRUE){
    plot(mst(g1,weights = edge_attr(g1, "weight") ), vertex.size=5,main="Spanning Tree from Poly")
    plot(mst(g2,weights = edge_attr(g2, "weight") ), vertex.size=5,main="Spanning Tree from Hamming")
    plot(mst(g3,weights = edge_attr(g3, "weight") ), vertex.size=5,main="Spanning Tree from IM")
    plot(mst(g4,weights = edge_attr(g4, "weight") ), vertex.size=5,main="Spanning Tree from HIM")
    plot(mst(g5,weights = edge_attr(g5, "weight") ), vertex.size=5,main="Spanning Tree from ST Sim")
  }
  
  
  return(list(dist1=dist1,dist2=dist2,dist3=dist3,dist4=dist4,dist5=dist5,mst1=mst1,mst2=mst2,mst3=mst3,mst4=mst4,mst5=mst5))
}
### --------------------------------------------------------------------














##########################################################################################################
############################             Distributions               #####################################
##########################################################################################################



### ------------------------------------------------------------------------------------------------------
############################################      distributions        ###################################
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




### ------------------------------------------------------------------------------------------------------
############################################      test functions        ##################################
### ------------------------------------------------------------------------------------------------------
test_get_dist<-function(){
    ##  Description
    ##  -------------
    ##  Compares histograms of distances
    N=30
    p_v=c(0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7)
    prop_v=c(0.05,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
    B=1000
    get_res=get_distribution(N,p_v,prop_v,B)
    save(get_res,p_v,prop_v,N,file="./saved_data/get_res.RData")
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


