source("./tests_synthetic_data/test_functions.R")
source("./spanning_trees.R")
source("./distances.R")

#### bootstrap tests for understanding the distance between two graphs of the same type
#### These functions were designed so as to be then used as part of shiny applictions
#### which might explain some of their layout



### --------------------------------------------------------------------
###########        First simple bootstrap test       ###################
### --------------------------------------------------------------------
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
      A<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
      A<-as.matrix(get.adjacency(A))
      A2<-generate_realistic_adjacency(N,args,opts, verbose=FALSE)
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
### --------------------------------------------------------------------



### --------------------------------------------------------------------
###########  bootstrap test for Spanning Tree distances    #############
### --------------------------------------------------------------------
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
### --------------------------------------------------------------------


### --------------------------------------------------------------------
###########  bootstrap test for HIM distances    ################
### --------------------------------------------------------------------
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
### --------------------------------------------------------------------
