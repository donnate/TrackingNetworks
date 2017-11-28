#### Run sequences of tests for clustering purposes:
#dir="/Users/cdonnat/Dropbox/Food_network/code"
#setwd(dir) 
source("~/Dropbox/Food_network/code/ADMM.R")
source("~/Dropbox/Food_network/code/tools.R")
source("~/Dropbox/Food_network/code/SC.R")
library(glmnet)
library(gplots)
library(stats)

test_ADMM2<-function(N,K){
  B<-matrix(rnorm(N^2,0,1),N^2,1)
  for (k in 1:(K-1)){
    B_k<-matrix(rnorm(N^2,0,1),N^2,1)
    B<-cbind(B,B_k)
  }
  A=3*B[,1]+5*B[,2]+0.3*B[,3]+2*B[,K]
  return(ADMM(A,B,rho=1,Z_alpha=100,Time=10,dir=1))
  
}






#### Comparison performance L1 vs L2 penalty
### L1 seems to achieve perfect results 

default_ADMM_args<-list(rho=10,eta=0.1,Z_alpha=100,Z_penalty=1,Time=50,dir=1,delta_1=10^(-4),delta_2=10^(-4))

test_recovery_performance<-function(N,K,M,args_ADMM=default_ADMM_args){
  coeff<-matrix(0,M,K)
  recovered_coeffsL1<-matrix(0,M,K)
  recovered_coeffsL2<-matrix(0,M,K)
  recovered_coeffsL3<-matrix(0,M,K)
  sparsity<-matrix(0,M,4)
  optimal_index<-matrix(0,M,2)
  timeit<-matrix(0, M,3)
  nbit<-matrix(0, M,2)
  w=1
  flag_eta=K/3
  weird_sample<-c()
  for( b in 1:M){
    ### pick a random nb of coeffs
    nb_coeff=sample(1:floor(K/2),1) 
    sparsity[b,1]=nb_coeff/K
    coeff[b,]<-sapply(1:K, FUN=function(i){
      if(i %in% spots){
        return(runif(1,-10,10))}
      else{
        return(0)
      }
    })
    B<-matrix(rnorm(N^2,0,1),N^2,1)
    A<-coeff[b,1]*B
    for (k in 1:(K-1)){
      B_k<-matrix(rnorm(N^2,0,1),N^2,1)
      B<-cbind(B,B_k)
      A<-A+coeff[b,k+1]*B_k
    }
    A<-A+matrix(rnorm(N^2,0,1),N^2,1) ## plus noise
    #print(paste("finished generating test sample nb ",b))
    L1_test=ADMM(A,B,args_ADMM,plot=FALSE)
    L2_test=ADMM(A,B,args_ADMM,plot=FALSE)
    start=Sys.time()
    L3_test<-cv.glmnet(B,A,alpha=1)
    best.lambda <- which(L3_test$lambda==L3_test$lambda.min)
    timeit[b,3]=Sys.time()-start
    optimal_index[b,1]=L1_test$suggested.indice
    optimal_index[b,2]=L2_test$suggested.indice
    timeit[b,1]=L1_test$time_taken
    timeit[b,2]=L2_test$time_taken
    nbit[b,1]=L1_test$nb_it
    nbit[b,2]=L2_test$nb_it
    recovered_coeffsL1[b,]=L1_test$lambda[,optimal_index[b,1]]
    sparsity[b,2]=length(which(recovered_coeffsL1[b,]!=0))/K
    recovered_coeffsL2[b,]=L2_test$lambda[,optimal_index[b,2]]
    sparsity[b,3]=length(which(recovered_coeffsL2[b,]!=0))/K
    recovered_coeffsL2[b,]=L2_test$lambda[,optimal_index[b,2]]
    sparsity[b,4]=L3_test$glmnet.fit$df[best.lambda]/K
    recovered_coeffsL3[b,]=L3_test$glmnet.fit$beta[,best.lambda]
    ### Store the weird samples that have high error for further investigation
    if (sum((recovered_coeffsL2[b,]-coeff[b,])^2)>flag_eta || sum((recovered_coeffsL1[b,]-coeff[b,])^2)>flag_eta  || sum((recovered_coeffsL3[b,]-coeff[b,])^2)>flag_eta ){
      weird_sample[[w]]=list(A=A,B=B,b=b)
      w<-w+1
    }
  }
  return(list(coeff=coeff,recovered_coeffsL1=recovered_coeffsL1,recovered_coeffsL2=recovered_coeffsL2,recovered_coeffsL3=recovered_coeffsL3,sparsity=sparsity,weird_sample=weird_sample, timeit=timeit))
}

#test_results=test_recovery_performance(100,30,100,Time=50)

analyze_results<-function(test_results){
  M=nrow(test_results$timeit)
  K=ncol(test_results$coeff)
  par(mfrow=c(2,2))
  tlim=max(test_results$timeit)+1.5
  plot(1:M,test_results$timeit[,1],col="red", xlab="trial",ylab="time", main="comparison running time",ylim=c(0,tlim),pch=19)
  points(1:M,test_results$timeit[,2],col="black",pch=19)
  points(1:M,test_results$timeit[,3],col="green",pch=19)
  legend(1,tlim,c("ADMM with L1-norm for Z penalty","ADMM with L2-norm for Z penalty","CV Lasso"),col=c("red","black","green"),cex=0.7,pch=19)
  diffL1<-(test_results$recovered_coeffsL1-test_results$coeff)^2
  diffL2<-(test_results$recovered_coeffsL2-test_results$coeff)^2
  diffL3<-(test_results$recovered_coeffsL3-test_results$coeff)^2
  perfL1<-apply(diffL1,1,sum)
  perfL2<-apply(diffL2,1,sum)
  perfL3<-apply(diffL3,1,sum)
  lim=max(c(perfL1,perfL2,perfL3))*1.05
  print(paste("max difference:",lim, " by  method number ", which.max(c(max(perfL1),max(perfL2),max(perfL3)))))
  print(c("max difference:",c(max(perfL1),max(perfL2),max(perfL3))))
  print(c("average difference:",c(mean(perfL1),mean(perfL2),mean(perfL3))))
  print(c("sd of difference:",c(sd(perfL1),sd(perfL2),sd(perfL3))))
  whos_best=sapply(1:length(perfL1),FUN=function(i){return(which.min(c(perfL1[i],perfL2[i],perfL3[i])))})
  hist(whos_best)
  plot(1:M,perfL1,col="red", xlab="trial",ylab="|| hat_lambda-lambda||^2", ylim=c(0.0,lim),main="Performance",pch=19 )
  points(1:M,perfL2,col="black",pch=19)
  points(1:M,perfL3,col="green",pch=19)
  legend(3,lim,c("ADMM with L1-norm for Z penalty","ADMM with L2-norm for Z penalty","CV Lasso"),col=c("red","black","green"),cex=0.7,pch=19)
  which_is_better=perfL1-perfL2
  #points(which(which_is_better>0),matrix(0,length(which(which_is_better>0)),1),col="yellow",pch=4)
  print(paste("L2 penalty outperforms L1 in ",  100*length(which(which_is_better>0))/M, "% of cases"))
  plot(1:M,test_results$sparsity[,2],col="red", xlab="trial",ylab="norm of || hat_lambda-lambda||^2", ylim=c(0,1),main=" Reg. Coeff. recovery",pch=19 )
  points(1:M,test_results$sparsity[,3],col="black",pch=19)
  points(1:M,test_results$sparsity[,4],col="green",pch=19)
  points(1:M,test_results$sparsity[,1],col="yellow",pch=14)
  legend(3,1,c("ADMM with L1-norm for Z penalty","ADMM with L2-norm for Z penalty","CV Lasso","ground truth"),col=c("red","black","green","yellow"),cex=0.7,pch=c(19,19,19,14))
  sparsity_summary=apply(test_results$sparsity,2,mean)
  print(c("Average true and recovered sparsity ",  sparsity_summary))
  
  
}
###test clusters
library(LCA)
generate_clusters<-function(N,K, nb_clusters,args_ADMM=default_args_ADMM, even_size=TRUE,plot_clst=FALSE){
  if (even_size==TRUE){
    size1<-rep(floor(K/nb_clusters),nb_clusters-1)
    size<-c(size1,K-sum(size1))
  }
  else{
    size<-rdirichlet(nb_clusters, rep(1/nb_clusters,nb_clusters))
  }
  clusters<-list()
  clusters_m<-matrix(0,K,N^2)
  coefficients<-matrix(0,nb_clusters,K)
  labels<-matrix(0,K,1)
  sparsity<-matrix(0,nb_clusters,1)
  a=1
  b=size[1]
  for (k in 1:nb_clusters){
    range_clust=a:b
    labels[a:b]=k
    nb_coeff=sample(1:size[k]/2,1)
    spots<-sample(a:b,nb_coeff, replace = FALSE)
    sparsity[k]=nb_coeff/K
    coeff<-runif(nb_coeff,-10,10)
    coefficients[k,spots]<-coeff
    coeff_shrunk<-coefficients[k,a:b]
    B<-matrix(rnorm(N^2,0,1)*size[k] ,N^2,size[k])
    clusters[[k]]<-B
    for (m in 1:size[k]){
      clusters_m[a+m-1,]<-B%*%coeff_shrunk+matrix(rnorm(N^2,0,1),N^2,1)
    }
    if (k<K){
      a=b+1
      b=b+size[k+1]
    }
    
  }
  if(plot_clst==TRUE) heatmap.2(clusters_m,main="clusters")
  return(list(X=clusters_m,labels=labels,sparsity=sparsity,clusters=clusters))
  
  
}




get_convex_format<-function(target,data){
  indices=setdiff(1:nrow(data),target)
  return(list(A=data[target,], B=data[indices,],regressors=indices))
}

test_cluster<-function(N,K,nb_clusters,method="Lasso",even_size=TRUE,Z_penalty=1,args_ADMM=default_ADMM_args,plots=FALSE){
  Sim=matrix(0,K,K)
  data<-generate_clusters(N,K,nb_clusters,even_size=TRUE,plot_clst=plots)
  for (k in 1:nrow(data$X)){
    pb<-get_convex_format(k,data$X)
    if (method=="Lasso"){
      reg<-cv.glmnet(t(pb$B),pb$A,alpha=1)
      best.lambda <- which(reg$lambda==reg$lambda.min)
      Sim[k,pb$regressors]<-reg$glmnet.fit$beta[,best.lambda]
    }
    else{
      if(is.null(args_ADMM)){
        args_ADMM=default_ADMM_args
      }else{
        for (narg in names(default_ADMM_args)){
          if(is.null(args_ADMM[[narg]])){
            args_ADMM[[narg]]=default_ADMM_args[[narg]]
          }
        }
      }
      reg<-ADMM(pb$A,t(pb$B),args_ADMM,plot=TRUE)
      opt <- reg$suggested.indice
      Sim[k,pb$regressors]<-reg$lambda[,opt]
    }
  }
  return(Sim)
}

test_cluster_methods<-function(N,K,nb_clusters,argsADMM=default_ADMM_args,even_size=TRUE,Z_penalty=1,plots=FALSE){
  Sim1=matrix(0,K,K)
  Sim2=matrix(0,K,K)
  Sim3=matrix(0,K,K)
  data<-generate_clusters(N,K,nb_clusters,even_size=TRUE,plot_clst=plots)
  argsADMM1<-argsADMM
  argsADMM1[["Z_penalty"]]=1
  argsADMM1[["rho"]]=10
  argsADMM2<-argsADMM
  argsADMM2[["Z_penalty"]]=2
  argsADMM2[["rho"]]=10
  
  for (k in 1:nrow(data$X)){
    pb<-get_convex_format(k,data$X)
    reg3<-cv.glmnet(t(pb$B),pb$A,alpha=1)
    best.lambda <- which(reg3$lambda==reg3$lambda.min)
    Sim3[k,pb$regressors]<-reg3$glmnet.fit$beta[,best.lambda]
    
    reg1<-ADMM(pb$A,t(pb$B),argsADMM1,plot=FALSE)
    opt <- reg1$suggested.indice
    Sim1[k,pb$regressors]<-reg1$lambda[,opt]
    reg2<-ADMM(pb$A,t(pb$B),argsADMM2,plot=FALSE)
    opt2 <- reg2$suggested.indice
    Sim2[k,pb$regressors]<-reg2$lambda[,opt2]
    
  }
  #print("okidok")
  return(list(Sim1=Sim1,Sim2=Sim2,Sim3=Sim3,data=data))
}

post_processing<-function(Sim,thres=0,nb_clusters,laplacian="normalized",true_labels,go_for_plotting=TRUE){
  Adj<-hard_threshold(Sim,thres)
  ### Run spectral clustering 
  SC<-Spectral_Clustering(Adj,nb_clusters,laplacian=laplacian, plot=F)
  km <- kmeans(SC, centers=nb_clusters, nstart=5)
  if(go_for_plotting==TRUE){
    par(mfrow=c(1,2))
    plot(SC, col=true_labels, pch=20,main="Projections colored by true label")
    plot(SC, col=km$cluster, pch=20,main="Projections colored by kmeans label")
  }
  
  
  return(list(Z=SC,clusters=km$cluster))
}
  
  
  
go_clusters<-function(N,K,nb_clusters,plots=FALSE){
  test_clst=test_cluster_methods(N,K,nb_clusters)
  par(mfrow=c(2,2))
  heatmap.2(test_clst$Sim1+t(test_clst$Sim1),Rowv=FALSE,Colv=FALSE,trace="none",dendrogram = "none",main="ADMM L1")
  heatmap.2(test_clst$Sim2+t(test_clst$Sim2),Rowv=FALSE,Colv=FALSE,trace="none",dendrogram = "none",main="ADMM L2")
  heatmap.2(test_clst$Sim3+t(test_clst$Sim3),Rowv=FALSE,Colv=FALSE,trace="none",dendrogram = "none",main="ADMM Lasso")
  Sim1<-test_clst$Sim1+t(test_clst$Sim1)
  Sim2<-test_clst$Sim2+t(test_clst$Sim2)
  Sim3<-test_clst$Sim3+t(test_clst$Sim3)
  spar1=1-sum(sapply(Sim1,FUN=function(x){
    delta(x,0)}))/K^2
  spar2=1-sum(sapply(Sim2,FUN=function(x){
    delta(x,0)}))/K^2
  spar3=1-sum(sapply(Sim3,FUN=function(x){
    delta(x,0)}))/K^2
  print("Sparsity of the adjacency matrices:")
  test_clst$data$sparsity
  par(mfrow=c(1,3))
  hist(Sim1)
  hist(Sim2)
  hist(Sim2)
  p1<-post_processing(Sim1,nb_clusters=nb_clusters,true_labels=test_clst$data$labels)
  #cluster.stats(d, fit1$cluster, test_clst$data$labels)
  
  #p2<-post_processing(Sim2,nb_clusters)
  p3<-post_processing(Sim3,nb_clusters=nb_clusters,true_labels=test_clst$data$labels)
  labels=test_clst$data$labels
  #cluster.stats(p1$clusters,labels)
  #pur1<-purity(p1$clusters,labels)
  #cluster.stats(Sim1, p1$clusters,labels)
  ### Run clustering on the inferred similarity matrix
  
}
####  Test more simple: standard LASSO form

library(fpc)

#cluster.stats(Sim1, p1$clusters,labels)


##### test ability to see smooth evolution over time

generate_temporal_evolution<-function(N,K,lambda=1){
  data<-matrix(0,K,N^2)
  data[1,]<-matrix(rnorm(N^2,3,1),N^2,1)  ### define the first signal  ##wth a stronger signal than the rest?
  for (k in 2:K){
    data[k,]<-lambda*data[k-1,]+matrix(rnorm(N^2,3,1),N^2,1)  ### innovation in that case is simply random noise ## so it is in particular isotropic
  }
  return(data)
}



get_Similarity<-function(data,method="Lasso",args_ADMM=NULL,plot_green_light=TRUE,verbose=FALSE){
  
  Sim=matrix(0, nrow(data),nrow(data))
  
  for (k in 1:nrow(data)){
    if (verbose==TRUE) print(k/nrow(data))
    pb<-get_convex_format(k,data)
    if (method=="Lasso"){
      if(sum(abs(pb$A))>0){
        reg<-cv.glmnet(t(pb$B),pb$A,alpha=1)
        best.lambda <- which(reg$lambda==reg$lambda.min)
        Sim[k,pb$regressors]<-reg$glmnet.fit$beta[,best.lambda]
      }
    }
    else{
      if(is.null(args_ADMM)){
        args_ADMM=default_ADMM_args
      }else{
        for (narg in names(default_ADMM_args)){
          if(is.null(args_ADMM[[narg]])){
            args_ADMM[[narg]]=default_ADMM_args[[narg]]
          }
        }
      }
      reg<-ADMM(pb$A,t(pb$B),args_ADMM,plot=FALSE)
      opt <- reg$suggested.indice
      Sim[k,pb$regressors]<-reg$lambda[,opt]
    }
  }
  return(Sim)
  
}

test_evolution<-function(N,K,data=NULL,method="Lasso",args_ADMM=NULL,plot_green_light=TRUE){
  Sim=matrix(0,K,K)
  
  if (is.null(data)){
    print("Random Weighted, directed graph is generated and evolves as G_{t+1}=G_{t}+Z_{t+1}")
    data<-generate_temporal_evolution(N,K)
  }
  else{
    print("Attempts to recover the evolution of the graph sequence passed as input data")
  }
  
  for (k in 1:nrow(data)){
    pb<-get_convex_format(k,data)
    if (method=="Lasso"){
      reg<-cv.glmnet(t(pb$B),pb$A,alpha=1)
      best.lambda <- which(reg$lambda==reg$lambda.min)
      Sim[k,pb$regressors]<-reg$glmnet.fit$beta[,best.lambda]
    }
    else{
      if(is.null(args_ADMM)){
        args_ADMM=default_ADMM_args
      }else{
        for (narg in names(default_ADMM_args)){
          if(is.null(args_ADMM[[narg]])){
            args_ADMM[[narg]]=default_ADMM_args[[narg]]
          }
        }
      }
      reg<-ADMM(pb$A,t(pb$B),args_ADMM,plot=FALSE)
      opt <- reg$suggested.indice
      Sim[k,pb$regressors]<-reg$lambda[,opt]
    }
  }
  Sim_sym=0.5*(Sim +t(Sim))  ## make it a symmetric matrix
  norm_sym<-apply(Sim^2,2,sum)
  Sim_trans=sapply(1:K,FUN=function(i){
    return(Sim[i,]^2/norm_sym[i])
  })
  ### Sym transform: best relects how close the matrices are from one another
  if (plot_green_light==TRUE){
    
    heatmap.2(Sim_sym,Rowv=FALSE, Colv=FALSE, dendrogram="none",trace="none",main=paste("Symmetrized similarity using ",method))
    #heatmap.2(Sim_trans,Rowv=FALSE, Colv=FALSE, dendrogram="none" ,trace="none",main=paste("Weights of each coefficient for ",method))
  }
  return(list(Sim=Sim, data=data,Sim_trans=Sim_trans))
}



comparison_performance<-function(N,K){
  data=generate_temporal_evolution(N,K)
  argsADMM1<-default_ADMM_args
  argsADMM1[["Z_penalty"]]=1
  argsADMM1[["rho"]]=5
  argsADMM2<-default_ADMM_args
  argsADMM2[["Z_penalty"]]=2
  testL1<-test_evolution(N,K,data=data, method="ADMM1",args_ADMM=argsADMM1)
  testL2<-test_evolution(N,K,data=data, method="ADMM2",args_ADMM=argsADMM2)
  testL3<-test_evolution(N,K,data=data, method="Lasso")
  
}







#### Questions quise posent: on a fait des tests, amis seulement avec des donnees normales gaussiennes... La question se pose: est-ce que ca marche pour autre chose?
#### Especially in terms of interaction graphs, where we mostly handle count data... I guess the approximation is justified if we take count data. normalized b volume and centered...
#### But it does require a little bit of tweaking around the data

get_Similarity_par<-function(data,method="Lasso",args_ADMM=NULL,plot_green_light=TRUE){
  
  Sim=matrix(0, nrow(data),nrow(data))

  Sim<-sapply( 1:nrow(data),FUN =function(k){
    print(k)
    pb<-get_convex_format(k,data)
    reg<-cv.glmnet(t(pb$B),pb$A,alpha=1)
    best.lambda <- which(reg$lambda==reg$lambda.min)
    sim[pb$regressors]<<-reg$glmnet.fit$beta[,best.lambda]
    return(sim)
  })
  return(Sim)
  
}

investigate_weird_results<-function(test_coeff){
  investigate=test_coeff$weird_sample
  for (k in 1:length(investigate)){
    b=investigate[[k]]$b
    real_coeff=test_coeff$coeff[b,]
    rec_coeff1=test_coeff$recovered_coeffsL1[b,]
    M=100
    diff1<-matrix(0,1,M)
    for (m in 1:M){
      L1_test=ADMM(investigate[[k]]$A,investigate[[k]]$B,args_ADMM,plot=FALSE)
      diff1[m]=sum((L1_test$lambda[,L1_test$suggested.indice]-real_coeff)^2)
      print(diff1[m])
    }
    L1_test=ADMM(investigate[[k]]$A,investigate[[k]]$B,args_ADMM,plot=FALSE)
    L1_test=ADMM(investigate[[k]]$A,investigate[[k]]$B,args_ADMM,plot=TRUE)
    diff1=sum((L1_test$lambda[,L1_test$suggested.indice]-real_coeff)^2); print(diff1)
    print(L1_test$suggested.indice)
    diff1=sum((L1_test$lambda[,L1_test$suggested.indice]-real_coeff)^2)
    rec_coeff2=test_coeff$recovered_coeffsL2[b,]
    rec_coeff3=test_coeff$recovered_coeffsL3[b,]
    diff2=(rec_coeff2-real_coeff)^2
    diff3=(rec_coeff3-real_coeff)^2
  }
  
}

