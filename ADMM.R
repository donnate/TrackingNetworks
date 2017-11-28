####  Analyse des proprietes du graphe:
#library(PBSmodelling)
library(tictoc)



ADMM<-function(A,B,args_ADMM=default_ADMM_args,plot=TRUE){
  # A is a  n \times 1 vector
  # B is a $K \times n$ matrix
  # lambda should be a K dimensional vector
  
  if(is.null(args_ADMM)){
    args_ADMM=default_ADMM_args
  }else{
    for (narg in names(default_ADMM_args)){
      if(is.null(args_ADMM[[narg]])){
        args_ADMM[[narg]]=default_ADMM_args[[narg]]
      }
    }
  }
  ### unpack arguments
  unpackList(args_ADMM, scope="L")
  
  stopwatch=Sys.time()
  if (is.null(eta)){
    eta=rho  ## default learning rate is the penalty for the augmented Lagrangian
  }
  normA<-sum(A^2)
  K<-ncol(B)
  N<-length(A)
  BtB<-t(B)%*%B
  lambda<-matrix(0, K, Time)
  Z<-matrix(0, N, Time)
  ### randomly initialize z and y
  lambda_dummy<-rnorm(K,1,1)
  Z[,1]<-A-B%*%lambda_dummy
  y<-matrix(0, N, Time)
  y[,1]<-rnorm(N,1,1)
  primal_res<-matrix(0,Time,1)
  primal_res[1]<-1
  dual_res<-matrix(0,Time,1)
  dual_res[1]<-1
  norm_res<-matrix(0,Time,1)
  t=1
  lambda[,1]<-lambda_dummy
  norm_res[1]<-sqrt(sum(sapply(Z[,1], FUN=function(x){
    return(x^2)
  })))
  R2<-matrix(0,Time,1)
  R2[1]<-norm_res[1]^2/normA
  converge=FALSE
  while (t <Time && converge==FALSE ){
    ### Update lambda
    if (dir==2){
      if (Z_penalty==1){
        Z[,t+1]<-y[,t]/rho+A-B%*%lambda[,t]
        Z[,t+1]<-sapply(Z[,t+1],FUN=function(x){
          return(soft_threshold(x,(Z_alpha)/rho))})
      }
      else{
        Z[,t+1]<-1/(rho+Z_alpha)*(y[,t]+rho*(A-B%*%lambda[,t]))
      }
      
      lambda[,t+1]<-solve(BtB)%*%t(B)%*%(y[,t]/rho-Z[,t+1]+A)
      lambda[,t+1]<-sapply(lambda[,t+1],FUN=function(x){
        return(soft_threshold(x,1/rho))})
    }
    else{
      lambda[,t+1]<-solve(BtB)%*%t(B)%*%(y[,t]/rho-Z[,t]+A)
      lambda[,t+1]<-sapply(lambda[,t+1],FUN=function(x){
        return(soft_threshold(x,1/rho))})
      if (Z_penalty==1){
        Z[,t+1]<-y[,t]/rho+A-B%*%lambda[,t+1]
        Z[,t+1]<-sapply(Z[,t+1],FUN=function(x){
          return(soft_threshold(x,Z_alpha/rho))})
      }
      else{
        Z[,t+1]<-1/(rho+Z_alpha)*(y[,t]+rho*(A-B%*%lambda[,t+1]))
      }
      
    }
    
    
    y[,t+1]<-y[,t]+eta*(A-B%*%lambda[,t+1]-Z[,t+1])
    primal_res[t+1]<-sum((lambda[,t+1]-lambda[,t])^2)##
    dual_res[t+1]<-sum((y[,t+1]-y[,t])^2)
    norm_res[t+1]<-sqrt(sum((Z[,t+1]-Z[,t])^2))
    R2[t+1]<-sum((A-B%*%lambda[,t+1])^2)/normA
    t=t+1
    
    #print(t)
    converge=( primal_res[t]<delta_1 && dual_res[t]/N<delta_2 && R2[t]<0.1)
    #print(converge)
  }
  #print(paste("stopped at iteration:",t))
  if (plot==TRUE){
    par(mfrow=c(2,2))
    plot(1:Time,primal_res,xlab="iteration number",ylab="value",main="primal residual",type="l")
    plot(1:Time,dual_res/N,xlab="iteration number",ylab="value",main="dual residual",type="l")
    plot(1:Time,norm_res,xlab="iteration number",ylab="value",main="Norm of difference in residual",type="l")
    plot(1:Time,R2,xlab="iteration number",ylab="value",main="R2",type="l")
  }

  time_taken=Sys.time()-stopwatch
  if (t<Time){
    sug.indice=t
  }
  else{
    gap_R2<-sapply(1:(Time-1), FUN=function(t){
      return(R2[t+1]-R2[t])
    })
    normlambda=apply(abs(lambda),2,sum)
    nn_zeros<-which(normlambda>10^(-3))
    explained<-which.min(R2)
    #print(explained)
    #print(nn_zeros)
    sug.indice=min(explained,max(which.max(gap_R2),5))
    ### Check to see that the lambda that we have selected is not entriely 0
   
  }

  
  return(list(lambda=lambda, Z=Z,y=y,primal_res=primal_res,R2=R2,dual_res=dual_res,suggested.indice=sug.indice,time_taken=time_taken,nb_it=t))
  
  
}


ADMM_with_elastic_net<-function(A,B,rho,alpha=0.1,eta=NULL,Time=10,dir=2,delta_1=0.01,delta_2=0.01,plot=TRUE){
  # A is a  n \times 1 vector
  # B is a $K \times n$ matrix
  # lambda should be a K dimensional vector
  if (is.null(eta)){
    eta=rho
  }
  K<-ncol(B)
  N<-length(A)
  normA<-sum(A^2)
  print(normA)
  BtB<-t(B)%*%B
  lambda<-matrix(0, K, Time)
  z<-matrix(0, N, Time)
  ### randomly initialize z and y
  lambda_dummy<-rnorm(K,1,1)
  z[,1]<-A-B%*%lambda_dummy
  y<-matrix(0, N, Time)
  y[,1]<-rnorm(N,1,1)
  primal_res<-matrix(0,Time,1)
  primal_res[1]<-1
  dual_res<-matrix(0,Time,1)
  norm_res<-matrix(0,Time,1)
  diff_res<-matrix(0,Time,1)
  R2<-matrix(0,Time,1)
  dual_res[1]<-1
  t=1
  lambda[,1]<-lambda_dummy
  diff_res[1]<-sum(z[,1]^2)
  norm_res[1]<-sum(z[,1]^2)
  R2[1]<-norm_res[1]/normA
  while (t <=(Time-1) ||( primal_res[t]>delta_1 && dual_res[t]>delta_2 ) ){
    ### Update lambda
    if (dir==2){
      z[,t+1]<-1/(rho+alpha)*(y[,t]+rho*(A-B%*%lambda[,t]))
      lambda[,t+1]<-solve(BtB)%*%t(B)%*%(y[,t]/rho-z[,t+1]+A)
      lambda[,t+1]<-sapply(lambda[,t+1],FUN=function(x){
        return(soft_threshold(x,1/rho))})
    }
    else{
      lambda[,t+1]<-solve(BtB)%*%t(B)%*%(y[,t]/rho-z[,t]+A)
      lambda[,t+1]<-sapply(lambda[,t+1],FUN=function(x){
        return(soft_threshold(x,1/rho))})
      z[,t+1]<-1/(rho+alpha)*(y[,t]+rho*(A-B%*%lambda[,t+1]))
      
    }
    
    
    y[,t+1]<-y[,t]+(A-B%*%lambda[,t+1]-z[,t+1])
    primal_res[t+1]<-sum((lambda[,t+1]-lambda[,t])^2)##
    dual_res[t+1]<-sum((y[,t+1]-y[,t])^2)
    norm_res[t+1]<-sum((z[,t+1]-z[,t])^2)
    diff_res[t+1]<-sqrt(sum(z[,t+1]^2))
    R2[t+1]<-sum((A-B%*%lambda[,t+1])^2)/normA
    t=t+1
    print(t)
  }
  print(paste("stopped at iteration:",t))
  if (plot==TRUE){
    par(mfrow=c(2,2))
    plot(1:Time,primal_res,xlab="iteration number",ylab="value",main="primal residual",type="l")
    plot(1:Time,dual_res,xlab="iteration number",ylab="value",main="dual residual",type="l")
    plot(1:Time,R2,xlab="iteration number",ylab="value",main="R2",type="l")
    plot(1:Time,diff_res,xlab="iteration number",ylab="value",main="Norm of difference in residual",type="l")
  }
  return(list(lambda=lambda, z=z,y=y,primal_res=primal_res,dual_res=dual_res,R2=R2,norm_res=norm_res,diff_res=diff_res,suggested.indice= which.max(norm_res),nb_it=t))
  
  
}





############################# Tests for the ADMM functions ######################################
test_ADMM<-function(N,K){
  A<-matrix(rnorm(N^2,0,1),N^2,1)
  B<-matrix(rnorm(N^2,0,1),N^2,1)
  for (k in 1:(K-1)){
    B_k<-matrix(rnorm(N^2,0,1),N^2,1)
    B<-cbind(B,B_k)
  }
  B_last=A
  index<-sample(N^2,ceiling(N^2/10)) ### sample 10% of the nodes and change them
  B_last[index]<-rnorm(length(index),0,1)
  B<-cbind(B,B_last)
  return(ADMM(A,B,rho=1,Time=10))
}



