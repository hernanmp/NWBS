rm(list = ls())

working_directory = "/Users/oscar/Documents/GitHub/NWBS"
setwd(working_directory)
#setwd(paste(working_directory,"/Code_J",sep=""))
source("utils_functions.R")
source('utils.R')

library(cumSeg)
library(strucchange)
library(Segmentor3IsBack)
library(changepoint)
library(stepR)
library(wbs)
library(changepoint.np)
M =  120
NMC = 100

n = 1
T_grid = c(8000,4000,1000)
###  Generate data
sigma = 1/sqrt(3)
tau_grid  =    c(seq(0.9,3.5,length = 15),seq(3.6,4.5,length = 5),100)

#####

hausdorff_cp =   matrix(0,NMC,  length(T_grid))
hausdorff_np_bs =   matrix(0,NMC,  length(T_grid))
hausdorff_np =   matrix(0,NMC,  length(T_grid))
hausdorff_wbs =   matrix(0,NMC,  length(T_grid))
hausdorff_smuce =   matrix(0,NMC,  length(T_grid))
hausdorff_hsmuce =   matrix(0,NMC,  length(T_grid))
hausdorff_segmentor =   matrix(0,NMC,  length(T_grid))
hausdorff_breakpoints =   matrix(0,NMC,  length(T_grid))
hausdorff_pelt =   matrix(0,NMC,  length(T_grid))

hausdorff2_cp =   matrix(0,NMC,  length(T_grid))
hausdorff2_np_bs =   matrix(0,NMC,  length(T_grid))
hausdorff2_np =   matrix(0,NMC,  length(T_grid))
hausdorff2_wbs =   matrix(0,NMC,  length(T_grid))
hausdorff2_smuce =   matrix(0,NMC,  length(T_grid))
hausdorff2_hsmuce =   matrix(0,NMC,  length(T_grid))
hausdorff2_segmentor =   matrix(0,NMC,  length(T_grid))
hausdorff2_breakpoints =   matrix(0,NMC,  length(T_grid))
hausdorff2_pelt =   matrix(0,NMC,  length(T_grid))

error_cp =   matrix(0,NMC,  length(T_grid))
error_np_bs =   matrix(0,NMC,  length(T_grid))
error_np =   matrix(0,NMC,  length(T_grid))
error_wbs =   matrix(0,NMC,  length(T_grid))
error_smuce =   matrix(0,NMC,  length(T_grid))
error_hsmuce =   matrix(0,NMC,  length(T_grid))
error_segmentor =   matrix(0,NMC,  length(T_grid))
error_breakpoints =   matrix(0,NMC,  length(T_grid))
error_pelt =   matrix(0,NMC,  length(T_grid))

T_ind =1
iter= 1
for(T_ind in 1: length(T_grid) )
{
  print("T_ind")
  for(iter in 1:NMC)
  {
    print("iter")
    #n =  fixed_n_grid[n_ind]
    T =   T_grid[T_ind]
    
    
    K = 5#floor(T^{1/2}  /  sqrt(log(T)))
    n = 1
    #T = 1000
    theta =  rep(.2,K+1);
    theta[2*(1:floor((K+1)/2))] =  1
    
    # temp = radiological_example(K, celsium  = 1,T)
    # background_train =  temp$background_train
    # f0  = temp$f0
    # p0 =  background_train # conventional background
    # f0[K+1,] =f0[1,]
    
    v =  rep(0,K)
    for(i in 1:K)
    {
      v[i] = i*floor(T/(K+1))
    } 
    ############generate data
    N  =  rep(n,T)   #rpois(T, 10)
    n =  max(N)
    y =  matrix(NA, T,n)
    
    for(t in 1:T)
    {
      if( t<= v[1])
      {
        y[t,1:N[t]]  =  theta[1]*rnorm(1)  
        #sigma*rt(N[t],3)
        #sample(x= (1:length(background_train))/length(background_train),size= N[t], prob = f0[1,]) 
      }
      if(K >1)
      {
        for(j in 1:(K-1))
        {
          if(  t >  v[j]   &&  t <=  v[j+1])
          {
            y[t,1:N[t]]  =   theta[j+1]*rnorm(1)
            #sample(x= (1:length(background_train))/length(background_train),size=  N[t], prob = f0[j+1,]) 
          }
        }
      }
      if(t > v[K])
      {
        y[t,1:N[t]]  = theta[K+1]*rnorm(1)
        #sample(x= (1:length(background_train))/length(background_train),size= N[t], prob = f0[K+1,]) 
      }
    }###############
    y_raw =  c()
    
    for(t in 1:T)
    {
      y_raw =  c(y_raw, y[t,1:N[t]])
    }
    
    ###########
    ###  hsmuce
    
    fit <- stepFit(as.vector(y), x = 1:dim(y)[1], alpha = 0.5, jumpint = TRUE, confband = TRUE,family = "hsmuce")
    est_hsmuce =  sort(setdiff(fit$rightEnd,c(1,T)))
    hausdorff_hsmuce[iter,T_ind] =  dist_change_points(est_hsmuce,v)
    hausdorff2_hsmuce[iter,T_ind] =  dist_change_points(v,est_hsmuce)
    error_hsmuce[iter,T_ind] =  K - length(est_hsmuce) 
    ##############
    
    
    out <- cpt.np(as.vector(y),method="PELT",nquantiles =10)
    temp = cpts(out)
    
    hausdorff_cp[iter,T_ind] =  dist_change_points(temp,v)
    hausdorff2_cp[iter,T_ind] =  dist_change_points(v,temp)
    error_cp[iter,T_ind] = length(v) - length(temp)
    ##############33333  
    
    min_y = min(y_raw)
    max_y = max(y_raw)
    grid =   seq(min_y,max_y, length = M)
    
    
    s= 1
    e =  T
    flag = 0 
    
    S =  NULL
    
    ######################################
    
    
    m = floor(T/2)
    y_new =  cbind(y[2*(1:m)-1],y[2*(1:m)])
    
    S = NBS_full(as.matrix(y_new[,1]),as.matrix(y_new[,2]),gam=2,rep(1,m))
    S = 2*S
    
    S =  sort(unique(S))
    hausdorff_np_bs[iter,T_ind] =  dist_change_points(S,v)
    hausdorff2_np_bs[iter,T_ind] =  dist_change_points(v,S)
    error_np_bs[iter,T_ind] =  K - length(S)
    
    #################
    #### NWBS
    
    
    M =   120
    alpha =  sample.int(size =M  , n = m,replace = TRUE)
    beta =   sample.int(size =M  , n = m,replace = TRUE)#alpha + floor((T- alpha)*runif(M))
    #
    for(j in 1:M)
    {
      aux =  alpha[j]
      aux2 =  beta[j]
      #
      alpha[j] = min(aux,aux2)
      beta[j] = max(aux,aux2)
    }
    
    S2 = NWBS_full(as.matrix(y_new[,1]),as.matrix(y_new[,2]),gam=2,rep(1,m),alpha,beta)
    S2 = 2*S2 
    S2 =  sort(unique(S2))
    
    
    hausdorff_np[iter,T_ind] =  dist_change_points(S2,v)
    hausdorff2_np[iter,T_ind] =  dist_change_points(v,S2)
    error_np[iter,T_ind] =  K - length(S2)
    
    
    ############################################3
    
    w <- wbs(y)
    w.cpt <- changepoints(w,penalty="bic.penalty")
    wbs_est = sort( w.cpt$cpt.ic$bic.penalty)
    
    
    hausdorff_wbs[iter,T_ind] =  dist_change_points(wbs_est,v)
    hausdorff2_wbs[iter,T_ind] =  dist_change_points(v,wbs_est)
    error_wbs[iter,T_ind] =  K - length(wbs_est)
    
    ##########################
    bp.nile <- breakpoints(y ~ 1)
    #bp.nile$breakpoints
    if(is.na(bp.nile$breakpoints[1])==TRUE)
    {
      bp.nile$breakpoints = NULL
    }
    breakpoints_est =  bp.nile$breakpoints
    hausdorff_breakpoints[iter,T_ind] =  dist_change_points(breakpoints_est,v)
    hausdorff2_breakpoints[iter,T_ind] =  dist_change_points(v,breakpoints_est)
    
    error_breakpoints[iter,T_ind] =  K - length(breakpoints_est)
    
    ######################
    
    z <- Segmentor(as.vector(y), model=2)
    ind  = SelectModel(z)
    temp  = as.numeric(z@breaks[ind,])
    temp  =  temp[intersect(which(temp>0) ,which(temp < T))]
    est_S =   c()
    if(length(temp) >0)
    {
      est_S =  c(est_S,temp[1])
      if(length(temp)>1)
      {
        for(j in 2:length(temp))
        {
          if(abs(temp[j] -  temp[j-1]  ) > 10)
          {
            est_S =  c(est_S, temp[j])
          }
        }
      }
    }
    
    hausdorff_segmentor[iter,T_ind] =  dist_change_points(est_S,v)
    hausdorff2_segmentor[iter,T_ind] =  dist_change_points(v,est_S)
    error_segmentor[iter,T_ind] =  K - length(est_S)
    ###############################
    
    
    
    
    temp  = cpt.mean(as.vector(y)/mad(diff(as.vector(y))/sqrt(2)), method="PELT")@cpts
    temp  =  temp[intersect(which(temp>0) ,which(temp < T))]
    est_pelt =   c()
    if(length(temp) >0)
    {
      est_pelt =  c(est_pelt,temp[1])
      if(length(temp)>1)
      {
        for(j in 2:length(temp))
        {
          if(abs(temp[j] -  temp[j-1]  ) > 10)
          {
            est_pelt =  c(est_pelt, temp[j])
          }
        }
      }
    }
    
    hausdorff_pelt[iter,T_ind] =  dist_change_points(est_pelt,v)
    hausdorff2_pelt[iter,T_ind] =  dist_change_points(v,est_pelt)
    error_pelt[iter,T_ind] =  K - length(est_pelt)
    ###################################3
    ################
    
    temp = smuceR(y, 1:length(y), family="gauss")
    indices =  which(abs(diff(fitted(temp)))>0)
    est_smuce =  c()
    if(length(indices)>0)
    {
      est_smuce= c(indices[1])
    }
    if(length(indices)>1)
    {
      for(j in 2:length(indices))
      {
        if(indices[j] -   indices[j-1] > 1  )
        {
          est_smuce =  c(est_smuce,indices[j])
        }
      }
    }
    
    hausdorff_smuce[iter,T_ind] =  dist_change_points(est_smuce,v)
    hausdorff2_smuce[iter,T_ind] =  dist_change_points(v,est_smuce)
    error_smuce[iter,T_ind] =  K - length(est_smuce)   
    ###################
    
    
  }###  close loop for iter
  
  
  print("cp")
  print(median(hausdorff_cp[,T_ind]))
  print(median(hausdorff2_cp[,T_ind]))
  print(mean(abs(error_cp[,T_ind])))
  
  print("NP")
  print(median(hausdorff_np[,T_ind]))
  print(median(hausdorff2_np[,T_ind]))
  print(mean(abs(error_np[,T_ind])))
  
  print("NP_bs" )
  print(median(hausdorff_np_bs[,T_ind]))
  print(median(hausdorff2_np_bs[,T_ind]))
  print(mean(abs(error_np_bs[,T_ind])))
  print("wbs")
  print(median(hausdorff_wbs[,T_ind]))
  print(median(hausdorff2_wbs[,T_ind]))
  print(mean(abs(error_wbs[,T_ind])))
  print("bp")
  print(median(hausdorff_breakpoints[,T_ind]))
  print(median(hausdorff2_breakpoints[,T_ind]))
  print(mean(abs(error_breakpoints[,T_ind])))
  print("segmentor")
  print(median(hausdorff_segmentor[,T_ind]))
  print(median(hausdorff2_segmentor[,T_ind]))
  print(mean(abs(error_segmentor[,T_ind])))
  print("pelt")
  print(median(hausdorff_pelt[,T_ind]))
  print(median(hausdorff2_pelt[,T_ind]))
  print(mean(abs(error_pelt[,T_ind])))
  print("smuce")
  print(median(hausdorff_smuce[,T_ind]))
  print(median(hausdorff2_smuce[,T_ind]))
  print(mean(abs(error_smuce[,T_ind])))
  print("hsmuce")
  print(median(hausdorff_hsmuce[,T_ind]))
  print(median(hausdorff2_hsmuce[,T_ind]))
  print(mean(abs(error_hsmuce[,T_ind])))
}