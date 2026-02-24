#mc_sam: MCMC object from JAGS model, containing the variance of Y and variance of random effects
#resp: response vector (dim: nobs x 1)
#clust: vector with id's of clusters (dim: nobs x 1)
#X: design matrix for fixed part  (dim: nobs x nbeta)
#Z: design matrix for random effects (dim: nobs x nzeta)
#betas: name of the regression parameters in the MCMC object (string)
#var_y: Name of variance of error parameter in MCMC object (string)
#varZ: name of the covariance matrix of the random effects in the JAGS file (string)

LMM_Bayesselect <- function(mc_sam, resp, clust, X, Z, betas, var_y, varZ)
{
  resp <- as.matrix(resp)
  nobs <- nrow(resp)*ncol(resp)
  Z <- as.matrix(Z)
  if(nrow(Z)==1)
    Z <- as.matrix(rep(1,nobs*ncol(resp)))
  
  names_mcmc <- dimnames(mc_sam[[1]])[[2]]
  
  if(length(betas)>1)
  {
    sam_beta <- as.matrix(mc_sam[,betas])
  }else{
    sam_beta <- as.matrix(mc_sam[,grep(betas,names_mcmc)])
  }
  
  sam_mean <- sam_beta %*% t(X)
  if(ncol(resp)>1)
    clust <- rep(1:nrow(resp),ncol(resp))
  
  param_in1 <- names_mcmc[names_mcmc%in%var_y]
  if(length(param_in1)==0)
    print(paste("Please save in MCMC also",var_y))
  
  sam_vary <- as.matrix(mc_sam[,var_y])
  
  sam_meanb <- as.matrix(mc_sam[,grep("mu",names_mcmc)])
  param_in2 <- grep(varZ,names_mcmc)
  if(length(param_in2)==0)
    print(paste("Please save in MCMC also",varZ))
  
  if(ncol(Z)>1)
  {
    sam_varz <- as.matrix(mc_sam[,grep(paste0(varZ,"\\["),names_mcmc)])
  }else{
    sam_varz <- as.matrix(mc_sam[,varZ])
  }
  nK <- nrow(sam_beta) #Final length of MCMC with all chains
  
  clustl <- split(seq(nobs), clust) #cluster ID list
  
  nclus <- length(clustl)
  
  #intitialize vectors
  
  ### Deviance over posterior means
  
  sam_varz_p <- colMeans(sam_varz)
  sam_vary_p <- mean(sam_vary)
  
  logp_yi<- rep(0,nK)
  logp_yib<- rep(0,nK)
  
  loopi <- function(i,sam_mean,sam_meanb,nK,clustl,Z,resp,sam_varz_p,sam_vary_p,sam_varz,sam_vary){ #Loop for observations
    mu_yi <- matrix(sam_mean[,unlist(clustl[i])], nrow=nK)
    mub_yi<-matrix(sam_meanb[,unlist(clustl[i])], nrow=nK)
    Z_i <- matrix(Z[unlist(clustl[i]),], ncol=ncol(Z))
    resp_i <- resp[unlist(clustl[i])]
    indav <- which(!is.na(resp_i)) #Indicator of available response
    
    #for posterior means
    mu_yi_p <- colMeans(sam_mean)[unlist(clustl[i])]
    mu_yi_pb<- colMeans(sam_meanb)[unlist(clustl[i])]
    
    marvar_p <- Z_i %*% matrix(sam_varz_p, ncol=ncol(Z)) %*% t(Z_i) + diag(sam_vary_p, nrow(Z_i))
    logL_pm <- dmvnorm(resp_i[indav], mu_yi_p[indav], sigma=as.matrix(marvar_p[indav,indav]), log=TRUE)
    logL_pmb <- dmvnorm(resp_i[indav], mu_yi_pb[indav], sigma=as.matrix(diag(sam_vary_p, nrow(Z_i))[indav,indav]), log=TRUE)
    
    like<-function(k,Z_i,sam_varz,Z,sam_vary,resp_i,mu_yi,mub_yi,indav){ #Loop for MCMC iterations
      marvar <- Z_i %*% matrix(sam_varz[k,], ncol=ncol(Z)) %*% t(Z_i) + diag(sam_vary[k], nrow(Z_i))
      p_yi <-dmvnorm(resp_i[indav], mu_yi[k,indav], sigma=as.matrix(marvar[indav,indav]))
      p_yib<-dmvnorm(resp_i[indav], mub_yi[k,indav],sigma=as.matrix(diag(sam_vary[k], nrow(Z_i))[indav,indav]))
      c(p_yi,p_yib,log(p_yi),log(p_yib),1/p_yi,1/p_yib)
    }
    
    result <- sapply(c(seq(1,nK,1)),like,Z_i=Z_i,sam_varz=sam_varz,Z=Z,sam_vary=sam_vary,resp_i=resp_i,mu_yi=mu_yi,mub_yi=mub_yi,indav=indav)
    logp_yi <- apply(cbind(logp_yi,-2*result[3,]),1,sum,na.rm = TRUE)              
    logp_yib <- apply(cbind(logp_yib,-2*result[4,]),1,sum,na.rm = TRUE)  
    resultmean <- apply(result,1,mean) # vector with means of c(p_yi, p_yib, log(p_yi), log(p_yib), 1/p_yi, 1/p_yib)
    resultvar <- apply(result,1,var)# vector with variances of c(p_yi, p_yib, log(p_yi), log(p_yib), 1/p_yi, 1/p_yib)
    
    logCPOi <- log(1/resultmean[5])
    logCPObi <- log(1/resultmean[6])
    logLm <- resultmean[3]
    logLmb <- resultmean[4]
    log_mLm <- log(resultmean[1])
    log_mLmb <- log(resultmean[2])
    varp_yi <- resultvar[3] #variance log p_yi
    varp_yib <- resultvar[4] #variance log p_yib
    c(logLm,logLmb,log_mLm,log_mLmb,logL_pm,logL_pmb,varp_yi,varp_yib,logCPOi,logCPObi,logp_yi,logp_yib)
  }
  
  resloopi <- sapply(c(seq(1,nclus,1)),loopi,sam_mean=sam_mean,sam_meanb=sam_meanb,nK=nK,clustl=clustl,Z=Z,resp=resp,sam_varz_p=sam_varz_p,sam_vary_p=sam_vary_p,sam_varz=sam_varz,sam_vary=sam_vary)
  
  sumres <- apply(resloopi,1,sum) # vector with means of c(logLm,logLmb,log_mLm,log_mLmb,logL_pm,logL_pmb,varp_yi,varp_yib,logCPOi,logCPObi,logp_yi,logp_yib) where logp_yi and lop_yib are vectors of 1 X niteration size, i.e, probabilities for each sample of \theta.
  
  mid <- (length(sumres)-10)/2
  
  #mean deviance
  mean(sumres[11:(10+mid)])
  mean(sumres[(11+mid):(mid*2)])
  
  #variance deviance for Gelman DIC calculation
  pv <- var(sumres[11:(10+mid)])/2
  pvb <- var(sumres[(11+mid):(10+mid*2)])/2
  
  Devm <- -2*sumres[1] #LogLm
  Dev_pmW <- -2*sumres[3] #log_mLm
  Dev_pm <- -2*sumres[5] #sum(logL_pm)
  
  Devmb <- -2*sumres[2] #LogLmb
  Dev_pmWb <- -2*sumres[4]#log_mLmb
  Dev_pmb <- -2*sumres[6] #sum(logL_pmb)
  
  pDm <- Devm - Dev_pm
  DICm1 <- Devm + pDm #Spiegelhalter
  DICm2<-Devm+pv#Gelman
  DICm3<-Devm+pDm*log(nobs)
  DICm4<-Devm+pv*log(nobs)
  DICm5<-Devm+pDm*log(nclus)
  DICm6<-Devm+pv*log(nclus)
  
  
  pDmb <- Devmb - Dev_pmb
  DICmb1 <- Devmb + pDmb
  DICmb2<-Devmb + pvb
  DICmb3<- Devmb + pDmb*log(nobs)
  DICmb4<- Devmb + pvb*log(nobs)
  DICmb5<- Devmb + pDmb*log(nclus)
  DICmb6<- Devmb + pvb*log(nclus)
  
  
  pDW2<-sumres[7] #varp_yi
  pDW <- Devm - Dev_pmW
  WAICm <- Dev_pmW+2*pDW
  WAICm2 <- Dev_pmW+2*pDW2
  
  pDWb2<-sumres[8] #varp_yib
  pDWb <- Devmb - Dev_pmWb
  WAICmb <- Dev_pmWb+2*pDWb
  WAICmb2 <- Dev_pmWb+2*pDWb2 
  
  oPSBF=sumres[9] # sum logCPOi
  #marginal
  lppd<- -2*sumres[9]
  mplppd<-lppd-Dev_pmW
  mlppd<-lppd+mplppd
  
  PSBF=sumres[10] #sum logCPObi
  #conditional
  lppdb <- -2*sumres[10]
  cplppd <- lppdb-Dev_pmWb
  clppd <- lppdb+cplppd
  
  return(list( mpD_DIC=pDm,mDIC1=DICm1, mDIC2=DICm2, mDIC3=DICm3, mDIC4=DICm4,mDIC5=DICm5, mDIC6=DICm6, mWAIC=WAICm, mpD_WAIC=pDW, mWAIC2=WAICm2, mpD_WAIC2=pDW2, pv=pv,mPSBFn=-oPSBF,cpD_DIC=pDmb,cDIC1=DICmb1,cDIC2=DICmb2,cDIC3=DICmb3,cDIC4=DICmb4,cDIC5=DICmb5,cDIC6=DICmb6,cWAIC=WAICmb, cpD_WAIC=pDWb, cWAIC2=WAICmb2, cpD_WAIC2=pDWb2, pvb=pvb,cPSBFn=-PSBF,mplppd= mplppd/2, mlppd= mlppd,cplppd=cplppd/2,clppd=clppd))
}
