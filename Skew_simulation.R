library(sn)
library(Matrix)
library(parallel)
library(doSNOW)  # To compute oDIC
library(foreach) # in parallel
library(dclone)
library(MCMCpack)
library(Matrix)
library(MASS)
library("mvtnorm")
library("rjags")
library(matlib)
library(expm)
library(mvtnorm)
setwd("I:\My work\Research materials\leuvenPhD\Paper2\proposal\ModelCopare\JAStemplate\JAS new\Simulation_Study\Skew-normal code\JAS_fuctions")
##Call function that calculates the criteria
source("SN_bayeselect.R")
N=200
n=5
simu=1
Sigmae=0.5
Sigmab=0.2
deltab=4 #1.3
deltae=1.6 
#male <- rep(rbinom(N,1,0.6),each=4)
#age <- scale(rep(c(8,10,12,14),N))
i=rep(c(1:N),each=n)
j=rep(c(1:n),N)
t=j-3
x=c(rep(1,N*n/2),rep(0,N*n/2))
Beta1=c(4,2,1)
set.seed(simu*15+5)
#b=rsn(N,0,Sigmab^2,deltab)
b=rsn(N,0,Sigmab^2)
vb=rep(b,each=n)
mean.y <- cbind(rep(1,N*n),x,t)%*%Beta1
mean.yc <- mean.y+vb
id <- sort(rep(1:N,n))
set.seed(simu*9+7)
ygen <- rsn(N*n, mean.yc, sqrt(Sigmae),deltae)  #[1:(n*N)]
data=data.frame(y=ygen,x=x,t=t,i=i,j=j,id=id)

flags <- rep(0,2)
criters <- matrix(0, ncol=2, nrow=15)

#monitor:
model=list()
model[[1]]="model.txt"
model[[2]]="model2.txt"
model[[3]]="model3.txt"

monitor=list()
monitor[[1]]= c("beta","taub","tau","deltab","deltae","sigma_e","var_e","sigma_b","mu")
monitor[[2]]= c("beta","taub","tau","deltab","deltae","sigma_e","var_e","var_b","mu")
monitor[[3]]= c("beta","taub","tau","deltab","deltae","sigma_e","var_e","var_b","mu")


#betas:
betas=list()
betas[[1]]=c("beta[1]","beta[2]","beta[3]")
betas[[2]]=c("beta[1]","beta[2]","beta[3]")
betas[[3]]=c("beta[1]","beta[2]","beta[3]")


#varZ
varZ=list()
varZ[[1]]="sigma_b"
varZ[[2]]="var_b"
varZ[[3]]="var_b"

X=list()
X[[1]]=cbind(1,data$x,data$t)
X[[2]]=cbind(1,data$x,data$t)
X[[3]]=cbind(1,data$x,data$t)

#Z:
Z=list()
Z[[1]]= cbind(1,data$t)
Z[[2]]=cbind(rep(1,length(data$t)))
Z[[3]]=cbind(rep(0,length(data$t)))


modeli <- 1
for (em in c(1,2,3)){
  if(em==1){
    jags.inits <-function(){list(beta=rnorm(3),taub=1,tau=1,deltae=3,sigma_e=2, deltab=c(1,1),v = structure(.Data = c(0.1, 0, 0, 0.1), .Dim = c(2, 2)))}
                                 
                                 #deltab=c(1,1),w=2,v = structure(.Data = c(0.1, 0, 0, 0.1), .Dim = c(2, 2)))
    }else{jags.inits <- function(){list(beta=rnorm(3),taub1=runif(1,0.5,1),tau=1, deltae=1,deltab=1)}
       }
    
  jdata <- list(y=ygen, x=x, t=t, N=N, n=n,id=id,K=N*n, R=diag(c(0.001,0.001)))  
  
  results <- jags.model(file=model[[em]], data=jdata,inits, n.chains=3)

  update(results, n.iter=100, n.burnin=50, n.thin=10)
  results_s <- coda.samples(results, monitor[[em]], n.iter=100, thin=10)
 
    bgrM <- gelman.diag(results_s, multivariate=FALSE)$psrf[,1]
  
  
  
  RA1 <- SN_bayeselect(results_s, resp=data$y, clust=data$id, X=X[[em]],Z=Z[[em]], betas=betas[[em]], var_y="var_e", varZ=varZ[[em]])
  RA1=unlist(RA1)
  RA=c(RA1)

  criters[,modeli] <- RA
  flags[modeli] <- flag
  modeli <- modeli + 1
}   



