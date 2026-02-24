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
source("SN_bayesselect.R")

N=200
n=5
Beta=c(4,2,1)
simu=1
Sigmae=0.5
Sigmab=0.2
deltab=0.5#1.3
deltae=0.6  #1.6
vI= 3.2668
i=rep(c(1:N),each=n)
j=rep(c(1:n),N)
t=j-3
#t<- scale(rep(c(8,10,12,14,16),N))
#x=c(rep(1,N*n/2),rep(0,N*n/2))
set.seed(simu*15+5)
x <- rep(rbinom(N,1,0.5),each=n)
#set.seed(simu*15+5)
#b=rsn(N,0,Sigmab^2,deltab)
#b=rsn(N,0,Sigmab^2)
#vb=rep(b,each=n)
set.seed(simu*10+4)
#b0= rnorm(N,0,sqrt(vI))
b0=rsn(N,0,Sigmab^2,deltab)
vb0 <- matrix(rep(b0,each=n),ncol=1)
mean.y <- cbind(rep(1,N*n),x,t)%*%Beta
mean.yc <- mean.y+vb
id <- sort(rep(1:N,n))
set.seed(simu*9+7)
ygen <- rnorm(N*n, mean.yc, sqrt(Sigmae))
data=data.frame(y=ygen,x=x,t=t,i=i,j=j,id=id)



#################################
#library(Hmisc)
#rcorr(as.matrix(data)) 
#cor(data$y,data$t)
################################
#simdata=data.frame()
inits1=list(beta1=0, beta2=0, beta3=0, taub1=runif(1,0.5,1),tau=1,deltae=0)
           # deltab=c(1,1) #v = structure(.Data = c(0.1, 0, 0, 0.1), .Dim = c(2, 2))
inits2=list(beta1=0, beta2=0, beta3=0, taub1=runif(1,0.5,1),tau=1, deltae=0)
           # v = structure(.Data = c(0.3, 0, 0, 0.2), .Dim = c(2, 2))
inits3=list(beta1=0, beta2=0, beta3=0, taub1=runif(1,0.5,1),tau=1, deltae=0)
jags.inits <-function(){list(beta=rnorm(3),tau=1,deltae=3,sigma_e=2, deltab=c(1,1),v = structure(.Data = c(0.1, 0, 0, 0.1), .Dim = c(2, 2)))}         #v = structure(.Data = c(0.1, 0, 0, 0.1), .Dim = c(2, 2)),delta=0)

inits= list(inits1,inits2,inits3)
jdata <- list(y=ygen, x=x, t=t, N=N, n=n,id=i,K=N*n,R=diag(c(0.001,0.001)))
results <- jags.model(file="model.txt", data=jdata,inits,n.chain=3,quiet = TRUE)
update(results, n.iter=70000, n.burnin=3000, n.thin=10)
results_s <- coda.samples(results, c("beta","tau","deltae","sigma_e","var_e","sigma_b","mu","deltab"), n.iter =150)  #,"deltab"


### Compute criteria
SN_bayeselect(results_s,resp=data$y,clust=data$id, X=cbind(1,data$x,data$t),
            Z=cbind(1,data$t), betas=c("beta[1]","beta[2]","beta[3]"), var_y="var_e", varZ="sigma_b",deltae="deltae",deltab="deltab")#,deltab="deltab"

