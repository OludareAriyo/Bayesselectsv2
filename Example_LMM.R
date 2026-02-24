
##Libraries
library("mvtnorm")
library("rjags")

setwd("JAS_fuctions")

##Call the function that calculates the criteria
source("func_bayesselect.R")

###
### Simulating dataset
###

N <- 50     #Number of patients
simu <- 10  #Iteration number in simulation

vI <- 2.9932 #Variance of random intercept
Sigmae <- 2.0242 #Residual variance
Beta1 <- c(24.9688, -2.3210,  1.4831) #Regression coefficients

#Simulate covariates
set.seed(simu*15+5)
male <- rep(rbinom(N,1,0.6),each=4)
age <- scale(rep(c(8,10,12,14),N))

#Simulate random effects
set.seed(simu*10+4)
b0 <- rnorm(N,0,sqrt(vI))
vb0 <- matrix(rep(b0,each=4),ncol=1)

#Simulate response
mean.y <- cbind(1,male,age)%*%Beta1
mean.yc <- mean.y+vb0
set.seed(simu*15+5)
ygen <- rnorm(N*4, mean.yc, sqrt(Sigmae))
id <- sort(rep(1:N,4))

#Data for JAGS
simdata <- list(y=ygen, sex=male, age=age[,1], id=id, K=length(ygen), n=max(id))

### Run JAGS model
results <- jags.model(file="JAGS_Model1.txt", data=simdata, n.chains=3)
update(results, n.iter=5000)
results_s <- coda.samples(results, c("beta","sigma_b","sigma_e","var_e","mu"), n.iter=10000, thin=10)

### Compute criteria
Bayesselect(results_s, resp=simdata$y, clust=simdata$id, X=cbind(1,simdata$sex,simdata$age),
            Z=cbind(1,simdata$age), betas=c("beta[1]","beta[2]","beta[3]"), var_y="var_e", varZ="sigma_b")


