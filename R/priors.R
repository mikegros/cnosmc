# Functions in this file:
#
# Functions to evaluate a log prior:
# LogpriorBstring   # Bernoulli
# Logpriorg         # Uniform(0,10)  because it should be positive (inhibitors are dealt with in other ways)
# Logpriorn         # exponential
# Logpriork         # Uniform(0,10)  somewhat arbitrarily set just because it should be positive.  Maybe Chiquared(3) would be ok?
#
# Functions to sample from the priors:
# SamplepriorBstring
# Samplepriorg
# Samplepriorn
# Samplepriork
#
#
# NOTE:
# All of these can easily be replaced to change the priors
#
#


############Evaluate the priors#######


############
LogpriorGstring = function(initGstring,p=rep(.9,length(initGstring))){
  sum(initGstring*log(p) + (1-initGstring)*log(1-p))
}
############
Logpriorg = function(gCube,lower=0,upper=1){
  # uniform prior
  ifelse(min(gCube)>=lower & max(gCube)<=upper, 0, -Inf)
}
############
Logpriorn = function(nCube,lambda=2){
  # exponential prior with mean lambda is proportional to
  -sum(nCube)/lambda
}
############
Logpriork = function(kCube,lower=0,upper=1){
  # uniform prior
  ifelse(min(kCube)>=lower & max(kCube)<=upper, 0, -Inf)
}
############
Logpriorsigsq = function(sigsq,alpha,beta){
  # inverse gamma prior
  -dgamma(1/sigsq,alpha,beta,log = TRUE)
}
############




############Sample from the priors#######


############
SamplepriorGstring = function(n,p=rep(0.9,n)){
  rbinom(n,1,prob=p)
}
############
Samplepriorg = function(n,lower=0,upper=10){
  # uniform prior
  runif(n,min=lower,max=upper)
}
############
Samplepriorn = function(n,lambda=.2){
  # exponential prior with mean lambda is proportional to
  rexp(n,1/lambda)
}
############
Samplepriork = function(n,lower=0,upper=10){
  # uniform prior
  runif(n,min=lower,max=upper)
}
############




