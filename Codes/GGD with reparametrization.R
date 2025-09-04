rm(list = ls()); gc()
require(VGAM)
require(msm)
require(numDeriv)
library(posterior)

#################################
####### Datos sintéticos #######
###############################
n = 100
real.phi = 3
real.mu = 2
real.alpha = 1
cen = 0.2 #Censura
set.seed(2023)
t <- rgengamma.stacy(n, scale = (1/real.mu),
                     d = real.alpha, k = real.phi) #Tiempo (t ~ GGD)
delta<-rep(1,n) #Será la indicadora de no censura
#Adiciona censura tipo I
y<-sample(1:n,n*cen) #Será el tiempo observado
for (j in 1:length(y)){
  t[y[j]]<-runif(1,0,t[j])
  delta[y[j]]<-0
}
d<-sum(delta) #Cantidad de datos no censurados

#################################
######### Funciones ############
###############################

####### log-likelihood
loglikeh <- function (theta){
  alpha = theta[1]
  mu = theta[2]
  phi = theta[3]
  epsilon  = 1e-19 #Allows evaluate the log function in 0's computational values
  likeh <- d*log(alpha)+(d*alpha*phi)*log(mu)-n*lgamma(phi)+(alpha*phi-1)*sum(delta*log(t))-(mu^alpha)*sum(delta*(t^alpha)) + 
    sum((1-delta)*(log(gamma(phi))+log(epsilon + pgamma((mu*t)^alpha,phi, lower=FALSE))))
  return(likeh)
}
#loglikeh(c(1,2,3))

####### log-posterior
logpost <- function (theta){
  alpha = theta[1]
  mu = theta[2]
  phi = theta[3]
  logprior <- -log(alpha)-log(mu)+0.5*log((phi^2)*(trigamma(phi)^2)-trigamma(phi)-1)-0.5*log(phi+(phi^2)*trigamma(phi)-1)
  logposterior = logprior + loglikeh(theta)
  return(logposterior)
}
#logpost(c(1,2,3))
grad(logpost,c(1,2,3))

####### Reparametrization lambda = phi/(mu)^alpha
logpost.rep <- function(theta){
  alpha = theta[1]
  lambda = theta[2]
  phi = theta[3]
  logposterior.rep = logpost(c(alpha,(phi/lambda)^(1/alpha),phi)) + (1/alpha)*log(phi/lambda) - log(alpha) - log(lambda)
  return(logposterior.rep)
}
#logpost.rep(c(1,3/2,3))
grad(logpost.rep,c(1,3/2,3))

####### log-trasnformation eta_j = log(theta_j)
logpost.eta <- function(eta){
  eta1 = eta[1]
  eta2 = eta[2]
  eta3 = eta[3]
  logposterior.eta = logpost.rep(c(exp(eta1),exp(eta2),exp(eta3))) + eta1 + eta2 + eta3
  return(logposterior.eta)
}
#logpost.eta(c(log(1),log(3/2),log(3)))
grad(logpost.eta,c(log(1),log(3/2),log(3)))


#################################
##########   MALA   ############
###############################
h = 0.009 #Propose step size (variance of prop)
T = 3000 #Iterations
eta1sim = rep(0,T+1) #Vector for simulations of eta1
eta2sim = rep(0,T+1) #Vector for simulations of eta2
eta3sim = rep(0,T+1) #Vector for simulations of eta3
# Initial values
eta1sim[1] = log(1) 
eta2sim[1] = log(3/2)
eta3sim[1] = log(3)
#Rejection counters
reject = rep(0,T+1)

for (j in 2:(T+1)){
  etasim = c(eta1sim[j-1], eta2sim[j-1], eta3sim[j-1])
  current_gradient <- grad(logpost.eta, etasim)
  prop1 <- rnorm(1,eta1sim[j-1] + (0.5*h)*current_gradient[1],sqrt(h))
  prop2 <- rnorm(1,eta2sim[j-1] + (0.5*h)*current_gradient[2],sqrt(h))
  prop3 <- rnorm(1,eta3sim[j-1] + (0.5*h)*current_gradient[3],sqrt(h))
  prop = c(prop1,prop2,prop3)
  prop_gradient <- grad(logpost.eta, prop)
  ratio <- logpost.eta(prop) - logpost.eta(etasim) + 
    log(dnorm(eta1sim[j-1], prop1 + (0.5*h)*prop_gradient[1], sqrt(h))) +
    log(dnorm(eta2sim[j-1], prop2 + (0.5*h)*prop_gradient[2], sqrt(h))) +
    log(dnorm(eta3sim[j-1], prop3 + (0.5*h)*prop_gradient[3], sqrt(h))) -
    log(dnorm(prop1,eta1sim[j-1] + (0.5*h)*current_gradient[1], sqrt(h))) -
    log(dnorm(prop2,eta2sim[j-1] + (0.5*h)*current_gradient[2], sqrt(h))) -
    log(dnorm(prop3,eta3sim[j-1] + (0.5*h)*current_gradient[3], sqrt(h)))
  prob <- min(1,exp(ratio))
  u <- runif(1)
  if(u < prob){
    eta1sim[j] <- prop1
    eta2sim[j] <- prop2
    eta3sim[j] <- prop3
  }else{
    eta1sim[j] <- eta1sim[j-1]
    eta2sim[j] <- eta2sim[j-1]
    eta3sim[j] <- eta3sim[j-1]
    reject[j] <- 1
  }
}
1-sum(reject)/T #Acceptation rate
#head(cbind(eta1sim,eta2sim,eta3sim), 10)
par(mfrow=c(3,2))
ts.plot(eta1sim, ylab = ~log(alpha));abline(h=c(log(1),mean(eta1sim)),col=c("green","darkred"));acf(eta1sim)
ts.plot(eta2sim, ylab = ~log(lambda));abline(h=c(log(3/2),mean(eta2sim)),col=c("green","darkred"));acf(eta2sim)
ts.plot(eta3sim, ylab = ~log(phi));abline(h=c(log(3),mean(eta3sim)),col=c("green","darkred"));acf(eta3sim)
par(mfrow=c(1,3))
plot(eta1sim,eta2sim, xlab = ~log(alpha), ylab = ~log(lambda))
plot(eta1sim,eta3sim, xlab = ~log(alpha), ylab = ~log(phi))
plot(eta2sim,eta3sim, xlab = ~log(lambda), ylab = ~log(phi))

## Cleaning up
burn = T/3
thin = 40
etasim = as_draws(cbind(eta1sim,eta2sim,eta3sim)[(burn+1):T,])
muestra = thin_draws(exp(etasim), thin = thin)
par(mfrow=c(3,2))
ts.plot(muestra[,1], ylab = ~alpha);abline(h=c(1,mean(muestra[,1])),col=c("green","darkred"));acf(muestra[,1])
ts.plot(muestra[,2], ylab = ~lambda);abline(h=c(3/2,mean(muestra[,2])),col=c("green","darkred"));acf(muestra[,2])
ts.plot(muestra[,3], ylab = ~phi);abline(h=c(3,mean(muestra[,3])),col=c("green","darkred"));acf(muestra[,3])
par(mfrow=c(1,3))
plot(muestra[,1],muestra[,2], xlab = ~alpha, ylab = ~lambda)
plot(muestra[,1],muestra[,3], xlab = ~alpha, ylab = ~phi)
plot(muestra[,2],muestra[,3], xlab = ~lambda, ylab = ~phi)
par(mfrow=c(1,3))
hist(muestra[,1], freq = F, xlab = ~alpha, xlim = c(0,5))
hist(muestra[,2], freq = F, xlab = ~lambda, xlim = c(0,5))
hist(muestra[,3], freq = F, xlab = ~phi, xlim = c(0,5))


#################################
##########   RWMH   ############
###############################
sigma2p = 0.01 #propose variance
T = 3000 #Iterations
eta1sim = rep(0,T+1) #Vector for simulations of eta1
eta2sim = rep(0,T+1) #Vector for simulations of eta2
eta3sim = rep(0,T+1) #Vector for simulations of eta3
# Initial values
eta1sim[1] = log(1) 
eta2sim[1] = log(3/2)
eta3sim[1] = log(3)
#Rejection counters
reject = rep(0,T+1)

for (j in 2:(T+1)){
  etasim = c(eta1sim[j-1], eta2sim[j-1], eta3sim[j-1])
  prop1 <- rnorm(1,eta1sim[j-1],sqrt(sigma2p))
  prop2 <- rnorm(1,eta2sim[j-1],sqrt(sigma2p))
  prop3 <- rnorm(1,eta3sim[j-1],sqrt(sigma2p))
  prop = c(prop1,prop2,prop3)
  ratio <- logpost.eta(prop) - logpost.eta(etasim)
  prob <- min(1,exp(ratio))
  u <- runif(1)
  if(u < prob){
    eta1sim[j] <- prop1
    eta2sim[j] <- prop2
    eta3sim[j] <- prop3
  }else{
    eta1sim[j] <- eta1sim[j-1]
    eta2sim[j] <- eta2sim[j-1]
    eta3sim[j] <- eta3sim[j-1]
    reject[j] <- 1
  }
}
1-sum(reject)/T #Acceptation rate
#head(cbind(eta1sim,eta2sim,eta3sim), 10)
par(mfrow=c(3,2))
ts.plot(eta1sim[1:200000], ylab = ~log(alpha));abline(h=c(log(1),mean(eta1sim[1:200000])),col=c("green","darkred"));acf(eta1sim[1:200000])
ts.plot(eta2sim[1:200000], ylab = ~log(lambda));abline(h=c(log(3/2),mean(eta2sim[1:200000])),col=c("green","darkred"));acf(eta2sim[1:200000])
ts.plot(eta3sim[1:200000], ylab = ~log(phi));abline(h=c(log(3),mean(eta3sim[1:200000])),col=c("green","darkred"));acf(eta3sim[1:200000])
par(mfrow=c(1,3))
plot(eta1sim,eta2sim, xlab = ~log(alpha), ylab = ~log(lambda))
plot(eta1sim,eta3sim, xlab = ~log(alpha), ylab = ~log(phi))
plot(eta2sim,eta3sim, xlab = ~log(lambda), ylab = ~log(phi))


#################################
#########  Diagnostic  #########
###############################

## Resumen
summary(muestra)
## Split-R
#alpha, lambda, phi
rhat_basic(muestra[,1], split=T);rhat_basic(muestra[,2], split=T);rhat_basic(muestra[,3], split=T)

## R-Vehtari (máximo entre Rbulk y Rpleg)
#alpha, lambda, phi
rhat(muestra[,1]);rhat(muestra[,2]);rhat(muestra[,3])

## Neff Gelman
#alpha, lambda, phi
ess_basic(muestra[,1], split=T);ess_basic(muestra[,2], split=T);ess_basic(muestra[,3], split=T)

## Neff Bulk
#alpha, lambda, phi
ess_bulk(muestra[,1]);ess_bulk(muestra[,2]);ess_bulk(muestra[,3])

## Neff Tail
#alpha, lambda, phi
ess_tail(muestra[,1]);ess_tail(muestra[,2]);ess_tail(muestra[,3])

## Neff Median
#alpha, lambda, phi
ess_quantile(muestra[,1],probs=0.5);ess_quantile(muestra[,2],probs=0.5);ess_quantile(muestra[,3],probs=0.5)

## Neff MAD
#alpha, lambda, phi
ess_sd(muestra[,1]);ess_sd(muestra[,2]);ess_sd(muestra[,3])


#################################
#########  Appendix  ###########
###############################

#####Grafics
#log-posterior and gradient with original parameters
graph1 = graph2 = graph3 = numeric()
x = seq(0.001,10, length = 1000)
for(i in 1:1000){
  graph1[i] = logpost(c(x[i],2,3))
  graph2[i] = logpost(c(1,x[i],3))
  graph3[i] = logpost(c(1,2,x[i]))
}
par(mfrow=c(2,3))
plot(x, exp(graph1),type="l",lwd=2,col="darkred",xlim=c(0,4),main=~alpha)
plot(x, exp(graph2),type="l",lwd=2,col="darkred",xlim=c(0,4),main=~mu)
plot(x, exp(graph3),type="l",lwd=2,col="darkred",xlim=c(0,4),main=~phi)
for(i in 1:1000){
  graph1[i] = grad(logpost,c(x[i],2,3))[1]
  graph2[i] = grad(logpost,c(1,x[i],3))[2]
  graph3[i] = grad(logpost,c(1,2,x[i]))[3]
}
plot(x,graph1,type="l",lwd=2,col="darkred",xlim=c(0,4),ylim=c(-200,200),main=~alpha);abline(h=0)
plot(x,graph2,type="l",lwd=2,col="darkred",xlim=c(0,4),ylim=c(-200,200),main=~mu);abline(h=0)
plot(x,graph3,type="l",lwd=2,col="darkred",xlim=c(0,4),ylim=c(-200,200),main=~phi);abline(h=0)


## log-posterior and gradient of reparametrization
for(i in 1:1000){
  graph1[i] = logpost.rep(c(x[i],3/2,3))
  graph2[i] = logpost.rep(c(1,x[i],3))
  graph3[i] = logpost.rep(c(1,3/2,x[i]))
}
par(mfrow=c(2,3))
plot(x,exp(graph1),type="l",lwd=2,col="darkred",xlim=c(0,5),main=~alpha)
plot(x,exp(graph2),type="l",lwd=2,col="darkred",xlim=c(0,5),main=~lambda)
plot(x,exp(graph3),type="l",lwd=2,col="darkred",xlim=c(0,5),main=~phi)
for(i in 1:1000){
  graph1[i] = grad(logpost.rep,c(x[i],3/2,3))[1]
  graph2[i] = grad(logpost.rep,c(1,x[i],3))[2]
  graph3[i] = grad(logpost.rep,c(1,3/2,x[i]))[3]
}
plot(x,graph1,type="l",lwd=2,col="darkred",xlim=c(0,5),ylim=c(-200,200),main=~alpha);abline(h=0)
plot(x,graph2,type="l",lwd=2,col="darkred",xlim=c(0,5),ylim=c(-200,200),main=~lambda);abline(h=0)
plot(x,graph3,type="l",lwd=2,col="darkred",xlim=c(0,5),ylim=c(-200,200),main=~phi);abline(h=0)

## log-posterior and gradient of log-transform
for(i in 1:1000){
  graph1[i] = logpost.eta(c(log(x[i]),log(3/2),log(3)))
  graph2[i] = logpost.eta(c(log(1),log(x[i]),log(3)))
  graph3[i] = logpost.eta(c(log(1),log(3/2),log(x[i])))
}
par(mfrow=c(2,3))
plot(log(x),exp(graph1),type="l",lwd=2,col="darkred",xlim=c(-0.5,1.5),main=~log(alpha))
plot(log(x),exp(graph2),type="l",lwd=2,col="darkred",xlim=c(-0.5,1.5),main=~log(lambda))
plot(log(x),exp(graph3),type="l",lwd=2,col="darkred",xlim=c(-0.5,1.5),main=~log(phi))
for(i in 1:1000){
  graph1[i] = grad(logpost.eta,c(log(x[i]),log(3/2),log(3)))[1]
  graph2[i] = grad(logpost.eta,c(log(1),log(x[i]),log(3)))[2]
  graph3[i] = grad(logpost.eta,c(log(1),log(3/2),log(x[i])))[3]
}
plot(log(x),graph1,type="l",lwd=2,col="darkred",xlim=c(-0.5,1.5),ylim=c(-200,300),main=~log(alpha));abline(h=0)
plot(log(x),graph2,type="l",lwd=2,col="darkred",xlim=c(-0.5,1.5),ylim=c(-200,300),main=~log(lambda));abline(h=0)
plot(log(x),graph3,type="l",lwd=2,col="darkred",xlim=c(-0.5,1.5),ylim=c(-200,300),main=~log(phi));abline(h=0)


### FIN.