p.pois <- function(K, Ne, Eb){ #Model with Poisson distirbuted effective parasite number 
  pc=10^(-15)
  ##N is number of sporozoites
  t=48 ##hours
  pn= exp(-1*Ne*exp(-1*K*Eb*t))
  for(i in 1:length(pn)){
    if(pn[i]<=pc){
      pn[i] = pc
    }
    if(pn[i]>1-pc){
      pn[i] = 1-pc
    }}
  return(pn)
}

L.pois <- function(K, Ne, data){
  Eb=data$PBL
  pinf=data$prot   #prob of infection
  
  neg.log.lik=-sum(log((p.pois(K,Ne,Eb)^(pinf))*((1-p.pois(K,Ne,Eb))^(1-pinf))))
  
  return(neg.log.lik)
}

library(bbmle)

fit.pb = mle2(minuslogl=L.pois, start = list(K=0.1,Ne=5), 
                data=list(data=pb),
                control = list(ndeps=c(1e-7,1e-7),maxit=10000)) ## two parameter mle using maxit = 10000 for complete convergence

summary(fit.pb)

pb0 <- profile(fit.pb)

confint(pb0)

confint(pb0,method="quad")

confint(pb0,method="uniroot")

par(mfrow=c(1,2))
plot(pb0,plot.confstr=TRUE)



############### model for py numsp=10 ########################
fit.py10 = mle2(minuslogl=L.pois, start = list(K=0.1,Ne=0.1), 
                 data=list(data=py10),
                 control = list(ndeps=c(1e-7,1e-7),maxit=100000))




summary(fit.py10)

fpy10 <- profile(fit.py10)

confint(fpy10)

par(mfrow=c(1,2))
plot(fpy10,plot.confstr=TRUE)

############### model for py numsp=50 ########################
fit.py50 = mle2(minuslogl=L.pois, start = list(K=0.1,Ne=5), 
              data=list(data=py50),
              control = list(ndeps=c(1e-7,1e-7),maxit=100000))




summary(fit.py50)

fpy50 <- profile(fit.py50)

confint(fpy50)

par(mfrow=c(1,2))
plot(fpy50,plot.confstr=TRUE)



############### model for py numsp=100 ########################
fit.py100 = mle2(minuslogl=L.pois, start = list(K=0.1,Ne=5), 
                data=list(data=py100),
                control = list(ndeps=c(1e-7,1e-7),maxit=100000))




summary(fit.py100)

fpy100 <- profile(fit.py100)

confint(fpy100)

par(mfrow=c(1,2))
plot(fpy100,plot.confstr=TRUE)




############### model for py numsp=1000 ########################
fit.py1000 = mle2(minuslogl=L.pois, start = list(K=0.1,Ne=3), 
                 data=list(data=py1000),
                 control = list(ndeps=c(1e-7,1e-7),maxit=100000))




summary(fit.py1000)

fpy1000 <- profile(fit.py1000)

confint(fpy1000)

par(mfrow=c(1,2))
plot(fpy1000,plot.confstr=TRUE)




############### model for py numsp=10000 ########################
fit.py10000 = mle2(minuslogl=L.pois, start = list(K=0.1,Ne=3), 
                  data=list(data=py10000),
                  control = list(ndeps=c(1e-7,1e-7),maxit=100000))




summary(fit.py10000)

fpy10000 <- profile(fit.py10000)

confint(fpy10000)

par(mfrow=c(1,2))
plot(fpy10000,plot.confstr=TRUE)


############### model for py entire ########################
fit.py = mle2(minuslogl=L.pois, start = list(K=0.1,Ne=3), 
                   data=list(data=py),
                   control = list(ndeps=c(1e-7,1e-7),maxit=100000))




summary(fit.py)

fpy <- profile(fit.py)

confint(fpy)

par(mfrow=c(1,2))
plot(fpy,plot.confstr=TRUE)









