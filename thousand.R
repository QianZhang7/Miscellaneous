library(nlme)
library(MASS)
library(boot)
library(Matrix)
library(BRugs)
library(lattice)
library(lme4)
library(R2WinBUGS)
library(coda)
library(lattice)
 
H = 1000
NP = 100
A=63
TIMES = 100
it=5000
# Max Number of Clusters.
K=30
M=5500
N=500
lowq = 0.05
uppq=0.95





rbeta<-array(NA,dim=c(it,1,TIMES))
rclass<-array(NA,dim=c(it,6,TIMES))
rcluster<-array(NA,dim=c(it,1,TIMES))
rclusterhop<-array(NA,dim=c(it,H,TIMES))
renum<-array(NA,dim=c(it,H,TIMES))
rgamma<-array(NA,dim=c(it,H,TIMES))
rintercept<-array(NA,dim=c(it,1,TIMES))
rmedianhop<-array(NA,dim=c(H,1,TIMES))
rpk<-array(NA,dim=c(it,K,TIMES))
rpkhop<-array(NA,dim=c(it,H,TIMES))
rqhop<-array(NA,dim=c(it,H,TIMES))
rRand<-array(NA,dim=c(H,1,TIMES))
rrank<-array(NA,dim=c(it,H,TIMES))
rrankhop<-array(NA,dim=c(it,H,TIMES))
rtheta<-array(NA,dim=c(it,K,TIMES))
rthetahop<-array(NA,dim=c(it,H,TIMES))
rthetahopDec<-array(NA,dim=c(it,H,TIMES))
rZ<-array(NA,dim=c(it,H,TIMES))
rtheDec<-array(NA,dim=c(it,H,TIMES))
rtheclass<-array(NA,dim=c(H,6,TIMES))
rRandDec<-array(NA,dim=c(H,1,TIMES))

for(k in 1:TIMES)
{
set.seed(A+k)

hno = sort(rep(1:H,NP))
NT = H*NP

age = rnorm(NT)
beta0=-4
beta1=1.2
Rand = rgamma(H, 2,1)
pij = numeric(NT)
for(i in 1:H)    pij[hno == i] = beta0+beta1*age[hno == i] + Rand[i]
# To construct the response, and let \beta = 1.
P = exp(pij)/(1+exp(pij))
death = numeric(NT)                       
U = runif(NT)
for(i in 1:NT) death[i] = ifelse (U[i]<P[i],1,0)

exp1 = data.frame(death,age,hno)
colnames(exp1) = c("death","age","hno")


#initial values
gamma <- rep(0,H)
theta <- rep(0,K)
tauinv <- rep(1,K)
    r <- rep(0.05,K-1)
    Z <- rep(1,H)
n<-c(NP)




################################################################################################
### The following R are for preprocessing the data before calling the WinBUGS
### If you are not interested in R, you just copy them into R without understanding them
Tran2HopbyPat <- function(data,Hid,Var){ ### The function to reshape the data
   myhid<-data[,names(data)==Hid]
   myvar<-data[,names(data)==Var]
   N <- length(myhid)  ### Total # of hispitals
   H <- length(table(myhid))   ### The MAX patients # of hospitals 
   MaxPat=max(table(myhid))
   hospital=rep(0,H)
   num.patient=rep(0,H)
   conv=array(NA,dim=c(H,MaxPat))
   conv[1,1]=myvar[1]
   cnt=1
   temp=1
   hospital[1]=myhid[1]
   for (i in 2:N){
	temp=temp+1;
	if (!(myhid[i-1] == myhid[i])){
	num.patient[cnt]=temp-1	
	cnt=cnt+1
	hospital[cnt]=myhid[i]
	temp=1
	}
	conv[cnt,temp]=myvar[i]
  }
   return(conv)
 }  ###End of function "Tran2HopbyPat" definition
death<-Tran2HopbyPat(exp1,"hno","death")
age<-Tran2HopbyPat(exp1,"hno","age") ### number # of hospitals
H <- length(table(exp1$hno)) ### patient number # of each hospital
hnof<-table(exp1$hno)
n<-hnof[]







setwd("D:\\D") ### Change to your path containing the data and Lab4BugsModel.bug. 


################################################################################################
## This is the part of R2WinBUGS to call WinBUGS from R
################################################################################################
###Prepare data and parameters' initial values for R to call WinBUGS
datainput <- list(death=death,age=age,H=H,K=K,n=as.numeric(n))
inits <- function(){list(beta=0,intercept=0,gamma=gamma,theta=theta,tauinv=tauinv,Z=Z,r=r,basemu=0, sigmaF0=1, DPPalpha=1,b=1)}
###Call WinBUGS using the bugs() function in the R2WinBUGS package
sim = bugs(datainput,inits,
           model.file="Lab4RankHospRandom.bug",
           parameters.to.save=c("beta","intercept","gamma", "theta", "Cluster","emat","enum","rank","pk","Z"),
           n.chains=1, n.iter=M, n.burnin=N, n.thin=1,
           bugs.directory="D:\\WinBUGS14",
           codaPkg=FALSE,debug=FALSE,bugs.seed=21538765)
### After the WinBUGS indicates that the 1000 updates were done, YOU MAY NEED TO MANUALLY CLOSE THE WinBUGS WINDOW


rank=sim$sims.list$rank
 fn1=paste("D:\\D\\rank",k,".csv")
  write.csv(rank,file=fn1)
rrank[,,k]=rank


Z=sim$sims.list$Z
 fn2=paste("D:\\D\\Z",k,".csv")
  write.csv(Z,file=fn2)
rZ[,,k]=Z

beta=sim$sims.list$beta
 fn3=paste("D:\\D\\beta",k,".csv")
  write.csv(beta,file=fn3)
rbeta[,,k]=beta

 theta=sim$sims.list$theta
 fn4=paste("D:\\D\\theta",k,".csv")
 write.csv(theta,file=fn4)
 rtheta[,,k]=theta
 
 
 pk=sim$sims.list$pk
 fn5=paste("D:\\D\\pk",k,".csv")
  write.csv(pk,file=fn5)
rpk[,,k]=pk
  
  
  intercept=sim$sims.list$intercept
 fn6=paste("D:\\D\\intercept",k,".csv")
  write.csv(intercept,file=fn6)
  rintercept[,,k]=intercept

  Cluster=sim$sims.list$Cluster
 fn7=paste("D:\\D\\cluster",k,".csv")
  write.csv(Cluster,file=fn7)
   rcluster[,,k]=Cluster 
  
  
  emat=sim$sims.list$emat
 fn8=paste("D:\\D\\emat",k,".csv")
  write.csv(emat,file=fn8)

  
  
  
  enum=sim$sims.list$enum
 fn9=paste("D:\\D\\enum",k,".csv")
  write.csv(enum,file=fn9)
    renum[,,k]=enum
  

  gamma=sim$sims.list$gamma
 fn10=paste("D:\\D\\gamma",k,".csv")
  write.csv(gamma,file=fn10)
  rgamma[,,k]=gamma

 fn11=paste("D:\\D\\Rand",k,".csv")
  write.csv(Rand,file=fn11)

  rRand[,,k]=Rand

# find theta and pk for every hospital
thetahop <- rep(0,it*H)
dim(thetahop) <- c(it,H)

for(i in 1:it){ 
  for (j in 1:H){
n<-Z[i,j]
    thetahop[i,j]<-theta[i,n]
	}
}
  rthetahop[,,k]=thetahop
  
  
pkhop <- rep(0,it*H)
dim(pkhop) <- c(it,H)
for(i in 1:it){ 
  for (j in 1:H){
n<-Z[i,j]
    pkhop[i,j]<-pk[i,n]
	}
}
  rpkhop[,,k]=pkhop

#rank the hopital by theta
clusterhop <- rep(0,it*H)
dim(clusterhop) <- c(it,H)
for(i in 1:it){ 
clusterhop[i,]<-as.numeric(factor(thetahop[i,]))
}
  rclusterhop[,,k]=clusterhop

#qualify the ranks

qhop <- rep(0,it*H)
dim(qhop) <- c(it,H)
for(i in 1:it){ 
qhop[i,]<-(clusterhop[i,]+1)/(Cluster[i]+1)

}
  rqhop[,,k]=qhop

#rank hospital using median of qhop
medianhop=matrix(nrow=H,ncol=1)
for(i in 1:H)
{
medianhop[i]=quantile(qhop[,i],0.5)
}
  rmedianhop[,,k]=medianhop


rankhop<-rank(medianhop,ties.method= "min")

 fn12=paste("D:\\D\\thetahop",k,".csv")
 
  write.csv(thetahop,file=fn12)
  
   fn13=paste("D:\\D\\pkhop",k,".csv")

  write.csv(pkhop,file=fn13)
  
    fn14=paste("D:\\D\\clusterhop",k,".csv")

  write.csv(clusterhop,file=fn14)
   
     fn15=paste("D:\\D\\qhop",k,".csv")

  write.csv(qhop,file=fn15)
  
     fn16=paste("D:\\D\\medianhop",k,".csv")

  write.csv(medianhop,file=fn16)
  
     fn17=paste("D:\\D\\rankhop",k,".csv")

  write.csv(rankhop,file=fn17)
  

  



  mthetahop = numeric(H)
  lthetahop = numeric(H)
  uthetahop = numeric(H)
thetahopDec = character(it*H)
dim(thetahopDec)=c(it,H)
for(i in 1:it) 
{

   lthetahop[i]= quantile(thetahop[i,],lowq)
  uthetahop[i]= quantile(thetahop[i,],uppq)

}

for(i in 1:it){
for(j in 1:H)
{
thetahopDec[i,j] = ifelse(lthetahop[i]>thetahop[i,j],"Good", ifelse(uthetahop[i]<thetahop[i,j],"Bad","Normal"))

}
}

  fn18=paste("D:\\D\\thetahopDec",k,".csv")
  write.csv(thetahopDec,file=fn18)



class = matrix(nrow=H,ncol=6)
colnames(class) = c("Good","Normal","Bad","Good%","Normal%","Bad%")
for(i in 1:H)
{
class[i,1] = length(which(as.character(thetahopDec[,i])=="Good" ))
class[i,2] = length(which(as.character(thetahopDec[,i])=="Normal" ))
class[i,3] = length(which(as.character(thetahopDec[,i])=="Bad" ))
class[i,4]=class[i,1]/it*H
class[i,5]=class[i,2]/it*H
class[i,6]=class[i,3]/it*H
}

  fn19=paste("D:\\D\\class",k,".csv")
  write.csv(class,file=fn19)

# summarize theta in every loop
mthe = numeric(k)
lthe = numeric(k)
uthe = numeric(k)
mthe[k] = mean(rthetahop[,,k])
lthe[k] = quantile(rthetahop[,,k],lowq)
uthe[k] = quantile(rthetahop[,,k],uppq)


theDec = character(it*H)
dim(theDec)=c(it,H)

for(i in 1:it){
for(j in 1:H)
{
theDec[i,j] = ifelse(lthe[k]>rthetahop[i,j,k],"Good", ifelse(uthe[k]<rthetahop[i,j,k],"Bad","Normal"))

}
}

  fn20=paste("D:\\D\\theDec",k,".csv")
  write.csv(theDec,file=fn20)
 rtheDec[,,k]=theDec


theclass = matrix(nrow=H,ncol=6)
colnames(theclass) = c("Good","Normal","Bad","Good%","Normal%","Bad%")
for(i in 1:H)
{
theclass[i,1] = length(which(as.character(thetahopDec[,i])=="Good" ))
theclass[i,2] = length(which(as.character(thetahopDec[,i])=="Normal" ))
theclass[i,3] = length(which(as.character(thetahopDec[,i])=="Bad" ))
theclass[i,4]=class[i,1]/it*H
theclass[i,5]=class[i,2]/it*H
theclass[i,6]=class[i,3]/it*H
}

  fn21=paste("D:\\D\\theclass",k,".csv")
  write.csv(theclass,file=fn21)
 rtheclass[,,k]=theclass

#classification all loops
mRand = numeric(k)
lRand = numeric(k)
uRand = numeric(k)
mRand[k] = mean(rRand[,,k])
lRand[k] = quantile(rRand[,,k],lowq)
uRand[k] = quantile(rRand[,,k],uppq)

RandDec = character(H)
for(j in 1:H)
{
RandDec[j] = ifelse(lRand[k]>rRand[j,1,k],"Good", ifelse(uRand[k]<rRand[j,1,k],"Bad","Normal"))
}

fn22=paste("D:\\D\\RandDec",k,".csv")
 write.csv(RandDec,file=fn22)
 rRandDec[,,k]=RandDec
}



