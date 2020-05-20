##########               Simulation Results               ##########
library("geoR")
library("MaxPro")
library("mined")
library("MASS")
library("fields")
library("np")
library("doParallel")
library("foreach")

##########    Silmulation 1. Parameter Estimation Without Nugget Effect   ##########
PE4_nonug=function(daTa,cgrad,n)
{
  para=matrix(0,5,6)     # without nugget effect
  findi=function(x){idex=which(daTa$loc.x1==x[1] & daTa$loc.x2==x[2])}
  
  #MED1: MED we use 
  MEDres=mined::SelectMinED(cbind(daTa$loc.x1,daTa$loc.x2),log(cgrad),n,1,2)
  MEDloc=MEDres$points
  MEDi=apply(MEDloc,1,findi)
  MEDy=daTa$y[MEDi]
  MEDGPres=likfit(as.geodata(cbind(MEDloc,MEDy)),trend = "1st",ini.cov.pars = c(1,1),fix.nugget = TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))  
  para[1,]=c(MEDGPres$beta,MEDGPres$cov.pars,MEDGPres$phi/MEDGPres$sigmasq)
  
  #MaxPro
  iexist=sample(1:nrow(daTa),1)
  exist=daTa[iexist,1:2]
  MaxProloc=MaxProAugment(exist,daTa[,1:2],nNew=n-1)$Design
  MaxProi=apply(MaxProloc,1,findi)
  MaxProy=daTa$y[MaxProi]
  MAPGPres=likfit(as.geodata(cbind(MaxProloc,MaxProy)),trend = "1st",ini.cov.pars = c(1,1),fix.nugget = TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf))) 
  para[2,]=c(MAPGPres$beta,MAPGPres$cov.pars,MAPGPres$phi/MAPGPres$sigmasq)
  
  #Deep and Wide (one time)
  ic=sample(1:nrow(daTa),5)
  Dis=rdist(daTa[,-3])
  inb=matrix(0,5,round(n/5))
  for(i in 1:5)
  {
    rk=rank(Dis[ic[i],],ties.method = "first")
    inb[i,]=which(rk>1 & rk<=round(n/5)+1)
  }
  inr=unique(as.vector(inb))
  DaW=daTa[inr,]
  DaWGPres=likfit(as.geodata(DaW),trend = "1st",ini.cov.pars = c(1,1),fix.nugget=TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))   
  para[3,]=c(DaWGPres$beta,DaWGPres$cov.pars,DaWGPres$phi/DaWGPres$sigmasq)
  
  #Random
  ir=sample(1:nrow(daTa),n)
  R=daTa[ir,]
  RGPres=likfit(as.geodata(R),trend = "1st",ini.cov.pars = c(1,1),fix.nugget = TRUE,
                limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))   
  para[4,]=c(RGPres$beta,RGPres$cov.pars,RGPres$phi/RGPres$sigmasq)
  
  #MED4: MED emphasize gradient
  MEDres4=mined::SelectMinED(cbind(daTa$loc.x1,daTa$loc.x2),4*log(cgrad),n,1,2)
  MEDloc4=MEDres4$points
  MEDi4=apply(MEDloc4,1,findi)
  MEDy4=daTa$y[MEDi4]
  MEDGPres4=likfit(as.geodata(cbind(MEDloc4,MEDy4)),trend = "1st",ini.cov.pars = c(1,1),fix.nugget = TRUE,
                   limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))  
  para[5,]=c(MEDGPres4$beta,MEDGPres4$cov.pars,MEDGPres4$phi/MEDGPres4$sigmasq)
  
  return(para)
}

PEnonug4_main=function()
{
  m=1e4
  beta=c(1,2,3)
  sigma2=1
  phi=0.25
  nseq=c(30,50,100,150,200)
  ### generate the data
  geoData=grf(m,cov.pars = c(sigma2,phi))
  y=beta[1]+beta[2]*geoData$coords[,1]+beta[3]*geoData$coords[,2]+geoData$data
  Data=data.frame(loc.x1=geoData$coords[,1],loc.x2=geoData$coords[,2],y=y)
  
  ### gradient evaluation
  ib=sample(1:m,1e3)
  Dataib=Data[ib,]
  bw=npregbw(Dataib$y~Dataib$loc.x1+Dataib$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  regres=npreg(Data$y~Data$loc.x1+Data$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  
  ### parameter estimation
  PE1=function(n){return(PE4_nonug(Data,cgrad,n))}
  out=sapply(nseq,PE1,simplify = FALSE)
  return(out)
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
Result1= foreach(i=1:500,
                 .combine=cbind,
                 .packages=c("geoR","mined","np","MaxPro","fields","MASS")) %dopar% PEnonug4_main()
stopCluster(cl)

