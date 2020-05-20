##########               Real Data Results               ##########
library("geoR")
library("MaxPro")
library("mined")
library("MASS")
library("fields")
library("np")
library("doParallel")
library("foreach")

##########    Real Data 1. North American Rainfall (parameter estimation and prediction)   ##########
data("NorthAmericanRainfall")
loc=cbind(NorthAmericanRainfall$longitude,NorthAmericanRainfall$latitude)
y=NorthAmericanRainfall$precip
Sloc=apply(loc,2,function(x) (x-min(x))/max(x-min(x)))
Data=data.frame(loc.x1=Sloc[,1],loc.x2=Sloc[,2],y=y)
bw=npregbw(Data$y~Data$loc.x1+Data$loc.x2,regtype="ll",bwmethod="cv.aic")
regres=npreg(bw,gradients = TRUE)
grad=regres$grad
cgrad=sqrt(grad[,1]^2+grad[,2]^2)


PE_PA_nonug=function(daTa,test,tgrad,tloc,n)  #data=train
{
  res=matrix(0,5,7)
  findi=function(x){idex=which(daTa$loc.x1==x[1] & daTa$loc.x2==x[2])}
  
  #MED1
  MEDres=mined::SelectMinED(cbind(daTa$loc.x1,daTa$loc.x2),log(tgrad),n,1,2)
  MEDloc=MEDres$points
  MEDi=apply(MEDloc,1,findi)
  MEDloc1=tloc[MEDi,]
  MEDy=daTa$y[MEDi]
  MEDGPres=likfit(as.geodata(cbind(MEDloc1,MEDy)),trend = "1st",ini.cov.pars = c(1e6,10),fix.nugget = TRUE)
  res[1,-7]=c(MEDGPres$beta,MEDGPres$cov.pars,MEDGPres$phi/MEDGPres$sigmasq)
  MEDprid=krige.conv(as.geodata(cbind(MEDloc1,MEDy)),locations = test[,-3],
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = MEDGPres$beta,
                                         cov.pars = MEDGPres$cov.pars))
  MEDpred=MEDprid$predict
  res[1,7]=sqrt(sum((MEDpred-test[,3])^2)/720)
  
  #MaxPro
  iexist=sample(1:nrow(daTa),1)
  exist=daTa[iexist,1:2]
  MaxProloc=MaxProAugment(exist,daTa[,1:2],nNew=n-1)$Design
  MaxProi=apply(MaxProloc,1,findi)
  MaxProloc1=tloc[MaxProi,]
  MaxProy=daTa$y[MaxProi]
  MAPGPres=likfit(as.geodata(cbind(MaxProloc1,MaxProy)),trend = "1st",ini.cov.pars = c(1e6,10),fix.nugget = TRUE)     #ini value
  res[2,-7]=c(MAPGPres$beta,MAPGPres$cov.pars,MAPGPres$phi/MAPGPres$sigmasq)
  MAPprid=krige.conv(as.geodata(cbind(MaxProloc1,MaxProy)),locations = test[,-3],
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = MAPGPres$beta,
                                         cov.pars = MAPGPres$cov.pars))
  MAPpred=MAPprid$predict
  res[2,7]=sqrt(sum((MAPpred-test[,3])^2)/720)
  
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
  DaWloc1=tloc[inr,]
  DaWy=daTa$y[inr]
  DaWGPres=likfit(as.geodata(cbind(DaWloc1,DaWy)),trend = "1st",ini.cov.pars = c(1e6,10),fix.nugget = TRUE)   
  res[3,-7]=c(DaWGPres$beta,DaWGPres$cov.pars,DaWGPres$phi/DaWGPres$sigmasq)
  DaWprid=krige.conv(as.geodata(cbind(DaWloc1,DaWy)),locations = test[,-3],
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = DaWGPres$beta,
                                         cov.pars = DaWGPres$cov.pars))
  DaWpred=DaWprid$predict
  res[3,7]=sqrt(sum((DaWpred-test[,3])^2)/720)
  
  #Random
  ir=sample(1:nrow(daTa),n)
  Rloc=tloc[ir,]
  Ry=daTa$y[ir]
  RGPres=likfit(as.geodata(cbind(Rloc,Ry)),trend = "1st",ini.cov.pars = c(1e6,10),fix.nugget = TRUE)     #ini value
  res[4,-7]=c(RGPres$beta,RGPres$cov.pars,RGPres$phi/RGPres$sigmasq)
  Rprid=krige.conv(as.geodata(cbind(Rloc,Ry)),locations = test[,-3],
                   krige=krige.control(trend.d = "1st",trend.l = "1st",beta = RGPres$beta,
                                       cov.pars = RGPres$cov.pars))
  Rpred=Rprid$predict
  res[4,7]=sqrt(sum((Rpred-test[,3])^2)/720)
  
  #MED4
  MEDres4=mined::SelectMinED(cbind(daTa$loc.x1,daTa$loc.x2),4*log(tgrad),n,1,2)
  MEDloc4=MEDres4$points
  MEDi4=apply(MEDloc4,1,findi)
  MEDloc14=tloc[MEDi4,]
  MEDy4=daTa$y[MEDi4]
  MEDGPres4=likfit(as.geodata(cbind(MEDloc14,MEDy4)),trend = "1st",ini.cov.pars = c(1e6,10),fix.nugget = TRUE)
  res[5,-7]=c(MEDGPres4$beta,MEDGPres4$cov.pars,MEDGPres4$phi/MEDGPres4$sigmasq)
  MEDprid4=krige.conv(as.geodata(cbind(MEDloc14,MEDy4)),locations = test[,-3],
                      krige=krige.control(trend.d = "1st",trend.l = "1st",beta = MEDGPres4$beta,
                                          cov.pars = MEDGPres4$cov.pars))
  MEDpred4=MEDprid4$predict
  res[5,7]=sqrt(sum((MEDpred4-test[,3])^2)/720)
  
  return(res)
}


PEPAnonug4_main=function()
{
  nseq=c(50,100,200,300,400,500)
  it=sample(1:1720,1000)
  train=Data[it,]
  tgrad=cgrad[it]
  tloc=loc[it,]
  test=cbind(loc,y)[-it,]
  
  MLEres=likfit(as.geodata(cbind(loc,y)[it,]),trend = "1st",ini.cov.pars = c(1e6,10),fix.nugget = TRUE)
  MLEpara=c(MLEres$beta,MLEres$cov.pars,MLEres$phi/MLEres$sigmasq)
  predres=krige.conv(as.geodata(cbind(loc,y)[it,]),locations =test[,-3],
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = MLEres$beta,
                                         cov.pars = MLEres$cov.pars))
  pred=sqrt(sum((predres$predict-test[,3])^2)/720)
  MlE=c(MLEpara,pred)
  
  PE_PA1=function(n){return(PE_PA_nonug(train,test,tgrad,tloc,n))}
  out=sapply(nseq,PE_PA1,simplify = FALSE)
  
  res=list(MLE=MlE,subsampling=out)
  return(res)
}

cl<- makeCluster(2) 
registerDoParallel(cl) 
ResultN = foreach(i=1:500,
                  .combine=cbind,
                  .packages=c("geoR","mined","np","MaxPro","fields","MASS")) %dopar% PEPAnonug4_main()
stopCluster(cl)
