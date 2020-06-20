library(np)
library(mined)
library(randtoolbox)
library(nabor)
library(lhs)
library(geoR)
library(gss)
library(foreach)
library(doParallel)
library(ggplot2)
#source("ssa.R")
#source("esps-m.R")
#source("ABS_select.R")

##########    Real Data 2. Ozone    ##########
# Raw data and plot the figure
df.zoom=Ozone_dat[,c(2:3,6)]
fig.raw<-ggplot(df.zoom, aes(lon,lat,colour=ozone))
fig.raw+geom_point()+
  scale_color_gradient(low = "black", high = "white")+
  ggtitle("Raw image")
names(df.zoom)<-c("x","y","value")

Data=data.frame(loc.x1=df.zoom$lon,loc.x2=df.zoom$lat,y=df.zoom$ozone)
ib=sample(1:nrow(Data),1000)
Dataib=Data[ib,]
bw=npregbw(Dataib$y~Dataib$loc.x1+Dataib$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
regres=npreg(Data$y~Data$loc.x1+Data$loc.x2,bws=bw,gradients=TRUE)
grad=regres$grad
cgrad=sqrt(grad[,1]^2+grad[,2]^2)
Data$grad=cgrad
nseq=round(c(5,10,15,20)*(15e4)^(2/9))
#  71 141 212 283

##########    Universal Kriging Model and Smoothing Spline Model   ##########
PAO=function(daTa,test,n,lambda.pro,theta,seed)
{
  res=matrix(0,4,2)     
  findi=function(x){idex=which(daTa$loc.x1==x[1] & daTa$loc.x2==x[2])}
  
  #GBMED 
  MEDres=mined::SelectMinED(cbind(daTa$loc.x1,daTa$loc.x2),log(daTa$grad),n,1,2)
  MEDloc=MEDres$points
  MEDi=apply(MEDloc,1,findi)
  MEDy=daTa$y[MEDi]
  
  MEDGPres=likfit(as.geodata(cbind(MEDloc,MEDy)),trend ="1st",ini.cov.pars = c(5e3,50),fix.nugget = TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))  
  MEDpred=krige.conv(as.geodata(cbind(MEDloc,MEDy)),locations = cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d ="1st",trend.l = "1st",beta = MEDGPres$beta,
                                         cov.pars = MEDGPres$cov.pars))
  MEDyhat=MEDpred$predict
  res[1,1]=log(mean((test$y-MEDyhat)^2))
  
  MEDfit=ssa(y~loc.x1*loc.x2, data=daTa[,-4], id.basis=MEDi,
             lambda = lambda.pro, theta=theta)
  MEDyhats=predict(MEDfit,newdata=test[,-3], se.fit=F)
  res[1,2]=log(mean((test$y-MEDyhats)^2))
  
  # SBS
  set.seed(seed)
  design=sobol(n,2)
  design[,1]=design[,1]*179.99*2-179.99
  design[,2]=design[,2]*(80.53+87.48)-87.48
  SBSi=nabor::knn(daTa[,1:2],design,k = 1)$nn.idx
  SBS=daTa[SBSi,1:3]
  SBSGPres=likfit(as.geodata(SBS),trend = "1st",ini.cov.pars = c(5e3,50),fix.nugget = TRUE,
         limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))   
  SBSpred=krige.conv(as.geodata(SBS),locations = cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = SBSGPres$beta,
                                         cov.pars = SBSGPres$cov.pars))
  SBSyhat=SBSpred$predict
  res[2,1]=log(mean((test$y-SBSyhat)^2))
  
  SBSfit=ssa(y~loc.x1*loc.x2, data=daTa[,-4], id.basis=SBSi,
             lambda = lambda.pro, theta=theta)
  SBSyhats=predict(SBSfit,newdata=test[,-3], se.fit=F)
  res[2,2]=log(mean((test$y-SBSyhats)^2))
  
  #ABS
  set.seed(seed)
  ABSi=ABS_select(daTa$y,n)
  ABS=daTa[ABSi,1:3]
  ABSGPres=likfit(as.geodata(ABS),trend = "1st",ini.cov.pars = c(5e3,50),fix.nugget = TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))   
  ABSpred=krige.conv(as.geodata(ABS),locations = cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = ABSGPres$beta,
                                         cov.pars = ABSGPres$cov.pars))
  ABSyhat=ABSpred$predict
  res[3,1]=log(mean((test$y-ABSyhat)^2))
  
  ABSfit=ssa(y~loc.x1*loc.x2, data=daTa[,-4], id.basis=ABSi,
             lambda = lambda.pro, theta=theta)
  ABSyhats=predict(ABSfit,newdata=test[,-3], se.fit=F)
  res[3,2]=log(mean((test$y-ABSyhats)^2))
  
  #Uniform
  set.seed(seed)
  UNIi=sample(1:nrow(daTa),n)
  UNI=daTa[UNIi,1:3]
  UNIGPres=likfit(as.geodata(UNI),trend = "1st",ini.cov.pars = c(5e3,50),fix.nugget = TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))   
  UNIpred=krige.conv(as.geodata(UNI),locations = cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = UNIGPres$beta,
                                         cov.pars = UNIGPres$cov.pars))
  UNIyhat=UNIpred$predict
  res[4,1]=log(mean((test$y-UNIyhat)^2))
  
  UNIfit=ssa(y~loc.x1*loc.x2, data=daTa[,-4], id.basis=UNIi,
             lambda = lambda.pro, theta=theta)
  UNIyhats=predict(UNIfit,newdata=test[,-3], se.fit=F)
  res[4,2]=log(mean((test$y-UNIyhats)^2))
  
  return(res)
}


PAO_main=function(nseq,seed)
{
  set.seed(seed)
  it=sample(1:nrow(Data),15e4)
  Dtrain=Data[it,]
  Dtest=Data[-it,-4]

  sam.size=ceiling(50*((15e4)^{1/4}))
  esps.re=esps(y~loc.x1*loc.x2, data=Dtrain[,-4], sam.size,r=3)
  r=3
  p=esps.re$p
  pro.gcv.lamb=esps.re$lambda*((15e4/sam.size)^(-r/(p*r+1)))
  lambda.pro=log10(15e4*pro.gcv.lamb)
  theta=esps.re$theta
  
  PAO1=function(n){return(PAO(Dtrain,Dtest,n,lambda.pro,theta,seed))}
  out=sapply(nseq,PAO1,simplify = FALSE)
  return(out)
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultO= foreach(i=1:200,
                  .combine=cbind,
                  .packages=c("geoR","mined","gss","nabor","randtoolbox")) %dopar% PAO_main(nseq,666*i+888)
stopCluster(cl)
