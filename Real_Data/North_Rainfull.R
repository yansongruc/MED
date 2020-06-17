library("geoR")
library("MaxPro")
library("mined")
library("MASS")
library("fields")
library("np")
library("doParallel")
library("foreach")
library("gss")
#source("ssa.R")
#source("esps-m.R")
#source("ABS_select.R")

##########    Real Data 1. North American Rainfall    ##########
data("NorthAmericanRainfall")
loc=cbind(NorthAmericanRainfall$longitude,NorthAmericanRainfall$latitude)
y=NorthAmericanRainfall$precip
Data=data.frame(loc.x1=loc[,1],loc.x2=loc[,2],y=y)
bw=npregbw(Data$y~Data$loc.x1+Data$loc.x2,regtype="ll",bwmethod="cv.aic")
regres=npreg(bw,gradients = TRUE)
grad=regres$grad
cgrad=sqrt(grad[,1]^2+grad[,2]^2)
Data$grad=cgrad
nseq=round(c(10,30,50,70,100)*1000^(2/9))
#46 139 232 325 464

##########    Universal Kriging Model    ##########
PEN_nonug=function(daTa,test,n)
{
  res=matrix(0,4,5)     
  findi=function(x){idex=which(daTa$loc.x1==x[1] & daTa$loc.x2==x[2])}
  
  #GBMED 
  MEDres=mined::SelectMinED(cbind(daTa$loc.x1,daTa$loc.x2),log(daTa$grad),n,1,2)
  MEDloc=MEDres$points
  MEDi=apply(MEDloc,1,findi)
  MEDy=daTa$y[MEDi]
  MEDGPres=likfit(as.geodata(cbind(MEDloc,MEDy)),trend ="1st",ini.cov.pars = c(1e6,10),fix.nugget = TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))  
  res[1,-5]=c(MEDGPres$beta,MEDGPres$phi/MEDGPres$sigmasq)
  MEDpred=krige.conv(as.geodata(cbind(MEDloc,MEDy)),locations = cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d ="1st",trend.l = "1st",beta = MEDGPres$beta,
                                         cov.pars = MEDGPres$cov.pars))
  MEDyhat=MEDpred$predict
  res[1,5]=sqrt(mean((test$y-MEDyhat)^2))
  
  # MaxPro
  iexist=sample(1:nrow(daTa),1)
  exist=daTa[iexist,1:2]
  MAPloc=MaxProAugment(exist,daTa[,1:2],nNew=n-1)$Design
  MAPi=apply(MAPloc,1,findi)
  MAPy=daTa$y[MAPi]
  MAPGPres=likfit(as.geodata(cbind(MAPloc,MAPy)),trend = "1st",ini.cov.pars = c(1e6,10),fix.nugget = TRUE)     #ini value
  res[2,-5]=c(MAPGPres$beta,MAPGPres$phi/MAPGPres$sigmasq)
  MAPpred=krige.conv(as.geodata(cbind(MAPloc,MAPy)),locations =cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = MAPGPres$beta,
                                         cov.pars = MAPGPres$cov.pars))
  MAPyhat=MAPpred$predict
  res[2,5]=sqrt(mean((test$y-MAPyhat)^2))
  
  #ABS
  ABSi=ABS_select(daTa$y,n)
  ABS=daTa[ABSi,1:3]
  ABSGPres=likfit(as.geodata(ABS),trend = "1st",ini.cov.pars = c(1e6,10),fix.nugget = TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))   
  res[3,-5]=c(ABSGPres$beta,ABSGPres$phi/ABSGPres$sigmasq)
  ABSpred=krige.conv(as.geodata(ABS),locations = cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = ABSGPres$beta,
                                         cov.pars = ABSGPres$cov.pars))
  ABSyhat=ABSpred$predict
  res[3,5]=sqrt(mean((test$y-ABSyhat)^2))
  
  #Uniform
  UNIi=sample(1:nrow(daTa),n)
  UNI=daTa[UNIi,1:3]
  UNIGPres=likfit(as.geodata(UNI),trend = "1st",ini.cov.pars = c(1e6,10),fix.nugget = TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))   
  res[4,-5]=c(UNIGPres$beta,UNIGPres$phi/UNIGPres$sigmasq)
  UNIpred=krige.conv(as.geodata(UNI),locations = cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = UNIGPres$beta,
                                         cov.pars = UNIGPres$cov.pars))
  UNIyhat=UNIpred$predict
  res[4,5]=sqrt(mean((test$y-UNIyhat)^2))
  
  return(res)
}


##########    Smoothing Splines Model    ##########
PEN_splines=function(daTa,test,n)
{
  res=rep(0,4)    
  findi=function(x){idex=which(daTa$loc.x1==x[1] & daTa$loc.x2==x[2])}
  
  #GBMED 
  MEDres=mined::SelectMinED(cbind(daTa$loc.x1,daTa$loc.x2),log(daTa$grad),n,1,2)
  MEDloc=MEDres$points
  MEDi=apply(MEDloc,1,findi)
  MEDfit=ssa(y~.+.^2, data=daTa[,-4], id.basis=MEDi,
              lambda = lambda.pro, theta=theta)
  MEDyhat=predict(MEDfit,newdata=test[,-3], se.fit=F)
  res[1]=sqrt(mean((test$y-MEDyhat)^2))
  
  # MaxPro
  iexist=sample(1:nrow(daTa),1)
  exist=daTa[iexist,1:2]
  MAPloc=MaxProAugment(exist,daTa[,1:2],nNew=n-1)$Design
  MAPi=apply(MAPloc,1,findi)
  MAPfit=ssa(y~.+.^2, data=daTa[,-4], id.basis=MAPi,
             lambda = lambda.pro, theta=theta)
  MAPyhat=predict(MAPfit,newdata=test[,-3], se.fit=F)
  res[2]=sqrt(mean((test$y-MAPyhat)^2))
  
  #ABS
  ABSi=ABS_select(daTa$y,n)
  ABSfit=ssa(y~.+.^2, data=daTa[,-4], id.basis=ABSi,
             lambda = lambda.pro, theta=theta)
  ABSyhat=predict(ABSfit,newdata=test[,-3], se.fit=F)
  res[3]=sqrt(mean((test$y-ABSyhat)^2))
  
  #Uniform
  UNIi=sample(1:nrow(daTa),n)
  UNIfit=ssa(y~.+.^2, data=daTa[,-4], id.basis=UNIi,
             lambda = lambda.pro, theta=theta)
  UNIyhat=predict(UNIfit,newdata=test[,-3], se.fit=F)
  res[4]=sqrt(mean((test$y-UNIyhat)^2))
  
  return(res)
}


PEN_main=function(nseq)
{
  it=sample(1:1720,1000)
  Dtrain=Data[it,]
  Dtest=Data[-it,-4]
  
  ### Universal kriging full data
  MLEres=likfit(as.geodata(Dtrain[,-4]),trend = "1st",ini.cov.pars = c(1e6,10),fix.nugget = TRUE)
  MLEpara=c(MLEres$beta,MLEres$cov.pars,MLEres$phi/MLEres$sigmasq)
  predres=krige.conv(as.geodata(Dtrain[,-4]),locations =Dtest[,-3],
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = MLEres$beta,
                                         cov.pars = MLEres$cov.pars))
  pred=sqrt(sum((predres$predict-Dtest[,3])^2)/nrow(Dtest))
  MlE=c(MLEpara,pred)
  
  ### Smoothing splines full data
  fit=ssa(y~.+.^2, data=Dtrain[,-4], id.basis=1:1000,
             lambda = lambda.pro, theta=theta)
  smoothpred=predict(fit,newdata=Dtest[,-3], se.fit=F)
  smooth=sqrt(mean((Dtest$y-smoothpred)^2))
  
  ### Universal kriging subsamples
  PEN1=function(n){return(PEN_nonug(Dtrain,Dtest,n))}
  krigout=sapply(nseq,PEN1,simplify = FALSE)
  
  ### Smoothing spline subsamples
  sam.size=ceiling(50*(1000^{1/4}))
  esps.re=esps(y~.+.^2, data=Dtrain[,-4], sam.size,r=3)
  r=3
  p=esps.re$p
  pro.gcv.lamb=esps.re$lambda*((1000/sam.size)^(-r/(p*r+1)))
  lambda.pro=log10(1000*pro.gcv.lamb)
  theta=esps.re$theta
  
  PEN2=function(n){return(PEN_splines(Dtrain,Dtest,n))}
  splineout=sapply(nseq,PEN2)
  
  res=list(Kriging=MlE,Splines=smooth,SubKriging=krigout,Subspline=splineout)
  return(res)
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultNS= foreach(i=1:200,
                 .combine=cbind,
                 .packages=c("geoR","mined","np","gss","fields","MaxPro")) %dopar% PEN_main(nseq)
stopCluster(cl)

MED_k_pa=matrix(0,5,200)
SBS_k_pa=matrix(0,5,200)
ABS_k_pa=matrix(0,5,200)
UNI_k_pa=matrix(0,5,200)
MED_s_pa=matrix(0,5,200)
SBS_s_pa=matrix(0,5,200)
ABS_s_pa=matrix(0,5,200)
UNI_s_pa=matrix(0,5,200)
Ful_k_pa=rep(0,200)
Ful_s_pa=rep(0,200)

for(i in 1:5){
  for(j in 1:200){
    MED_k_pa[i,j]=ResultNS[3,j][[1]][[i]][1,5]
    SBS_k_pa[i,j]=ResultNS[3,j][[1]][[i]][2,5]
    ABS_k_pa[i,j]=ResultNS[3,j][[1]][[i]][3,5]
    UNI_k_pa[i,j]=ResultNS[3,j][[1]][[i]][4,5]
    
    MED_s_pa[i,j]=ResultNS[4,j][[1]][1,i]
    SBS_s_pa[i,j]=ResultNS[4,j][[1]][2,i]
    ABS_s_pa[i,j]=ResultNS[4,j][[1]][3,i]
    UNI_s_pa[i,j]=ResultNS[4,j][[1]][4,i]
    
    Ful_k_pa[j]=ResultNS[1,j][[1]][7]
    Ful_s_pa[j]=ResultNS[2,j][[1]]
  }
}

data_NS=data.frame(kmse=c(apply(MED_k_pa,1,mean),apply(SBS_k_pa,1,mean),apply(ABS_k_pa,1,mean),apply(UNI_k_pa,1,mean),rep(mean(Ful_k_pa),5)),
                  ksd=c(apply(MED_k_pa,1,sd), apply(SBS_k_pa,1,sd),apply(ABS_k_pa,1,sd),apply(UNI_k_pa,1,sd),rep(sd(Ful_k_pa),5)),
                  smse=c(apply(MED_s_pa,1,mean),apply(SBS_s_pa,1,mean),apply(ABS_s_pa,1,mean),apply(UNI_s_pa,1,mean),rep(mean(Ful_s_pa),5)),
                  ssd=c(apply(MED_s_pa,1,sd), apply(SBS_s_pa,1,sd),apply(ABS_s_pa,1,sd),apply(UNI_s_pa,1,sd),rep(sd(Ful_s_pa),5)),
                  Method=factor(rep(c("GBMED","SBS","ABS","UNIF","Full"),each=5)),
                  k=rep(c(46,139,232,325,464),times=5))

pd=position_dodge(10)
#brewer.pal(8,'Set1')
p_K=ggplot(data_NS,aes(x=k,y=kmse,group=Method,colour=Method))+
  theme_light()+theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_rect(colour = "black"),
                      axis.text=element_text(size=12),
                      axis.title=element_text(size=14),
                      legend.position = "right",
                      legend.title = element_blank(),
                      legend.key.width=unit(2,"line"),
                      legend.key.height=unit(1,"line"))+
  geom_errorbar(aes(ymin=kmse-ksd, ymax=kmse+ksd), width=1, position=pd) +
  geom_line(aes(linetype=Method), position=pd)+
  geom_point(position=pd,size=1)+scale_shape_manual(values=c(1,1,1,1))+
  scale_linetype_manual(values=c(3,5,1,2,4))+
  scale_size_manual(values=c(1,1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A","#A65628","#E41A1C","#377EB8","#984EA3"))+
  xlab("k")+ylab("MSE")
p_K








