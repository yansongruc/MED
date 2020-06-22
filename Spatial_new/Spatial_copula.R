##########               Spatial Simulation 2: copula               ##########
library("geoR")
library("MaxPro")
library("mined")
library("fields")
library("np")
library("doParallel")
library("foreach")
library("mvtnorm")

##########    Parameter Estimation Without Nugget Effect   ##########
gene_y = function(x){
  t1=x[1];t2=x[2]
  f1=sign(t1)*t1^2
  f2=sign(t2)*(abs(t2))^3
  Sigma=0.5*matrix(1,2,2)+0.5*diag(1,2)
  res=dmvnorm(c(f1,f2),c(0,0),Sigma)*6*abs(t1)*t2*t2
  return(res)
}


PE_nonug=function(daTa,test,n,seed)
{
  res=matrix(0,4,6)     
  findi=function(x){idex=which(daTa$loc.x1==x[1] & daTa$loc.x2==x[2])}
  
  #GBMED 
  res[1,6]=proc.time()[3]
  MEDres=mined::SelectMinED(cbind(daTa$loc.x1,daTa$loc.x2),log(daTa$grad),n,1,2)
  MEDloc=MEDres$points
  res[1,6]=proc.time()[3]-res[1,6]
  MEDi=apply(MEDloc,1,findi)
  MEDy=daTa$y[MEDi]
  MEDGPres=likfit(as.geodata(cbind(MEDloc,MEDy)),trend ="1st",ini.cov.pars = c(1,1),fix.nugget = TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))  
  res[1,1:4]=c(MEDGPres$beta,MEDGPres$phi/MEDGPres$sigmasq)
  MEDpred=krige.conv(as.geodata(cbind(MEDloc,MEDy)),locations = cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d ="1st",trend.l = "1st",beta = MEDGPres$beta,
                                         cov.pars = MEDGPres$cov.pars))
  MEDyhat=MEDpred$predict
  res[1,5]=log(mean((test$y-MEDyhat)^2))
  
  #SBS
  set.seed(seed)
  design=MaxProLHD(n,2)$Design
  design=4*design-2
  res[2,6]=proc.time()[3]
  SBSi=nabor::knn(daTa[,1:2],design,k = 1)$nn.idx
  res[2,6]=proc.time()[3]-res[2,6]
  SBS=daTa[SBSi,1:3]
  SBSGPres=likfit(as.geodata(SBS),trend = "1st",ini.cov.pars = c(1,1),fix.nugget = TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf))) 
  res[2,1:4]=c(SBSGPres$beta,SBSGPres$phi/SBSGPres$sigmasq)
  SBSpred=krige.conv(as.geodata(SBS),locations = cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = SBSGPres$beta,
                                         cov.pars = SBSGPres$cov.pars))
  SBSyhat=SBSpred$predict
  res[2,5]=log(mean((test$y-SBSyhat)^2))
  
  #MaxPro
  #set.seed(seed)
  #iexist=sample(1:nrow(daTa),1)
  #exist=daTa[iexist,1:2]
  #MAPloc=MaxProAugment(exist,daTa[,1:2],nNew=n-1)$Design
  #MAPi=apply(MAPloc,1,findi)
  #MAPy=daTa$y[MAPi]
  #MAPGPres=likfit(as.geodata(cbind(MAPloc,MAPy)),trend = "1st",ini.cov.pars = c(1,1),fix.nugget = TRUE,
  #                limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf))) 
  #res[2,-5]=c(MAPGPres$beta,MAPGPres$phi/MAPGPres$sigmasq)
  #MAPpred=krige.conv(as.geodata(cbind(MAPloc,MAPy)),locations = cbind(test$loc.x1,test$loc.x2),
  #                  krige=krige.control(trend.d = "1st",trend.l = "1st",beta = MAPGPres$beta,
  #                                       cov.pars = MAPGPres$cov.pars))
  #MAPyhat=MAPpred$predict
  #res[2,5]=log(mean((test$y-MAPyhat)^2))
  
  #ABS
  set.seed(seed)
  res[3,6]=proc.time()[3]
  ABSi=ABS_select(daTa$y,n)
  res[3,6]=proc.time()[3]-res[3,6]
  ABS=daTa[ABSi,1:3]
  ABSGPres=likfit(as.geodata(ABS),trend = "1st",ini.cov.pars = c(1,1),fix.nugget=TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))   
  res[3,1:4]=c(ABSGPres$phi,ABSGPres$phi/ABSGPres$sigmasq)
  ABSpred=krige.conv(as.geodata(ABS),locations = cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = ABSGPres$beta,
                                         cov.pars = ABSGPres$cov.pars))
  ABSyhat=ABSpred$predict
  res[3,5]=log(mean((test$y-ABSyhat)^2))
  
  #Uniform
  set.seed(seed)
  res[4,6]=proc.time()[3]
  UNIi=sample(1:nrow(daTa),n)
  res[4,6]=proc.time()[3]-res[4,6]
  UNI=daTa[UNIi,1:3]
  UNIGPres=likfit(as.geodata(UNI),trend = "1st",ini.cov.pars = c(1,1),fix.nugget = TRUE,
                  limits = pars.limits(phi = c(1e-3,+Inf),sigmasq = c(1e-3,+Inf)))   
  res[4,1:4]=c(UNIGPres$beta,UNIGPres$phi/UNIGPres$sigmasq)
  UNIpred=krige.conv(as.geodata(UNI),locations = cbind(test$loc.x1,test$loc.x2),
                     krige=krige.control(trend.d = "1st",trend.l = "1st",beta = UNIGPres$beta,
                                         cov.pars = UNIGPres$cov.pars))
  UNIyhat=UNIpred$predict
  res[4,5]=log(mean((test$y-UNIyhat)^2))
  
  return(res)
}


PEnonug_main=function(seed,snr)
{
  m=1e4
  ntest=5e3
  nseq=c(30,50,100,150,200)
  
  ### generate the data
  set.seed(seed)
  loc=matrix(runif(2*(m+ntest)),(m+ntest),2)
  loc=4*loc-2
  y1=apply(loc,1,gene_y)
  set.seed(seed)
  y2=grf((m+ntest),grid=loc,cov.pars = c(1,0.25))$data
  Data=data.frame(loc.x1=loc[,1],loc.x2=loc[,2],y=y1+sqrt(var(y1)/var(y2)/snr)*y2)
  set.seed(seed)
  itest=sample(1:(m+ntest),ntest)
  Dtrain=Data[-itest,]
  Dtest=Data[itest,]
  
  ### gradient evaluation
  set.seed(seed)
  ib=sample(1:nrow(Dtrain),1e3)
  Dataib=Dtrain[ib,]
  bw=npregbw(Dataib$y~Dataib$loc.x1+Dataib$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  regres=npreg(Dtrain$y~Dtrain$loc.x1+Dtrain$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  Dtrain$grad=cgrad
  
  ### parameter estimation
  PE1=function(n){return(PE_nonug(Dtrain,Dtest,n,seed))}
  out=sapply(nseq,PE1,simplify = FALSE)
  return(out)
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultC_5= foreach(i=1:200,
                 .combine=cbind,
                 .packages=c("geoR","mined","np","MaxPro","nabor","mvtnorm")) %dopar% PEnonug_main(666*i+888,5)
stopCluster(cl)

MED_c5_pa=matrix(0,5,200)
MAP_c5_pa=matrix(0,5,200)
ABS_c5_pa=matrix(0,5,200)
UNI_c5_pa=matrix(0,5,200)
for(i in 1:5){
  for(j in 1:200){
    MED_c5_pa[i,j]=ResultC_5[i,j][[1]][1,5]
    MAP_c5_pa[i,j]=ResultC_5[i,j][[1]][2,5]
    ABS_c5_pa[i,j]=ResultC_5[i,j][[1]][3,5]
    UNI_c5_pa[i,j]=ResultC_5[i,j][[1]][4,5]
  }
}

data_C5=data.frame(mse=c(apply(MED_c5_pa,1,mean),apply(MAP_c5_pa,1,mean),apply(ABS_c5_pa,1,mean),apply(UNI_c5_pa,1,mean)),
                    sd=c(apply(MED_c5_pa,1,sd), apply(MAP_c5_pa,1,sd),apply(ABS_c5_pa,1,sd),apply(UNI_c5_pa,1,sd)),
                    Method=factor(rep(c("GBMED","SBS","ABS","UNIF"),each=5)),
                    k=rep(c(30,50,100,150,200),times=4))

pd=position_dodge(5)
p_C5=ggplot(data_C5,aes(x=k,y=mse,group=Method,colour=Method))+
  theme_light()+theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_rect(colour = "black"),
                      axis.text=element_text(size=12),
                      axis.title=element_text(size=14),
                      legend.position = "right",
                      legend.title = element_blank(),
                      legend.key.width=unit(2,"line"),
                      legend.key.height=unit(1,"line"))+
  geom_errorbar(aes(ymin=mse-sd, ymax=mse+sd), width=1, position=pd) +
  geom_line(aes(linetype=Method), position=pd)+
  geom_point(position=pd, aes(shape=Method),size=2)+
  scale_linetype_manual(values=c(3,1,2,4))+
  scale_size_manual(values=c(1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A","#E41A1C", "#377EB8","#984EA3"))+
  xlab("k")+ylab("log(MSE)")
p_C5

MED_c2_tm=matrix(0,5,200)
MAP_c2_tm=matrix(0,5,200)
ABS_c2_tm=matrix(0,5,200)
UNI_c2_tm=matrix(0,5,200)
for(i in 1:5){
  for(j in 1:200){
    MED_c2_tm[i,j]=ResultC_2[i,j][[1]][1,6]
    MAP_c2_tm[i,j]=ResultC_2[i,j][[1]][2,6]
    ABS_c2_tm[i,j]=ResultC_2[i,j][[1]][3,6]
    UNI_c2_tm[i,j]=ResultC_2[i,j][[1]][4,6]
  }
}

apply(MED_c2_tm,1,mean)
apply(MED_c2_tm,1,sd)

