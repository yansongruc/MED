#######       Linear Model: Compare with IBOSS       #######
library("IBOSS")
library("np")
library("mined")
library("doParallel")
library("foreach")


#### Dtrain: training data; Dtest: testing data; cgrad: the L_2-norm of gradients of Dtrain; n: subsample size
#### OLS
Linear=function(Dtrain,Dtest,cgrad,n)
{
  res=matrix(0,3,4)
  findi=function(x){idex=which(Dtrain[,1]==x[1] & Dtrain[,2]==x[2])}
  
  #MED
  MEDres=mined::SelectMinED(as.matrix(Dtrain[,-3]),log(cgrad),n,1,2)
  MEDloc=MEDres$points
  MEDi=apply(MEDloc,1,findi)
  MEDsub=Dtrain[MEDi,]
  MEDfit=lm(y~loc.x1+loc.x2,MEDsub)
  res[1,-4]=MEDfit$coefficients
  MEDpred=predict(MEDfit,newdata=Dtest[,-3])
  res[1,4]=mean((Dtest[,3]-MEDpred)^2)
  
  #IBOSS
  IBOSSfit=iboss.od(cbind(Dtrain$loc.x1,Dtrain$loc.x2), Dtrain$y, n)
  res[2,-4]=IBOSSfit$beta
  IBOSSpred=cbind(1,as.matrix(Dtest[,-3]))%*%res[2,-4]
  res[2,4]=mean((Dtest[,3]-IBOSSpred)^2)
  
  #Uniform
  Ri=sample(1:nrow(Dtrain),n)
  Rsub=Dtrain[Ri,]
  Rfit=lm(y~loc.x1+loc.x2,Rsub)
  res[3,-4]=Rfit$coefficients
  Rpred=predict(Rfit,newdata=Dtest[,-3])
  res[3,4]=mean((Dtest[,3]-Rpred)^2)
  
  return(res)
}

#### GLS with known heteroscedasticity sigma_i^2=x_i1^2+x_i2^2
GLinear=function(Dtrain,Dtest,cgrad,n)
{
  res=matrix(0,3,4)
  findi=function(x){idex=which(Dtrain[,1]==x[1] & Dtrain[,2]==x[2])}
  
  #MED
  MEDres=mined::SelectMinED(as.matrix(Dtrain[,-3]),log(cgrad),n,1,2)
  MEDloc=MEDres$points
  MEDi=apply(MEDloc,1,findi)
  MEDsub=Dtrain[MEDi,]
  MEDweight=1/apply(MEDsub[,-3]^2,1,sum)
  MEDfit=lm(y~loc.x1+loc.x2,MEDsub,weights = MEDweight)
  res[1,-4]=MEDfit$coefficients
  MEDpred=predict(MEDfit,newdata=Dtest[,-3])
  res[1,4]=mean((Dtest[,3]-MEDpred)^2)
  
  #IBOSS
  IBOSSfit=iboss.od(cbind(Dtrain$loc.x1,Dtrain$loc.x2), Dtrain$y, n)
  IBOSSi=IBOSSfit$index
  IBOSSsub=Dtrain[IBOSSi,]
  IBOSSweight=1/apply(IBOSSsub[,-3]^2,1,sum)
  IBOSSfit1=lm(y~loc.x1+loc.x2,IBOSSsub,weights = IBOSSweight)
  res[2,-4]=IBOSSfit1$coefficients
  IBOSSpred=cbind(1,as.matrix(Dtest[,-3]))%*%res[2,-4]
  res[2,4]=mean((Dtest[,3]-IBOSSpred)^2)
  
  #Uniform
  Ri=sample(1:nrow(Dtrain),n)
  Rsub=Dtrain[Ri,]
  Rweight=1/apply(Rsub[,-3]^2,1,sum)
  Rfit=lm(y~loc.x1+loc.x2,Rsub,weights = Rweight)
  res[3,-4]=Rfit$coefficients
  Rpred=predict(Rfit,newdata=Dtest[,-3])
  res[3,4]=mean((Dtest[,3]-Rpred)^2)
  
  return(res)
}


####### Case 1: X~U[0,1]*U[0,1]; Fix N=5000; subsample size=c(100,200,400,600,800,1000)
Linear_Main=function()
{
  #generate the data, N=1e4, Dtrain=5e3, Dtest=5e3; beta=(1,1,1), sigma^2=9; random X; 
  beta=rep(1,3)
  X=matrix(runif(2e4),1e4,2)               # generate covariate X
  Y=cbind(1,X)%*%beta+rnorm(1e4, 0, 3)     # generate Y
  Data=data.frame(loc.x1=X[,1],loc.x2=X[,2],y=Y)
  itrain=sample(1:1e4,5e3)
  Dtrain=Data[itrain,]
  Dtest=Data[-itrain,]
  ibw=sample(1:5e3,1e3)
  Dataibw=Dtrain[ibw,]
  bw=npregbw(Dataibw$y~Dataibw$loc.x1+Dataibw$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  regres=npreg(Dtrain$y~Dtrain$loc.x1+Dtrain$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  
  Fullfit=lm(y~loc.x1+loc.x2,Dtrain)
  Fullpred=predict(Fullfit,newdata=Dtest[,-3])
  Fullout=c(Fullfit$coefficients,mean((Dtest[,3]-Fullpred)^2))
  
  nseq=c(100,200,400,600,800,1000)
  Linearn=function(n){return(Linear(Dtrain,Dtest,cgrad,n))}
  out=sapply(nseq,Linearn,simplify = FALSE)
  
  return(list(Full=Fullout,Sub=out))
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultLi1 = foreach(i=1:500,
                    .combine=cbind,
                    .packages=c("mined","np","IBOSS")) %dopar% Linear_Main()
stopCluster(cl)


####### Case 2: Nonuniform X; Fix N=5000; subsample size=c(100,200,400,600,800,1000)
Linear_Main2=function()
{
  #generate the data
  beta=rep(1,3)
  X1=matrix(runif(18000,0,0.5),9000,2)
  X2=cbind(runif(500,0,0.5),runif(500,0.5,1))
  X3=cbind(runif(250,0.5,1),runif(250,0.5,1))
  X4=cbind(runif(250,0.5,1),runif(250,0,0.5))
  X=rbind(X1,X2,X3,X4)
  Y=cbind(1,X)%*%beta+rnorm(1e4, 0, 3)
  Data=data.frame(loc.x1=X[,1],loc.x2=X[,2],y=Y)
  itrain=sample(1:1e4,5e3)
  Dtrain=Data[itrain,]
  Dtest=Data[-itrain,]
  ibw=sample(1:5e3,1e3)
  Dataibw=Dtrain[ibw,]
  bw=npregbw(Dataibw$y~Dataibw$loc.x1+Dataibw$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  regres=npreg(Dtrain$y~Dtrain$loc.x1+Dtrain$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  
  Fullfit=lm(y~loc.x1+loc.x2,Dtrain)
  Fullpred=predict(Fullfit,newdata=Dtest[,-3])
  Fullout=c(Fullfit$coefficients,mean((Dtest[,3]-Fullpred)^2))
  
  nseq=c(100,200,400,600,800,1000)
  Linearn=function(n){return(Linear(Dtrain,Dtest,cgrad,n))}
  out=sapply(nseq,Linearn,simplify = FALSE)
  
  return(list(Full=Fullout,Sub=out))
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultLi2 = foreach(i=1:500,
                    .combine=cbind,
                    .packages=c("mined","np","IBOSS")) %dopar% Linear_Main2()
stopCluster(cl)


######## Case 3: Nonuniform X; heteroscedasticity; Fix N; n(or k)=c(100,200,400,600,800,1000)
Linear_Main4=function()
{
  #generate the data
  beta=rep(1,3)
  X1=matrix(runif(18000,0,0.5),9000,2)
  X2=cbind(runif(500,0,0.5),runif(500,0.5,1))
  X3=cbind(runif(250,0.5,1),runif(250,0.5,1))
  X4=cbind(runif(250,0.5,1),runif(250,0,0.5))
  X=rbind(X1,X2,X3,X4)
  Ve=sqrt(apply(X^2,1,sum))
  e=rep(0,1e4)
  for(i in 1:1e4)
  {
    e[i]=rnorm(1,0,Ve[i])
  }
  Y=cbind(1,X)%*%beta+e
  Data=data.frame(loc.x1=X[,1],loc.x2=X[,2],y=Y)
  itrain=sample(1:1e4,5e3)
  Dtrain=Data[itrain,]
  Dtest=Data[-itrain,]
  ibw=sample(1:5e3,1e3)
  Dataibw=Dtrain[ibw,]
  bw=npregbw(Dataibw$y~Dataibw$loc.x1+Dataibw$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  regres=npreg(Dtrain$y~Dtrain$loc.x1+Dtrain$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  
  Fullfit=lm(y~loc.x1+loc.x2,Dtrain)
  Fullpred=predict(Fullfit,newdata=Dtest[,-3])
  Fullout=c(Fullfit$coefficients,mean((Dtest[,3]-Fullpred)^2))
  
  nseq=c(100,200,400,600,800,1000)
  GLinearn=function(n){return(GLinear(Dtrain,Dtest,cgrad,n))}
  out=sapply(nseq,GLinearn,simplify = FALSE)
  
  return(list(Full=Fullout,Sub=out))
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultLi4 = foreach(i=1:500,
                    .combine=cbind,
                    .packages=c("mined","np","IBOSS")) %dopar% Linear_Main4()
stopCluster(cl)


######## Case 4: Normal X; heteroscedasticity; Fix N; n(or k)=c(100,200,400,600,800,1000)
Linear_Main5=function()
{
  #generate the data
  beta=rep(1,3)
  mu=rep(0,2)
  Sigma=0.5*matrix(1,2,2)+0.5*diag(1,2)
  X=mvrnorm(1e4,mu,Sigma)
  Ve=sqrt(apply(X^2,1,sum))
  e=rep(0,1e4)
  for(i in 1:1e4)
  {
    e[i]=rnorm(1,0,Ve[i])
  }
  Y=cbind(1,X)%*%beta+e
  Data=data.frame(loc.x1=X[,1],loc.x2=X[,2],y=Y)
  itrain=sample(1:1e4,5e3)
  Dtrain=Data[itrain,]
  Dtest=Data[-itrain,]
  ibw=sample(1:5e3,1e3)
  Dataibw=Dtrain[ibw,]
  bw=npregbw(Dataibw$y~Dataibw$loc.x1+Dataibw$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  regres=npreg(Dtrain$y~Dtrain$loc.x1+Dtrain$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  
  Fullweight=1/apply(Dtrain[,-3]^2,1,sum)
  Fullfit=lm(y~loc.x1+loc.x2,Dtrain,weights = Fullweight)
  Fullpred=predict(Fullfit,newdata=Dtest[,-3])
  Fullout=c(Fullfit$coefficients,mean((Dtest[,3]-Fullpred)^2))
  
  nseq=c(100,200,400,600,800,1000)
  GLinearn=function(n){return(GLinear(Dtrain,Dtest,cgrad,n))}
  out=sapply(nseq,GLinearn,simplify = FALSE)
  
  return(list(Full=Fullout,Sub=out))
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultLi5 = foreach(i=1:500,
                    .combine=cbind,
                    .packages=c("mined","np","IBOSS","MASS")) %dopar% Linear_Main5()

stopCluster(cl)
