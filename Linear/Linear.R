#######       Linear Model: Compare with IBOSS       #######
library("IBOSS")
library("np")
library("mined")

####### Case 1: X~U[0,1]*U[0,1]; Fix N; n(or k)=c(100,200,400,600,800,1000)
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

Linear_Main=function()
{
  #generate the data, N=1e4, Dtrain=5e3, Dtest=5e3; beta=(1,1,1), sigma^2=9; random X; 
  beta=rep(1,3)
  X=matrix(runif(2e4),1e4,2)
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
ResultLi1 = foreach(i=1:500,
                    .combine=cbind,
                    .packages=c("mined","np","IBOSS")) %dopar% Linear_Main()
stopCluster(cl)


####### Case 2: Nonuniform X; Fix N; n(or k)=c(100,200,400,600,800,1000)
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

######## Nonuniform X; heteroscedasticity; Fix N; n(or k)=c(100,200,400,600,800,1000)
Linear_Main3=function()
{
  #generate the data
  beta=rep(1,3)
  X1=matrix(runif(18000,0,0.5),9000,2)
  X2=cbind(runif(500,0,0.5),runif(500,0.5,1))
  X3=cbind(runif(250,0.5,1),runif(250,0.5,1))
  X4=cbind(runif(250,0.5,1),runif(250,0,0.5))
  X=rbind(X1,X2,X3,X4)
  Y1=cbind(1,X1)%*%beta+rnorm(9e3, 0, 1)
  Y2=cbind(1,X2)%*%beta+rnorm(500, 0, 2)
  Y3=cbind(1,X3)%*%beta+rnorm(250, 0, 3)
  Y4=cbind(1,X4)%*%beta+rnorm(250, 0, 2)
  Y=c(Y1,Y2,Y3,Y4)
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
ResultLi3 = foreach(i=1:500,
                    .combine=cbind,
                    .packages=c("mined","np","IBOSS")) %dopar% Linear_Main3()
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
  e=apply(X^2,1,sum)
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
  Linearn=function(n){return(Linear(Dtrain,Dtest,cgrad,n))}
  out=sapply(nseq,Linearn,simplify = FALSE)
  
  return(list(Full=Fullout,Sub=out))
}
cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultLi4 = foreach(i=1:500,
                    .combine=cbind,
                    .packages=c("mined","np","IBOSS")) %dopar% Linear_Main4()
stopCluster(cl)
