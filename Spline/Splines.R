#######       Smoothing Splines: Compare with SBS by Meng       #######
library("mined")
library("MaxPro")
library("np")
library("gss")

######   True Surface 
# function g: from the first setting of Meng's (SBS) simulation. 
g=function(t1,t2)
{
  g1=t1*t2
  g2=(2*t2-1)^2
  g3=sin(2*pi*t1)/(2-sin(2*pi*t1))
  g4=0.1*sin(2*pi*t2)+0.2*cos(2*pi*t2)+0.3*(sin(2*pi*t2))^2+0.4*(cos(2*pi*t2))^3+0.5*(sin(2*pi*t2))^3
  g5=sin(pi*(t1+t2))/(2-sin(pi*(t1+t2)))
  res=g1+g2+g3+g4+g5
  return(res)
}

xs=seq(0,1,length.out = 100)
xg=expand.grid(xs,xs)
zmat=matrix(g(xg[,1],xg[,2]),100,100)
image(xs,xs,zmat,xlab="",ylab="")

# function \eta: from the second setting of Ma's (ABS) simulation. 
eta=function(t1,t2)
{
  f1=2*sign(t1)*t1^2
  f2=3*sign(t2)*(abs(t2))^2
  Sigma=0.5*matrix(1,2,2)+0.5*diag(1,2)
  res=dmvnorm(c(f1,f2),c(0,0),Sigma)*36*abs(t1)*t2^2
  return(res)
}

xs=seq(-1,1,length.out = 100)
zmat=matrix(0,100,100)
for(i in 1:100)
{
  for(j in 1:100)
  {
    zmat[i,j]=eta(xs[i],xs[j])
  }
}
image(xs,xs,zmat,xlab="",ylab="")



##### Prediction with various methods
# interaction model; by ssanova() from package "gss"
Ssplineinter=function(Dtrain,Dtest,cgrad,n)
{
  res=rep(0,3)
  findi=function(x){idex=which(Dtrain$loc.x1==x[1] & Dtrain$loc.x2==x[2])}
  
  #Uniform
  iR=sample(1:nrow(Dtrain),n)
  Rssintr=ssanova(y~loc.x1*loc.x2,data = Dtrain,type = list(loc.x1="cubic",loc.x2="cubic"),id.basis = iR)
  res[1]=mean((predict(Rssintr,Dtest[,-3])-Dtest$y)^2)
  
  #MaxPro(SBS)
  iexist=sample(1:nrow(Dtrain),1)
  exist=Dtrain[iexist,-3]
  MaxProloc=MaxProAugment(exist,Dtrain[,-3],nNew=n-1)$Design
  iMAP=apply(MaxProloc,1,findi)
  MAPssintr=ssanova(y~loc.x1*loc.x2,data = Dtrain,type = list(loc.x1="cubic",loc.x2="cubic"),id.basis = iMAP)
  res[2]=mean((predict(MAPssintr,Dtest[,-3])-Dtest$y)^2)
  
  #MED
  MEDres=mined::SelectMinED(cbind(Dtrain$loc.x1,Dtrain$loc.x2),log(cgrad),n,1,2)
  MEDloc=MEDres$points
  iMED=apply(MEDloc,1,findi)
  MEDssintr=ssanova(y~loc.x1*loc.x2,data = Dtrain,type = list(loc.x1="cubic",loc.x2="cubic"),id.basis = iMED)
  res[3]=mean((predict(MEDssintr,Dtest[,-3])-Dtest$y)^2)
  
  return(res)
}

# additive model; by Meng's ssa() 
# need to load file ssa.R and esps-m.R
Sspline=function(Dtrain,Dtest,cgrad,n)
{
  res=rep(0,3)
  findi=function(x){idex=which(Dtrain$loc.x1==x[1] & Dtrain$loc.x2==x[2])}
  pre=esps(y~loc.x1+loc.x2,data = Dtrain,sam.size = 1e3,r=3)
  
  #Uniform
  iR=sample(1:nrow(Dtrain),n)
  Rssadd=ssa(y~loc.x1+loc.x2,data = Dtrain,theta = log10(pre$theta),lambda = log10(pre$lambda),
             type = list(loc.x1="cubic",loc.x2="cubic"),id.basis = iR)
  res[1]=mean((predict(Rssadd,Dtest[,-3])-Dtest$y)^2)
  
  #MaxPro
  iexist=sample(1:nrow(Dtrain),1)
  exist=Dtrain[iexist,-3]
  MaxProloc=MaxProAugment(exist,Dtrain[,-3],nNew=n-1)$Design
  iMAP=apply(MaxProloc,1,findi)
  MAPssadd=ssa(y~loc.x1+loc.x2,data = Dtrain,theta = log10(pre$theta),lambda = log10(pre$lambda),
               type = list(loc.x1="cubic",loc.x2="cubic"),id.basis = iMAP)
  res[2]=mean((predict(MAPssadd,Dtest[,-3])-Dtest$y)^2)
  
  #MED
  MEDres=mined::SelectMinED(cbind(Dtrain$loc.x1,Dtrain$loc.x2),log(cgrad),n,1,2)
  MEDloc=MEDres$points
  iMED=apply(MEDloc,1,findi)
  MEDssadd=ssa(y~loc.x1+loc.x2,data = Dtrain,theta = log10(pre$theta),lambda = log10(pre$lambda),
               type = list(loc.x1="cubic",loc.x2="cubic"),id.basis = iMED)
  res[3]=mean((predict(MEDssadd,Dtest[,-3])-Dtest$y)^2)
  
  return(res)
}



##### Case 1: data are generated from function g of size n={2^10, 2^11, 2^12, 2^13, 2^14}; 
#             subsample size k1=5n^(2/9) and k2=10n^(1/9) (as Meng's); 
#             noise \sigma^2=0.4^2, SNR=6.5;
Sspline_Main1=function(N)
{
  ### generate the data : g function
  Data=data.frame(loc.x1=runif(N+5000),loc.x2=runif(N+5000))
  for(i in 1:nrow(Data))
  {
    Data$y[i]=g(Data[i,1],Data[i,2])
  }
  Data$y=Data$y+0.4*rnorm(nrow(Data))   
  itest=sample(1:nrow(Data),5e3)
  Dtrain=Data[-itest,]
  Dtest=Data[itest,]
  ib=sample(1:nrow(Dtrain),1e3)
  Dataib=Dtrain[ib,]
  bw=npregbw(Dataib$y~Dataib$loc.x1+Dataib$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  regres=npreg(Dtrain$y~Dtrain$loc.x1+Dtrain$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  
  res=matrix(0,2,3)
  res[1,]=Ssplineinter(Dtrain,Dtest,cgrad,round(5*(N)^(2/9)))
  res[2,]=Ssplineinter(Dtrain,Dtest,cgrad,round(10*(N)^(1/9)))
  return(res)
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultS1 = foreach(i=1:100,
                   .combine=cbind,
                   .packages=c("mined","np","MaxPro","gss")) %dopar% sapply(c(2^10,2^11,2^12,2^13,2^14),Sspline_Main1,simplify = FALSE)
stopCluster(cl)


##### Case 3: data are generated from function \eta of size n={2^10, 2^11, 2^12, 2^13, 2^14}; 
#             subsample size k1=5n^(2/9) and k2=10n^(1/9) (as Meng's); 
#             noise \sigma^2=0.1^2, SNR=3.6;
Sspline_Main3=function(N)
{
  ### generate the data : g function
  Data=data.frame(loc.x1=runif(N+5000),loc.x2=runif(N+5000))
  for(i in 1:nrow(Data))
  {
    Data$y[i]=eta(Data[i,1],Data[i,2])
  }
  Data$y=Data$y+0.1*rnorm(nrow(Data))   
  itest=sample(1:nrow(Data),5e3)
  Dtrain=Data[-itest,]
  Dtest=Data[itest,]
  ib=sample(1:nrow(Dtrain),1e3)
  Dataib=Dtrain[ib,]
  bw=npregbw(Dataib$y~Dataib$loc.x1+Dataib$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  regres=npreg(Dtrain$y~Dtrain$loc.x1+Dtrain$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  
  res=matrix(0,2,3)
  res[1,]=Ssplineinter(Dtrain,Dtest,cgrad,round(5*(N)^(2/9)))
  res[2,]=Ssplineinter(Dtrain,Dtest,cgrad,round(10*(N)^(1/9)))
  return(res)
}


cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultS3 = foreach(i=1:100,
                   .combine=cbind,
                   .packages=c("mined","np","MaxPro","gss","mvtnorm")) %dopar% sapply(c(2^10,2^11,2^12,2^13,2^14),Sspline_Main3,simplify = FALSE)
stopCluster(cl)



##### Case 2: data are generated from function \eta of size n=5000; 
#             subsample size k={20, 30, 40, 50}; 
#             noise \sigma^2=0.1^2, SNR=3.6;
Sspline_Main2=function()
{
  ### generate the data
  Data=data.frame(loc.x1=runif(1e4,-1,1),loc.x2=runif(1e4,-1,1))
  for(i in 1:nrow(Data))
  {
    Data$y[i]=eta(Data[i,1],Data[i,2])    # eta function from Ma
  }
  Data$y=Data$y+rnorm(1e4,mean=0,sd=0.1)   
  itrain=sample(1:1e4,5e3)
  Dtrain=Data[itrain,]
  Dtest=Data[-itrain,]
  ib=sample(1:5e3,1e3)
  Dataib=Dtrain[ib,]
  bw=npregbw(Dataib$y~Dataib$loc.x1+Dataib$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  regres=npreg(Dtrain$y~Dtrain$loc.x1+Dtrain$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  
  nseq=c(20,30,40,50)
  Ssplinen=function(n){return(Ssplineinter(Dtrain,Dtest,cgrad,n))}
  out=sapply(nseq,Ssplinen,simplify = FALSE)
  return(out)
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultS2 = foreach(i=1:100,
                   .combine=cbind,
                   .packages=c("mined","np","MaxPro","gss","mvtnorm")) %dopar% Sspline_Main2()
stopCluster(cl)



##### Case 4: data are generated from function g of size n=5000; 
#             subsample size k={20, 30, 40, 50}; 
#             noise \sigma^2=0.4^2, SNR=6.5;
Sspline_Main4=function()
{
  ### generate the data
  Data=data.frame(loc.x1=runif(1e4,-1,1),loc.x2=runif(1e4,-1,1))
  for(i in 1:nrow(Data))
  {
    Data$y[i]=g(Data[i,1],Data[i,2])    # g function from Ma
  }
  Data$y=Data$y+rnorm(1e4,mean=0,sd=0.4)   #SNR=5
  itrain=sample(1:1e4,5e3)
  Dtrain=Data[itrain,]
  Dtest=Data[-itrain,]
  ib=sample(1:5e3,1e3)
  Dataib=Dtrain[ib,]
  bw=npregbw(Dataib$y~Dataib$loc.x1+Dataib$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  regres=npreg(Dtrain$y~Dtrain$loc.x1+Dtrain$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  
  nseq=c(20,30,40,50)
  Ssplinen=function(n){return(Ssplineinter(Dtrain,Dtest,cgrad,n))}
  out=sapply(nseq,Ssplinen,simplify = FALSE)
  return(out)
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultS4 = foreach(i=1:100,
                   .combine=cbind,
                   .packages=c("mined","np","MaxPro","gss","mvtnorm")) %dopar% Sspline_Main4()
stopCluster(cl)


##### Case 5: data are generated from function g of size n=5000; 
#             subsample size k={20, 30, 40, 50}; 
#             noise \sigma^2=0.7^2, SNR=2;
Sspline_Main5=function()
{
  ### generate the data
  Data=data.frame(loc.x1=runif(1e4,-1,1),loc.x2=runif(1e4,-1,1))
  for(i in 1:nrow(Data))
  {
    Data$y[i]=g(Data[i,1],Data[i,2])    # g function from Ma
  }
  Data$y=Data$y+rnorm(1e4,mean=0,sd=0.7)   #SNR=2
  itrain=sample(1:1e4,5e3)
  Dtrain=Data[itrain,]
  Dtest=Data[-itrain,]
  ib=sample(1:5e3,1e3)
  Dataib=Dtrain[ib,]
  bw=npregbw(Dataib$y~Dataib$loc.x1+Dataib$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  regres=npreg(Dtrain$y~Dtrain$loc.x1+Dtrain$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  
  nseq=c(20,30,40,50)
  Ssplinen=function(n){return(Ssplineinter(Dtrain,Dtest,cgrad,n))}
  out=sapply(nseq,Ssplinen,simplify = FALSE)
  return(out)
}

cl<- makeCluster(8) 
registerDoParallel(cl) 
ResultS5 = foreach(i=1:100,
                   .combine=cbind,
                   .packages=c("mined","np","MaxPro","gss","mvtnorm")) %dopar% Sspline_Main5()

stopCluster(cl)

