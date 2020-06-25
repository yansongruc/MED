library(mvtnorm)
library(MASS)
library(np)
library(mined)
library(MaxPro)
library(nabor)
library(randtoolbox)
library(gss)
library(foreach)
library(doParallel)
library(ggplot2)
#source("ssa.R")
#source("esps-m.R")
#source("ABS_select.R")
#source("MED_select.R")

################################################################################
# Setting 2: function copula in Ma's second setting or Meng's supplement       #
#            Uniform locations                                                 #
#            Full sample sizes are c(2^12, 2^13, 2^14, 2^15, 2^16)             #
#            Subsample sizes are k1=10*N^(1/9) and k2=5*N^(2/9)                #
################################################################################
set.seed(100)
nloop=100
NN=2^(seq(12,16,1))
ntest=5000
SNR=2          # SNR=5 or 2
sub1_meta=round(10*NN^(1/9))
#sub1_meta=c(25,27,29,32,34)
sub2_meta=round(5*NN^(2/9))
#sub2_meta=c(32,37,43,50,59)

# function copula
pp=2
gene_y = function(x){
  t1=x[1];t2=x[2]
  f1=sign(t1)*t1^2
  f2=sign(t2)*(abs(t2))^3
  Sigma=0.5*matrix(1,2,2)+0.5*diag(1,2)
  res=dmvnorm(c(f1,f2),c(0,0),Sigma)*6*abs(t1)*t2*t2
  return(res)
}

# Generate testing data: uniform; on [-2,2]^2
set.seed(100)
X_test=as.data.frame(matrix(runif(2*ntest,-2,2),ntest,2)) 
y.t=apply(X_test,1,gene_y)

# using a part of data to find suitable p, lambda, theta by esps-m.R
sam.size=ceiling(50*(NN[5]^{1/4}))
set.seed(100)
X=as.data.frame(matrix(runif(2*NN[5],-2,2),NN[5],2))
y=apply(X,1,gene_y)
esps.re=esps(y~.+.^2, data=data.frame(y, X), sam.size,r=3)

# do parallel
cl <- makeCluster(8)
registerDoParallel(cl)

UNIF1_meta=UNIF2_meta=ABS1_meta=ABS2_meta=SBS1_meta=SBS2_meta=MED1_meta=MED2_meta=matrix(0,nloop,length(NN))
for(jj in 1:length(NN)){
  
  N=NN[jj]
  sub1=sub1_meta[jj]
  sub2=sub2_meta[jj]
  
  #Design matrix
  set.seed(100)
  X=as.data.frame(matrix(runif(2*N,-2,2),N,2)) 
  y.true=apply(X,1,gene_y)
  
  ###################################################################
  setseed=100+821*jj
  set.seed(setseed)
  
  # generate the data?
  v=var(y.true)
  y=y.true + rnorm(N, 0, sqrt(v/SNR))   
  
  r <- 3
  p <- esps.re$p
  pro.gcv.lamb <- esps.re$lambda*((N/sam.size)^(-r/(p*r+1)))
  lambda.pro <- log10(N*pro.gcv.lamb)
  theta <- esps.re$theta
  
  ####################################################################3
  aaa<- foreach(i = 1:nloop, .packages=c("gss","nabor","randtoolbox","np","mined","MaxPro"), .combine="cbind") %dopar% {
    
    setseed=100+821*jj+123*i
    set.seed(setseed)
    
    v=var(y.true)
    y=y.true + rnorm(N, 0, sqrt(v/SNR))   #generate the data
    
    ###############################
    # Unif
    unif.ind = sample.int(N,size=sub1,replace=F)
    fit1 <- ssa(y~.+.^2, data=X, id.basis=unif.ind,
                lambda = lambda.pro, theta=theta)
    unif_1.temp <- predict(fit1, newdata=X_test, se.fit=F)
    
    unif.ind = sample.int(N,size=sub2,replace=F)
    fit2 <- ssa(y~.+.^2, data=X, id.basis=unif.ind,
                lambda = lambda.pro, theta=theta)
    unif_2.temp <- predict(fit2, newdata=X_test, se.fit=F)
    
    ###############################
    # ABS
    abs.ind = ABS_select(y,sub1)
    fit1 <- ssa(y~.+.^2, data=X, id.basis=abs.ind,
                lambda = lambda.pro, theta=theta)
    abs_1.temp <- predict(fit1, newdata=X_test, se.fit=F)
    
    abs.ind = ABS_select(y,sub2)
    fit2 <- ssa(y~.+.^2, data=X, id.basis=abs.ind,
                lambda = lambda.pro, theta=theta)
    abs_2.temp <- predict(fit2, newdata=X_test, se.fit=F)
    
    ###############################3
    # SBS
    #design1<- sobol(n=sub1 ,dim = pp, init = T, scrambling = 1, seed=setseed)
    design1<- MaxProLHD(sub1, pp)$Design
    design1=4*design1-2
    lhd1.ind <- nabor::knn(X, design1, k = 1)$nn.idx
    fit1 <- ssa(y~.+.^2, data=X, id.basis= lhd1.ind,
                lambda = lambda.pro, theta=theta)
    SBS_1.temp <- predict(fit1, newdata=X_test, se.fit=F)
    
    #design2<- sobol(n=sub2 ,dim = pp, init = T, scrambling = 1, seed=setseed)
    design2<- MaxProLHD(sub2, pp)$Design
    design2=4*design2-2
    lhd2.ind <- nabor::knn(X, design2, k = 1)$nn.idx
    fit2 <- ssa(y~.+.^2, data=X, id.basis= lhd2.ind,
                lambda = lambda.pro, theta=theta)
    SBS_2.temp <- predict(fit2, newdata=X_test, se.fit=F)
    
    ###############################
    # GBMED
    med.ind = MED_select(X,y,sub1)
    fit1 <- ssa(y~.+.^2, data=X, id.basis=med.ind,
                lambda = lambda.pro, theta=theta)
    med_1.temp <- predict(fit1, newdata=X_test, se.fit=F)
    
    med.ind = MED_select(X,y,sub2)
    fit2 <- ssa(y~.+.^2, data=X, id.basis=med.ind,
                lambda = lambda.pro, theta=theta)
    med_2.temp <- predict(fit2, newdata=X_test, se.fit=F)
    
    result = cbind(unif_1.temp, unif_2.temp, abs_1.temp, abs_2.temp, SBS_1.temp, SBS_2.temp,
                   med_1.temp,med_2.temp)
    result
  }
  
  UNIF1_meta[,jj] = apply(abs(aaa[,1:nloop*8-7]-as.vector(y.t)),2,mean)
  UNIF2_meta[,jj] = apply(abs(aaa[,1:nloop*8-6]-as.vector(y.t)),2,mean)
  ABS1_meta[,jj] =  apply(abs(aaa[,1:nloop*8-5]-as.vector(y.t)),2,mean)
  ABS2_meta[,jj] =  apply(abs(aaa[,1:nloop*8-4]-as.vector(y.t)),2,mean)
  SBS1_meta[,jj] = apply(abs(aaa[,1:nloop*8-3]-as.vector(y.t)),2,mean)
  SBS2_meta[,jj] = apply(abs(aaa[,1:nloop*8-2]-as.vector(y.t)),2,mean)
  MED1_meta[,jj] = apply(abs(aaa[,1:nloop*8-1]-as.vector(y.t)),2,mean)
  MED2_meta[,jj] = apply(abs(aaa[,1:nloop*8-0]-as.vector(y.t)),2,mean)
  
}

UNIF1c2=log(UNIF1_meta)
UNIF2c2=log(UNIF2_meta)
ABS1c2=log(ABS1_meta)
ABS2c2=log(ABS2_meta)
SBS1c2=log(SBS1_meta)
SBS2c2=log(SBS2_meta)
MED1c2=log(MED1_meta)
MED2c2=log(MED2_meta)

mse.matc2<-data.frame(mse=c(apply(UNIF1c2,2,mean), apply(UNIF2c2,2,mean), apply(ABS1c2,2,mean), 
                           apply(ABS2c2,2,mean), apply(SBS1c2,2,mean), apply(SBS2c2,2,mean),
                           apply(MED1c2,2,mean),apply(MED2c2,2,mean)),
                     sd=c(apply(UNIF1c2,2,sd), apply(UNIF2c2,2,sd), apply(ABS1c2,2,sd), 
                          apply(ABS2c2,2,sd), apply(SBS1c2,2,sd), apply(SBS2c2,2,sd),
                          apply(MED1c2,2,sd),apply(MED2c2,2,sd)),
                     Method=factor(rep(c("UNIF1","UNIF2","ABS1","ABS2","SBS1","SBS2","GBMED1","GBMED2"), each=length(NN))),
                     n=rep(log(NN),8))

pd <- position_dodge(0.2)
p1 = ggplot(mse.matc2,aes(x=n,y=mse,group=Method,colour=Method))
p23 = p1+theme_light()+theme(panel.grid.major = element_blank(),
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
  scale_shape_manual(values=c(1,2,1,2,1,2,1,2))+
  scale_linetype_manual(values=c(3,3,1,1,2,2,4,4))+
  scale_size_manual(values=c(1,1,1,1,1,1,1,1))+
  scale_color_manual(values=c("#4DAF4A", "#4DAF4A", "#E41A1C", "#E41A1C","#377EB8","#377EB8","#984EA3","#984EA3"))+
  xlab("log(n)")+ylab("log(MSE)")
p23

stopCluster(cl)
#300*300
