library(mvtnorm)
library(MASS)
library(randtoolbox)
library(gss)
library(foreach)
library(doParallel)
library(ggplot2)
#source("ssa.R")
#source("esps-m.R")
<<<<<<< HEAD
=======
#source("ABS_select.R")
>>>>>>> f7fcb26... Splines
#source("MED_select.R")

################################################################################
# Setting 1: function g in Meng's first setting                                #
#            Sobol sequence locations                                          #
<<<<<<< HEAD
#            Full sample sizes are c(2^10, 2^11, 2^12, 2^13, 2^14)             #
=======
#            Full sample sizes are c(2^12, 2^13, 2^14, 2^15, 2^16)             #
>>>>>>> f7fcb26... Splines
#            Subsample sizes are k1=10*N^(1/9) and k2=5*N^(2/9)                #
################################################################################
set.seed(100)
nloop=100
<<<<<<< HEAD
NN=2^(seq(10,14,1))
ntest=5000
SNR=2                         #signal noise ratio 5 or 2
=======
NN=2^(seq(12,16,1))
ntest=5000
SNR=5                         #signal noise ratio 5 or 2
>>>>>>> f7fcb26... Splines
sub1_meta=round(10*NN^(1/9))
#sub1_meta=c(22,23,25,27,30)
sub2_meta=round(5*NN^(2/9))
#sub2_meta=c(23,27,32,37,43)

# function g
pp=2
g1 = function(x) x
g2 = function(x) (2*x-1)^2
g3 = function(x) sin(2*pi*x)/(2-sin(2*pi*x))
g4 = function(x) 0.1*sin(2*pi*x)+0.2*cos(2*pi*x)+0.3*sin(2*pi*x)^2+0.4*cos(2*pi*x)^3+0.5*sin(2*pi*x)^3
gene_y = function(x){
  x1=x[,1];x2=x[,2]
  g1(x1*x2)+g2(x2)+g3(x1)+g4(x2)+g3((x1+x2)/2)
}

# Generate testing data: sobol; on [0,1]^2
<<<<<<< HEAD
X_test<-as.data.frame(sobol(n=ntest ,dim = pp, init = T, scrambling = 1, seed=100))
y.t<-gene_y(X_test)

# using a part of data to find suitable p, lambda, theta by esps-m.R
sam.size <- ceiling(50*(NN[5]^{1/4}))
X<-as.data.frame(sobol(n=NN[5] ,dim = pp, init = T, scrambling = 1, seed=100))
=======
X_test<-as.data.frame(sobol(n=ntest ,dim = pp))
y.t<-gene_y(X_test)

# using a part of data to find suitable p, lambda, theta by esps-m.R
set.seed(100)
sam.size <- ceiling(50*(NN[5]^{1/4}))
X<-as.data.frame(sobol(n=NN[5] ,dim = pp))
>>>>>>> f7fcb26... Splines
y<-gene_y(X)
esps.re <- esps(y~.+.^2, data=data.frame(y, X), sam.size,r=3)

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
<<<<<<< HEAD
  X<-as.data.frame(sobol(n=N ,dim = pp, init = T, scrambling = 1, seed=100))
=======
  X<-as.data.frame(sobol(n=N ,dim = pp))
>>>>>>> f7fcb26... Splines
  y.true<-gene_y(X)
  
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
    lhd1.ind <- nabor::knn(X, design1, k = 1)$nn.idx
    fit1 <- ssa(y~.+.^2, data=X, id.basis= lhd1.ind,
                lambda = lambda.pro, theta=theta)
    SBS_1.temp <- predict(fit1, newdata=X_test, se.fit=F)
    
    #design2<- sobol(n=sub2 ,dim = pp, init = T, scrambling = 1, seed=setseed)
    design2<- MaxProLHD(sub2, pp)$Design
    lhd2.ind <- nabor::knn(X, design2, k = 1)$nn.idx
    fit2 <- ssa(y~.+.^2, data=X, id.basis= lhd2.ind,
                lambda = lambda.pro, theta=theta)
    SBS_2.temp <- predict(fit2, newdata=X_test, se.fit=F)

    ###############################3
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

<<<<<<< HEAD
UNIF12=log(UNIF1_meta)
UNIF22=log(UNIF2_meta)
ABS12=log(ABS1_meta)
ABS22=log(ABS2_meta)
SBS12=log(SBS1_meta)
SBS22=log(SBS2_meta)
MED12=log(MED1_meta)
MED22=log(MED2_meta)

mse.mat2<-data.frame(mse=c(apply(UNIF12,2,mean), apply(UNIF22,2,mean), apply(ABS12,2,mean), 
                          apply(ABS22,2,mean), apply(SBS12,2,mean), apply(SBS22,2,mean),
                          apply(MED12,2,mean),apply(MED22,2,mean)),
                    sd=c(apply(UNIF12,2,sd), apply(UNIF22,2,sd), apply(ABS12,2,sd), 
                         apply(ABS22,2,sd), apply(SBS12,2,sd), apply(SBS22,2,sd),
                         apply(MED12,2,sd),apply(MED22,2,sd)),
=======
UNIF1=log(UNIF1_meta)
UNIF2=log(UNIF2_meta)
ABS1=log(ABS1_meta)
ABS2=log(ABS2_meta)
SBS1=log(SBS1_meta)
SBS2=log(SBS2_meta)
MED1=log(MED1_meta)
MED2=log(MED2_meta)

mse.mat<-data.frame(mse=c(apply(UNIF1,2,mean), apply(UNIF2,2,mean), apply(ABS1,2,mean), 
                          apply(ABS2,2,mean), apply(SBS1,2,mean), apply(SBS2,2,mean),
                          apply(MED1,2,mean),apply(MED2,2,mean)),
                    sd=c(apply(UNIF1,2,sd), apply(UNIF2,2,sd), apply(ABS1,2,sd), 
                         apply(ABS2,2,sd), apply(SBS1,2,sd), apply(SBS2,2,sd),
                         apply(MED1,2,sd),apply(MED2,2,sd)),
>>>>>>> f7fcb26... Splines
                    Method=factor(rep(c("UNIF1","UNIF2","ABS1","ABS2","SBS1","SBS2","GBMED1","GBMED2"), each=length(NN))),
                    n=rep(log(NN),8))

pd <- position_dodge(0.2)
<<<<<<< HEAD
p1 = ggplot(mse.mat2,aes(x=n,y=mse,group=Method,colour=Method))
=======
p1 = ggplot(mse.mat,aes(x=n,y=mse,group=Method,colour=Method))
>>>>>>> f7fcb26... Splines
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
