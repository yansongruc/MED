library(mvtnorm)
library(MASS)
library(ggplot2)
library(gss)
library(MaxPro)
library(randtoolbox)
library(np)
library(mined)

c1="white"
c2="gray10"
pp=2
SNR=5
mtype=list("cubic",c(0,1))

copula=function(t)
{
  t1=t[1];t2=t[2]
  f1=sign(t1)*t1^2
  f2=sign(t2)*(abs(t2))^3
  Sigma=0.5*matrix(1,2,2)+0.5*diag(1,2)
  res=dmvnorm(c(f1,f2),c(0,0),Sigma)*6*abs(t1)*t2*t2
  return(res)
}

N=5000
sub2=floor(5*N^(2/9))    #sub2=33
#Design matrix
set.seed(100)
X<-as.data.frame(sobol(n=N ,dim = pp, init = T, scrambling = 1, seed=100))
X=4*X-2
y.true=apply(X,1,copula)
ssig=sd(y.true)
y <- y.true + rnorm(N, 0, ssig/SNR)    

sam.size <- ceiling(50*(N^{1/4}))
esps.re <- esps(y~.+.^2, data=data.frame(y, X), sam.size,r=3)
p <- esps.re$p
pro.gcv.lamb <- esps.re$lambda*((N/sam.size)^(-r/(p*r+1)))
lambda.pro <- log10(N*pro.gcv.lamb)
theta <- esps.re$theta

######################################################################
##### True surface
xs=seq(-2,2,length.out = 100)
df=data.frame(V1=rep(xs,each=100),V2=rep(xs,times=100))
df$y=apply(df,1,copula)

##############################
# Unif
set.seed(55)
unif.ind = sample.int(N,size=sub2,replace=F)
fit2 <- ssa(y~.+.^2, data=X, id.basis=unif.ind,
            lambda = lambda.pro, theta=theta)
unif_2.temp <- predict(fit2, newdata=df[,1:2], se.fit=F)
df$uniy=unif_2.temp

###############################
# ABS
set.seed(100)
abs.ind = ABS_select(y,sub2)
fit2 <- ssa(y~.+.^2, data=X, id.basis=abs.ind,
            lambda = lambda.pro, theta=theta)
abs_2.temp <- predict(fit2, newdata=df[,1:2], se.fit=F)
df$absy=abs_2.temp

###############################
# SBS
set.seed(786)
design2<- MaxProLHD(sub2, pp)$Design
design2=4*design2-2
lhd2.ind <- nabor::knn(X, design2, k = 1)$nn.idx
fit2 <- ssa(y~.+.^2, data=X, id.basis= lhd2.ind,
            lambda = lambda.pro, theta=theta)
SBS_2.temp <- predict(fit2, newdata=df[,1:2], se.fit=F)
df$sbsy=SBS_2.temp

###############################
# GBMED
set.seed(100)
med.ind = MED_select(X,y,sub2)
fit2 <- ssa(y~.+.^2, data=X, id.basis=med.ind,
            lambda = lambda.pro, theta=theta)

med_2.temp <- predict(fit2, newdata=df[,1:2], se.fit=F)
df$medy=med_2.temp

######################################################################
##### plot figures
df1=df
df1[1,3:7]=rep(-0.2,5)
df1[7475,3:7]=rep(max(df$y),5)   #rescale color

###############################
# True surface
p=ggplot(df1,aes(V1,V2))+geom_raster(aes(fill=y))
p1=p+theme_bw()+theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_rect(colour = "black"),
                      axis.text=element_text(size=12),
                      axis.title=element_text(size=14),
                      legend.position = "none",
                      legend.title = element_blank(),
                      legend.key.width=unit(1,"line"),
                      legend.key.height=unit(1,"line"))+
  geom_raster(data=df,aes(V1,V2,fill=y))+
  scale_fill_gradient(low = c1, high = c2)+
  xlab(expression(x[1]))+ylab(expression(x[2]))
p1

###############################
# UNI
p <- ggplot(df1, aes(V1,V2))+geom_raster(aes(fill=uniy))
p2<-p +theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"),
        plot.title = element_text(size=10, hjust = 0.5),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.width=unit(0.5,"line"),
        legend.key.height=unit(1,"line"))+
  scale_fill_gradient(low = c1, high = c2)+
  xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_raster(data = df,aes(V1,V2,fill=uniy))+
  geom_point(data=X[unif.ind,],aes(V1,V2),col="red")
p2

###############################
# ABS
p <- ggplot(df1, aes(V1,V2))+geom_raster(aes(fill=absy))
p3<-p +theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"),
        plot.title = element_text(size=10, hjust = 0.5),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.width=unit(0.5,"line"),
        legend.key.height=unit(1,"line"))+
  scale_fill_gradient(low = c1, high = c2)+
  xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_raster(data = df,aes(V1,V2,fill=absy))+
  geom_point(data=X[abs.ind,],aes(V1,V2),col="red")
p3

###############################
# SBS
p <- ggplot(df1, aes(V1,V2))+geom_raster(aes(fill=sbsy))
p4<-p +theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"),
        plot.title = element_text(size=10, hjust = 0.5),
        legend.position = "none",
        legend.title = element_blank(),
        legend.key.width=unit(0.5,"line"),
        legend.key.height=unit(1,"line"))+
  scale_fill_gradient(low = c1, high = c2)+
  xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_raster(data = df,aes(V1,V2,fill=sbsy))+
  geom_point(data=X[lhd2.ind,],aes(V1,V2),col="red")
p4

###############################
# GBMED
p <- ggplot(df1, aes(V1,V2))+geom_raster(aes(fill=medy))
p5<-p +theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black"),
        plot.title = element_text(size=10, hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank(),
        legend.key.width=unit(0.5,"line"),
        legend.key.height=unit(1,"line"))+
  scale_fill_gradient(low = c1, high = c2)+
  xlab(expression(x[1]))+ylab(expression(x[2]))+
  geom_raster(data = df,aes(V1,V2,fill=medy))+
  geom_point(data=X[med.ind,],aes(V1,V2),col="red")
p5

