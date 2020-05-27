uni_S3_pa1=matrix(0,5,100)
map_S3_pa1=matrix(0,5,100)
med_S3_pa1=matrix(0,5,100)

uni_S3_pa2=matrix(0,5,100)
map_S3_pa2=matrix(0,5,100)
med_S3_pa2=matrix(0,5,100)

for(i in 1:5){
  for(j in 1:100){
    uni_S3_pa1[i,j]=ResultS3[i,j][[1]][1,1]
    map_S3_pa1[i,j]=ResultS3[i,j][[1]][1,2]
    med_S3_pa1[i,j]=ResultS3[i,j][[1]][1,3]
    
    uni_S3_pa2[i,j]=ResultS3[i,j][[1]][2,1]
    map_S3_pa2[i,j]=ResultS3[i,j][[1]][2,2]
    med_S3_pa2[i,j]=ResultS3[i,j][[1]][2,3]
  }
}

data_S3=data.frame(pa1=c(as.vector(t(uni_S3_pa1)),as.vector(t(map_S3_pa1)),as.vector(t(med_S3_pa1))),
                   pa2=c(as.vector(t(uni_S3_pa2)),as.vector(t(map_S3_pa2)),as.vector(t(med_S3_pa2))),
                   n=rep(rep(c(2^10,2^11,2^12,2^13,2^14),each=100),times=3),
                   method=as.factor(rep(c("UNI","MaxPro","GBMED"),each=500)))
boxplot_S3_pa1=ggplot(data_S3,aes(x=as.factor(log(n)),y=log(pa1),fill=method))+geom_boxplot()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "log(n)", y = "log(MSE)")
boxplot_S3_pa1  

data_S3_mean=data.frame(pa1mean=c(apply(uni_S3_pa1,1,mean),apply(map_S3_pa1,1,mean),apply(med_S3_pa1,1,mean)),
                        pa2mean=c(apply(uni_S3_pa2,1,mean),apply(map_S3_pa2,1,mean),apply(med_S3_pa2,1,mean)),
                        n=rep(c(2^10,2^11,2^12,2^13,2^14),times=3),
                        method=as.factor(rep(c("UNI","MaxPro","GBMED"),each=5)))  
line_S3_pa1=ggplot(data_S3_mean,aes(x=log(n),y=log(pa1mean),color=method,linetype=method))+geom_line()+geom_point()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "log(n)", y = "log(MSE)")
line_S3_pa1

line_S3_pa2=ggplot(data_S3_mean,aes(x=log(n),y=log(pa2mean),color=method,linetype=method))+geom_line()+geom_point()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "log(n)", y = "log(MSE)")
line_S3_pa2



uni_S1_pa1=matrix(0,5,100)
map_S1_pa1=matrix(0,5,100)
med_S1_pa1=matrix(0,5,100)

uni_S1_pa2=matrix(0,5,100)
map_S1_pa2=matrix(0,5,100)
med_S1_pa2=matrix(0,5,100)

for(i in 1:5){
  for(j in 1:100){
    uni_S1_pa1[i,j]=ResultS1[i,j][[1]][1,1]
    map_S1_pa1[i,j]=ResultS1[i,j][[1]][1,2]
    med_S1_pa1[i,j]=ResultS1[i,j][[1]][1,3]
    
    uni_S1_pa2[i,j]=ResultS1[i,j][[1]][2,1]
    map_S1_pa2[i,j]=ResultS1[i,j][[1]][2,2]
    med_S1_pa2[i,j]=ResultS1[i,j][[1]][2,3]
  }
}

data_S1=data.frame(pa1=c(as.vector(t(uni_S1_pa1)),as.vector(t(map_S1_pa1)),as.vector(t(med_S1_pa1))),
                   pa2=c(as.vector(t(uni_S1_pa2)),as.vector(t(map_S1_pa2)),as.vector(t(med_S1_pa2))),
                   n=rep(rep(c(2^10,2^11,2^12,2^13,2^14),each=100),times=3),
                   method=as.factor(rep(c("UNI","MaxPro","GBMED"),each=500)))
boxplot_S1_pa1=ggplot(data_S1,aes(x=as.factor(log(n)),y=log(pa1),fill=method))+geom_boxplot()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "log(n)", y = "log(MSE)")
boxplot_S1_pa1  

data_S1_mean=data.frame(pa1mean=c(apply(uni_S1_pa1,1,mean),apply(map_S1_pa1,1,mean),apply(med_S1_pa1,1,mean)),
                        pa2mean=c(apply(uni_S1_pa2,1,mean),apply(map_S1_pa2,1,mean),apply(med_S1_pa2,1,mean)),
                        n=rep(c(2^10,2^11,2^12,2^13,2^14),times=3),
                        method=as.factor(rep(c("UNI","MaxPro","GBMED"),each=5)))  
line_S1_pa1=ggplot(data_S1_mean,aes(x=log(n),y=log(pa1mean),color=method,linetype=method))+geom_line()+geom_point()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "log(n)", y = "log(MSE)")
line_S1_pa1

line_S1_pa2=ggplot(data_S1_mean,aes(x=log(n),y=log(pa2mean),color=method,linetype=method))+geom_line()+geom_point()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "log(n)", y = "log(MSE)")
line_S1_pa2


uni_S2_pa=matrix(0,4,100)
map_S2_pa=matrix(0,4,100)
med_S2_pa=matrix(0,4,100)
for(i in 1:4){
  for(j in 1:100){
    uni_S2_pa[i,j]=ResultS2[i,j][[1]][1]
    map_S2_pa[i,j]=ResultS2[i,j][[1]][2]
    med_S2_pa[i,j]=ResultS2[i,j][[1]][3]
  }
}

data_S2=data.frame(pa=c(as.vector(t(uni_S2_pa)),as.vector(t(map_S2_pa)),as.vector(t(med_S2_pa))),
                   k=rep(rep(c(20,30,40,50),each=100),times=3),
                   method=as.factor(rep(c("UNI","MaxPro","GBMED"),each=400)))
data_S2_mean=data.frame(pamean=c(apply(uni_S2_pa,1,mean),apply(map_S2_pa,1,mean),apply(med_S2_pa,1,mean)),
                        k=rep(c(20,30,40,50),times=3),
                        method=as.factor(rep(c("UNI","MaxPro","GBMED"),each=4)))  
line_S2_pa=ggplot(data_S2_mean,aes(x=k,y=log(pamean),color=method,linetype=method))+geom_line()+geom_point()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "k", y = "log(MSE)")
line_S2_pa

uni_S4_pa=matrix(0,4,100)
map_S4_pa=matrix(0,4,100)
med_S4_pa=matrix(0,4,100)
for(i in 1:4){
  for(j in 1:100){
    uni_S4_pa[i,j]=ResultS4[i,j][[1]][1]
    map_S4_pa[i,j]=ResultS4[i,j][[1]][2]
    med_S4_pa[i,j]=ResultS4[i,j][[1]][3]
  }
}

data_S4=data.frame(pa=c(as.vector(t(uni_S4_pa)),as.vector(t(map_S4_pa)),as.vector(t(med_S4_pa))),
                   k=rep(rep(c(20,30,40,50),each=100),times=3),
                   method=as.factor(rep(c("UNI","MaxPro","GBMED"),each=400)))
data_S4_mean=data.frame(pamean=c(apply(uni_S4_pa,1,mean),apply(map_S4_pa,1,mean),apply(med_S4_pa,1,mean)),
                        k=rep(c(20,30,40,50),times=3),
                        method=as.factor(rep(c("UNI","MaxPro","GBMED"),each=4)))  
line_S4_pa=ggplot(data_S4_mean,aes(x=k,y=log(pamean),color=method,linetype=method))+geom_line()+geom_point()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "k", y = "log(MSE)")
line_S4_pa


uni_S5_pa=matrix(0,4,100)
map_S5_pa=matrix(0,4,100)
med_S5_pa=matrix(0,4,100)
for(i in 1:4){
  for(j in 1:100){
    uni_S5_pa[i,j]=ResultS5[i,j][[1]][1]
    map_S5_pa[i,j]=ResultS5[i,j][[1]][2]
    med_S5_pa[i,j]=ResultS5[i,j][[1]][3]
  }
}

data_S5=data.frame(pa=c(as.vector(t(uni_S5_pa)),as.vector(t(map_S5_pa)),as.vector(t(med_S5_pa))),
                   k=rep(rep(c(20,30,40,50),each=100),times=3),
                   method=as.factor(rep(c("UNI","MaxPro","GBMED"),each=400)))
data_S5_mean=data.frame(pamean=c(apply(uni_S5_pa,1,mean),apply(map_S5_pa,1,mean),apply(med_S5_pa,1,mean)),
                        k=rep(c(20,30,40,50),times=3),
                        method=as.factor(rep(c("UNI","MaxPro","GBMED"),each=4)))  
line_S5_pa=ggplot(data_S5_mean,aes(x=k,y=log(pamean),color=method,linetype=method))+geom_line()+geom_point()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "k", y = "log(MSE)")
line_S5_pa
