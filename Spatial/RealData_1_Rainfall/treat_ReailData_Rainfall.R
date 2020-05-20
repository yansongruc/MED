###### treatment of ResultN ######
MED1_N_beta1=matrix(0,6,500)
MED1_N_beta2=matrix(0,6,500)
MED1_N_beta3=matrix(0,6,500)
MED1_N_ratio=matrix(0,6,500)
MED1_N_preda=matrix(0,6,500)

MED2_N_beta1=matrix(0,6,500)
MED2_N_beta2=matrix(0,6,500)
MED2_N_beta3=matrix(0,6,500)
MED2_N_ratio=matrix(0,6,500)
MED2_N_preda=matrix(0,6,500)

MAP_N_beta1=matrix(0,6,500)
MAP_N_beta2=matrix(0,6,500)
MAP_N_beta3=matrix(0,6,500)
MAP_N_ratio=matrix(0,6,500)
MAP_N_preda=matrix(0,6,500)

DAW_N_beta1=matrix(0,6,500)
DAW_N_beta2=matrix(0,6,500)
DAW_N_beta3=matrix(0,6,500)
DAW_N_ratio=matrix(0,6,500)
DAW_N_preda=matrix(0,6,500)

R_N_beta1=matrix(0,6,500)
R_N_beta2=matrix(0,6,500)
R_N_beta3=matrix(0,6,500)
R_N_ratio=matrix(0,6,500)
R_N_preda=matrix(0,6,500)

for(i in 1:6)
{
  for(j in 1:500)
  {
    MED1_N_beta1[i,j]=ResultN[2,j][[1]][[i]][1,1]
    MED1_N_beta2[i,j]=ResultN[2,j][[1]][[i]][1,2]
    MED1_N_beta3[i,j]=ResultN[2,j][[1]][[i]][1,3]
    MED1_N_ratio[i,j]=ResultN[2,j][[1]][[i]][1,6]
    MED1_N_preda[i,j]=ResultN[2,j][[1]][[i]][1,7]
    
    MED2_N_beta1[i,j]=ResultN[2,j][[1]][[i]][5,1]
    MED2_N_beta2[i,j]=ResultN[2,j][[1]][[i]][5,2]
    MED2_N_beta3[i,j]=ResultN[2,j][[1]][[i]][5,3]
    MED2_N_ratio[i,j]=ResultN[2,j][[1]][[i]][5,6]
    MED2_N_preda[i,j]=ResultN[2,j][[1]][[i]][5,7]
    
    MAP_N_beta1[i,j]=ResultN[2,j][[1]][[i]][2,1]
    MAP_N_beta2[i,j]=ResultN[2,j][[1]][[i]][2,2]
    MAP_N_beta3[i,j]=ResultN[2,j][[1]][[i]][2,3]
    MAP_N_ratio[i,j]=ResultN[2,j][[1]][[i]][2,6]
    MAP_N_preda[i,j]=ResultN[2,j][[1]][[i]][2,7]
    
    DAW_N_beta1[i,j]=ResultN[2,j][[1]][[i]][3,1]
    DAW_N_beta2[i,j]=ResultN[2,j][[1]][[i]][3,2]
    DAW_N_beta3[i,j]=ResultN[2,j][[1]][[i]][3,3]
    DAW_N_ratio[i,j]=ResultN[2,j][[1]][[i]][3,6]
    DAW_N_preda[i,j]=ResultN[2,j][[1]][[i]][3,7]
    
    R_N_beta1[i,j]=ResultN[2,j][[1]][[i]][4,1]
    R_N_beta2[i,j]=ResultN[2,j][[1]][[i]][4,2]
    R_N_beta3[i,j]=ResultN[2,j][[1]][[i]][4,3]
    R_N_ratio[i,j]=ResultN[2,j][[1]][[i]][4,6]
    R_N_preda[i,j]=ResultN[2,j][[1]][[i]][4,7]
  }
}

MlEres_N=matrix(0,500,7)
for(i in 1:500)
{
  MlEres_N[i,]=ResultN[1,i][[1]]
}

data_N=data.frame(beta1=c(as.vector(t(MED1_N_beta1)),as.vector(t(MED2_N_beta1)),as.vector(t(MAP_N_beta1)),as.vector(t(DAW_N_beta1)),as.vector(t(R_N_beta1))),
                  beta2=c(as.vector(t(MED1_N_beta2)),as.vector(t(MED2_N_beta2)),as.vector(t(MAP_N_beta2)),as.vector(t(DAW_N_beta2)),as.vector(t(R_N_beta2))),
                  beta3=c(as.vector(t(MED1_N_beta3)),as.vector(t(MED2_N_beta3)),as.vector(t(MAP_N_beta3)),as.vector(t(DAW_N_beta3)),as.vector(t(R_N_beta3))),
                  ratio=c(as.vector(t(MED1_N_ratio)),as.vector(t(MED2_N_ratio)),as.vector(t(MAP_N_ratio)),as.vector(t(DAW_N_ratio)),as.vector(t(R_N_ratio))),
                  preda=c(as.vector(t(MED1_N_preda)),as.vector(t(MED2_N_preda)),as.vector(t(MAP_N_preda)),as.vector(t(DAW_N_preda)),as.vector(t(R_N_preda))),
                  k=factor(rep(rep(c(50,100,200,300,400,500),each=500),times=5)),
                  method=factor(rep(c("MED1","MED4","MaxPro","DaW","Random"),each=3000)))


boxplot_N_beta1=ggplot(data_N,aes(x=k,y=beta1,fill=method))+geom_boxplot()+             #boxplot
  scale_fill_brewer(palette="Pastel1")+                                         #color
  theme_light()+geom_hline(aes(yintercept=apply(MlEres_N,2,mean)[1]))+                                  #theme and line
  labs(x = "subsample size k", y = expression(widehat(beta[1])),title = expression(beta[1]))+  #labs
  theme(plot.title = element_text(hjust=0.5,size=16),legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15))                              #title location
boxplot_N_beta1

boxplot_N_beta2=ggplot(data_N,aes(x=k,y=beta2,fill=method))+geom_boxplot()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+geom_hline(aes(yintercept=apply(MlEres_N,2,mean)[2]))+
  labs(x = "subsample size k", y = expression(widehat(beta[2])),title = expression(beta[2]))+ 
  theme(plot.title = element_text(hjust=0.5,size=16),legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15))
boxplot_N_beta2

boxplot_N_beta3=ggplot(data_N,aes(x=k,y=beta3,fill=method))+geom_boxplot()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+geom_hline(aes(yintercept=apply(MlEres_N,2,mean)[3]))+
  labs(x = "subsample size k", y = expression(widehat(beta[3])),title = expression(beta[3]))+ 
  theme(plot.title = element_text(hjust=0.5,size=16),legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15))
boxplot_N_beta3

boxplot_N_ratio=ggplot(data_N,aes(x=k,y=ratio,fill=method))+geom_boxplot()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+geom_hline(aes(yintercept=apply(MlEres_N,2,mean)[6]))+
  labs(x = "subsample size k", y = expression(widehat(phi/sigma^2)),title = expression(phi/sigma^2))+ 
  theme(plot.title = element_text(hjust=0.5,size=16),legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15))
boxplot_N_ratio

boxplot_N_preda=ggplot(data_N,aes(x=k,y=preda,fill=method))+geom_boxplot()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+geom_hline(aes(yintercept=apply(MlEres_N,2,mean)[7]))+
  labs(x = "subsample size k", y = "mspe",title = "prediction accuracy")+ 
  theme(plot.title = element_text(hjust=0.5,size=16),legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15))+scale_y_log10()
boxplot_N_preda

meanpred_N=data.frame(meanpreda=c(apply(MED1_N_preda,1,mean),apply(MED2_N_preda,1,mean),apply(MAP_N_preda,1,mean),
                                  apply(DAW_N_preda,1,mean),apply(R_N_preda,1,mean)),
                      k=rep(c(50,100,200,300,400,500),times=5),
                      method=as.factor(rep(c("MED1","MED4","MaxPro","DaW","Random"),each=6)))
line_N_meanpreda=ggplot(meanpred_N,aes(x=k,y=meanpreda,color=method))+geom_line(size=1)+geom_point()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "subsample size k", y = "mean mspe",title = "prediction accuracy")+geom_hline(aes(yintercept=apply(MlEres_N,2,mean)[7]))+
  theme(plot.title = element_text(hjust=0.5,size = 16),legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15))
line_N_meanpreda



ResultNm=data.frame(k=rep(c(50,100,200,300,400,500),each=5),
                    method=as.factor(rep(c("MED1","MED4","MaxPro","DaW","Random"),times=6)))

mse_N_beta1=matrix(0,5,6)  # row:MED1,MED2,MaxPro,Deep and Wide,Random ; column: 50,100,200,300,400 and 500
for(r in 1:5)
{
  for(i in 1:6)
  {
    mse_N_beta1[r,i]=mean((data_N$beta1[(3000*(r-1)+500*(i-1)+1):(3000*(r-1)+500*i)]-apply(MlEres_N,2,mean)[1])^2)
  }
}
ResultNm$beta1=as.vector(mse_N_beta1)

mse_N_beta2=matrix(0,5,6)  # row:MED1,MED2,MaxPro,Deep and Wide,Random ; column: 30,50,100,150,200
for(r in 1:5)
{
  for(i in 1:6)
  {
    mse_N_beta2[r,i]=mean((data_N$beta2[(3000*(r-1)+500*(i-1)+1):(3000*(r-1)+500*i)]-apply(MlEres_N,2,mean)[2])^2)
  }
}
ResultNm$beta2=as.vector(mse_N_beta2)

mse_N_beta3=matrix(0,5,6)  # row:MED1,MED2,MaxPro,Deep and Wide,Random ; column: 30,50,100,150,200
for(r in 1:5)
{
  for(i in 1:6)
  {
    mse_N_beta3[r,i]=mean((data_N$beta3[(3000*(r-1)+500*(i-1)+1):(3000*(r-1)+500*i)]-apply(MlEres_N,2,mean)[3])^2)
  }
}
ResultNm$beta3=as.vector(mse_N_beta3)

mse_N_ratio=matrix(0,5,6)  # row:MED1,MED2,MaxPro,Deep and Wide,Random ; column: 30,50,100,150,200
for(r in 1:5)
{
  for(i in 1:6)
  {
    mse_N_ratio[r,i]=mean((data_N$ratio[(3000*(r-1)+500*(i-1)+1):(3000*(r-1)+500*i)]-apply(MlEres_N,2,mean)[6])^2)
  }
}
ResultNm$ratio=as.vector(mse_N_ratio)

ResultNm$betatotal=apply(ResultNm[,3:5],1,sum)

mse_N_beta=ggplot(ResultNm,aes(x=k,y=betatotal,color=method))+geom_line(size=1)+geom_point()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "subsample size k", y = "eMSE",title = "eMSE of beta")+ 
  theme(plot.title = element_text(hjust=0.5,size = 16),legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15))
mse_N_beta

mse_N_ratio=ggplot(ResultNm,aes(x=k,y=ratio,color=method))+geom_line(size=1)+geom_point()+scale_fill_brewer(palette="Pastel1")+
  theme_light()+labs(x = "subsample size k", y = "eMSE",title = "eMSE of ratio")+ 
  theme(plot.title = element_text(hjust=0.5,size=16),legend.title = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 15))
mse_N_ratio
