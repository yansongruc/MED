#output treatment
#Case 1
MED_Li1_beta1=matrix(0,6,500)
MED_Li1_beta2=matrix(0,6,500)
MED_Li1_beta3=matrix(0,6,500)
MED_Li1_pa=matrix(0,6,500)

IBO_Li1_beta1=matrix(0,6,500)
IBO_Li1_beta2=matrix(0,6,500)
IBO_Li1_beta3=matrix(0,6,500)
IBO_Li1_pa=matrix(0,6,500)

UNI_Li1_beta1=matrix(0,6,500)
UNI_Li1_beta2=matrix(0,6,500)
UNI_Li1_beta3=matrix(0,6,500)
UNI_Li1_pa=matrix(0,6,500)

for(i in 1:6)
{
  for(j in 1:500)
  {
    MED_Li1_beta1[i,j]=ResultLi1[2,j][[1]][[i]][1,1]
    MED_Li1_beta2[i,j]=ResultLi1[2,j][[1]][[i]][1,2]
    MED_Li1_beta3[i,j]=ResultLi1[2,j][[1]][[i]][1,3]
    MED_Li1_pa[i,j]=ResultLi1[2,j][[1]][[i]][1,4]
    
    IBO_Li1_beta1[i,j]=ResultLi1[2,j][[1]][[i]][2,1]
    IBO_Li1_beta2[i,j]=ResultLi1[2,j][[1]][[i]][2,2]
    IBO_Li1_beta3[i,j]=ResultLi1[2,j][[1]][[i]][2,3]
    IBO_Li1_pa[i,j]=ResultLi1[2,j][[1]][[i]][2,4]
    
    UNI_Li1_beta1[i,j]=ResultLi1[2,j][[1]][[i]][3,1]
    UNI_Li1_beta2[i,j]=ResultLi1[2,j][[1]][[i]][3,2]
    UNI_Li1_beta3[i,j]=ResultLi1[2,j][[1]][[i]][3,3]
    UNI_Li1_pa[i,j]=ResultLi1[2,j][[1]][[i]][3,4]
  }
}

data_Li1=data.frame(beta1=c(as.vector(t(MED_Li1_beta1)),as.vector(t(IBO_Li1_beta1)),as.vector(t(UNI_Li1_beta1))),
                    beta2=c(as.vector(t(MED_Li1_beta2)),as.vector(t(IBO_Li1_beta2)),as.vector(t(UNI_Li1_beta2))),
                    beta3=c(as.vector(t(MED_Li1_beta3)),as.vector(t(IBO_Li1_beta3)),as.vector(t(UNI_Li1_beta3))),
                    pa=c(as.vector(t(MED_Li1_pa)),as.vector(t(IBO_Li1_pa)),as.vector(t(UNI_Li1_pa))),
                    k=as.factor(rep(rep(c(100,200,400,600,800,1000),each=500),times=3)),
                    method=as.factor(rep(c("MED","IBOSS","UNI"),each=3000)))

data_Li1_mse=data.frame(beta1=c(apply((MED_Li1_beta1-1)^2,1,mean),apply((IBO_Li1_beta1-1)^2,1,mean),apply((UNI_Li1_beta1-1)^2,1,mean)),
                        beta2=c(apply((MED_Li1_beta2-1)^2,1,mean),apply((IBO_Li1_beta2-1)^2,1,mean),apply((UNI_Li1_beta2-1)^2,1,mean)),
                        beta3=c(apply((MED_Li1_beta3-1)^2,1,mean),apply((IBO_Li1_beta3-1)^2,1,mean),apply((UNI_Li1_beta3-1)^2,1,mean)),
                        k=rep(c(100,200,400,600,800,1000),times=3),
                        method=as.factor(rep(c("MED","IBOSS","UNI"),each=6))) 
data_Li1_mse$slope=data_Li1_mse$beta2+data_Li1_mse$beta3
data_Li1_mse$pa=c(apply(MED_Li1_pa,1,mean),apply(IBO_Li1_pa,1,mean),apply(UNI_Li1_pa,1,mean)) #only mean

Full_Li1=matrix(0,500,4)
for(i in 1:500)
{
  Full_Li1[i,]=ResultLi1[1,i][[1]]
}

beta1_Li1_mse=ggplot(data_Li1_mse,aes(x=k,y=beta1,color=method))+geom_line()+geom_point()+
 scale_fill_brewer(palette="Accent")+
  theme_light()+labs(x = "subsample size k", y = "eMSE",title = expression(beta[0]))+ 
  theme(plot.title = element_text(hjust=0.5))
beta1_Li1_mse

slope_Li1_mse=ggplot(data_Li1_mse,aes(x=k,y=slope,color=method))+geom_line()+geom_point()+
  scale_fill_brewer(palette="Accent")+
  theme_light()+labs(x = "subsample size k", y = "eMSE",title = "slope")+ 
  theme(plot.title = element_text(hjust=0.5))
slope_Li1_mse

pa_Li1_mse=ggplot(data_Li1_mse,aes(x=k,y=pa,color=method))+geom_line()+geom_point()+
  scale_fill_brewer(palette="Accent")+geom_hline(aes(yintercept=mean(Full_Li1[,4])))+
  theme_light()+labs(x = "subsample size k", y = "eMSE",title = "prediction accuracy")+ 
  theme(plot.title = element_text(hjust=0.5))
pa_Li1_mse

