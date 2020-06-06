##################################################################################################################
###     GBMED subsampling method                                                                               ###
###                                                                                                            ###
###  --INPUT--                                                                                                 ###
###    daTa: location X and response y                                                                         ###
###    n: subsample size                                                                                       ###
###                                                                                                            ###
###  --OUTPUT--                                                                                                ###
###    a vector contains the "id" for the selected basis                                                       ###  
###  
##################################################################################################################
library("np")
library("mined")
MED_select=function(X,y,n)
{
  Data=data.frame(loc.x1=X[,1],loc.x2=X[,2],y=y)
  if(nrow(X)>1000){
    ib=sample(1:nrow(X),1e3)
    Dataib=Data[ib,]
    bw=npregbw(Dataib$y~Dataib$loc.x1+Dataib$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  }
  else
    bw=npregbw(Data$y~Data$loc.x1+Data$loc.x2,regtype="ll",bwmethod="cv.aic")$bw
  
  regres=npreg(Data$y~Data$loc.x1+Data$loc.x2,bws=bw,gradients=TRUE)
  grad=regres$grad
  cgrad=sqrt(grad[,1]^2+grad[,2]^2)
  
  findi=function(x){idex=which(Data$loc.x1==x[1] & Data$loc.x2==x[2])}
  MEDres=mined::SelectMinED(cbind(Data$loc.x1,Data$loc.x2),log(cgrad),n,1,2)
  MEDloc=MEDres$points
  iMED=apply(MEDloc,1,findi)
  
  return(iMED)
}
