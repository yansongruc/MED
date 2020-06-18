##################################################################################################################
###  This function implemetns and adaptive basis selection (ABS) method                                        ###
###                                                                                                            ###
###  --INPUT--                                                                                                 ###
###    my.y: response vector                                                                                   ###
###    nbasis: number of basis                                                                                 ###
###                                                                                                            ###
###  --OUTPUT--                                                                                                ###
###    a vector contains the "id" for the selected basis                                                       ###                 ###
##################################################################################################################


ABS_select<-function(my.y, nbasis){
  
  my.y.slice <- hist(my.y,breaks="Scott", plot=F)$breaks
  my.sample <- integer()
  
  nslice <- length(my.y.slice)-1
  nobs.slice <- floor(nbasis/nslice)
  
  for(ii in 1:(nslice)){  
    my.which <- which( my.y.slice[ii]<= my.y   &  my.y  < my.y.slice[ii+1])
    my.sample <- union(my.sample, sample(my.which, min(nobs.slice, length(my.which)) ))
  }
  
  c(my.sample, sample(setdiff(1:length(my.y),my.sample), nbasis-length(my.sample)))
}
