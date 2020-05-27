esps = function(form, data, sam.size=NULL, r=NULL, iter.c=10, iter.p=5){
  mf = model.frame(form, data = data)
  sam.theta <- NULL
  sam.lamb <- NULL
  obs = dim(data)[1]
  for(t in 1:iter.c){
    set.seed(t)
    sam.indx <- sample(obs, sam.size)
    sam.dat <- data[sam.indx, ]
    sam.fit <- ssanova(form, data=sam.dat, seed=t, alpha=1.0)
    sam.theta <- rbind(sam.theta, sam.fit$theta)
    sam.lamb[t] <- sam.fit$nlambda
  }
  me.lamb <- median(sam.lamb)
  me.theta <- apply(sam.theta, 2, median)
  me.lamb.org = 10^me.lamb/sam.size
  if(is.null(r)){
    cat("r=4 for the univariate case and r=3 for the multivariate case!", "\n")
  }
  #choose the best p

  
  # resp = model.extract(mf, "response")
  # p.score <- matrix(nrow=iter.p, ncol=2)
  # p.size <- 5*sam.size
  # for(pp in 1:iter.p){
  #   df <- as.matrix(cut(resp, hist(resp,plot=F)$breaks))
  #   colnames(df) <- "invl"
  #   sam.indx.p <- strata(df, levels(df), p.size, method="srswor")[,1]
  #   sam.p <- data[sam.indx.p,]
  #   for(p in c(1:2)){
  #     pro.lamb.p <- me.lamb.org*((p.size/sam.size)^(-r/(p*r+1)))
  #     lambda.p <- log10(p.size*pro.lamb.p)
  #     theta.p <- me.theta
  #     sam.fit.p <-  ssa(form, data=sam.p, lambda = lambda.p, theta=theta.p, alpha=1.0)
  #     p.score[pp, p] <- sam.fit.p$score
  #   }
  # }
  # p.c <- as.numeric(names(which.max(table(apply(p.score, 1, which.min)))))
  p.c=1
list(lambda=me.lamb.org, theta=me.theta, p=p.c)
}
