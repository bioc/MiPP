get.mipp.logistic <- function(x.train, y.train, x.test, y.test){

	y.train <- factor(y.train); levels(y.train) <- c("1","0")
	y.test <- factor(y.test); levels(y.test) <- c("1","0")
      if(is.data.frame(x.train)) x.train <- as.matrix(x.train)
      if(is.data.frame(x.test))  x.test  <- as.matrix(x.test)

	fit <- glm(y.train ~ x.train, family="binomial")

	predx <- cbind(1, x.test)%*%t(matrix(fit$coef, nrow=1))
	prob <- 1/(1+exp(-predx))
	
	postdf <- data.frame(prob, y.test)
	post.prob <- ifelse(postdf$y.test=="1", 1-postdf$prob, postdf$prob)
	ind <- ifelse(post.prob > .5, 1, 0)

	N <- length(y.test)
	nMiss <- N - sum(ind)
	Er <- nMiss/N
	MiPP <- sum(post.prob)-nMiss
	sMiPP <- MiPP/N
	
	return(list(N.Miss=nMiss, ErrorRate=Er, MiPP=MiPP, sMiPP=sMiPP))
}
