get.mipp.svm.rbf <- function(x.train, y.train, x.test, y.test){

      if(is.data.frame(x.train)) x.train <- as.matrix(x.train)
      if(is.data.frame(x.test))  x.test  <- as.matrix(x.test)
	y.train <- factor(y.train)
	y.test <- factor(y.test)

	gammap <- 1/length(ncol(x.train))
	fit <- svm(x.train, y.train, kernel="radial", gamma=gammap)
	
	True.class <- y.test
	Pred.class <- predict(fit, x.test)

	fofx <- numeric(length(y.test))
	for(i in 1:length(y.test)){
		xin <- x.test[i,]
		fofx[i] <- rbfkernel.decision.function(xin, x.train, fit)
	}

      c <- 100 #optimal parameter?
	prob <- 1/(1+c*exp(-fofx))
	postdf <- data.frame(prob, True.class)
	post.prob <- ifelse(postdf$True.class==Pred.class, 1-postdf$prob, postdf$prob)

	N <- length(y.test)
	nMiss <- N - sum(True.class==Pred.class)
	Er <- nMiss/N
	MiPP <- sum(post.prob)-nMiss
	sMiPP <- MiPP/N

	return(list(N.Miss=nMiss, ErrorRate=Er, MiPP=MiPP, sMiPP=sMiPP))
}


rbfkernel.decision.function <- function(newx, oldx, svmobj) {
    # oldx is the original training data matrix
    # svmobj is the name of the svm object
    
    # Extract y*alpha:
    	svcoefs <- svmobj$coefs
    # Extract b:
    	svconstant <- -1*svmobj$rho
    # Extract gamma:
    	svgamma <- svmobj$gamma
    # Get the support vectors
    	svdata <- oldx[svmobj$index,]
    # How many support vectors?
    	numsv <- length(svmobj$index)
    # reformat newx
    	p <- length(oldx[1,])
    	xt <- matrix(0, nrow=numsv, ncol=p)
    	for(i in 1:p){
        	xt[,i] <- rep(newx[i], numsv)
    	}     
    # rbf kernel:
    	difs <- (svdata - xt)
    	difs2 <- apply(difs, 2, function(x)x^2)
    	difs3 <- apply(difs2, 1, sum)
    	ks <- exp(-1*svgamma*difs3) 
    # compute h(x):
    	h <- t(ks) %*% svcoefs
    # compute f(x):
    	#h + svconstant
    	return(h + svconstant)
}

