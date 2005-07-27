# Linear Kernel Decision Function
linearkernel.decision.function <-function(newx, oldx, svmobj) {
    # oldx is the original training data matrix
    # svmobj is the name of the svm object

    # Extract y*alpha:
    	svcoefs <- svmobj$coefs
    # Extract b:
    	svconstant <- -1*svmobj$rho   
    # Get the support vectors
    	svdata <- oldx[svmobj$index,]
    # Reformat the new x
    	xt <- newx
    	nrowxt <- length(oldx[1,])
    	dim(xt) <- c(nrowxt,1)   
    # linear kernel:
    	prods <- svdata %*% xt   
    # compute h(x):
    	h <- t(prods) %*% svcoefs 
    # compute f(x):
    	#h + svconstant    
    	return(h + svconstant)    
}


get.mipp.svm.linear <- function(x.train, y.train, x.test, y.test){

        x.train <- as.matrix(x.train)
        x.test  <- as.matrix(x.test)

	y.train <- factor(y.train)
	y.test <- factor(y.test)

	fit <- svm(x.train, y.train, kernel="linear")

	True.class <- y.test
	Pred.class <- predict(fit, x.test)

	fofx <- numeric(length(y.test))
	for(i in 1:length(y.test)){
		xin <- x.test[i,]
		fofx[i] <- linearkernel.decision.function(xin, x.train, fit)
	}

      c <-1 #optimal parameter
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

