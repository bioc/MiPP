#Computing Miss Error and MiPP after LDA
get.mipp.lda <- function(x.train, y.train, x.test, y.test){
     
     y <- y.train
     dat.train <- cbind(x.train, y)
     y <- y.test
     dat.test  <- cbind(x.test, y)

     if(is.data.frame(dat.train)==FALSE) dat.train <- data.frame(dat.train)
     if(is.data.frame(dat.test)==FALSE) dat.test <- data.frame(dat.test)
     colnames(dat.train) <- c(1:ncol(x.train), "y")
     colnames(dat.test)  <- c(1:ncol(x.test), "y")

     fit <- lda(y ~ ., dat.train)
     out <- predict(fit, dat.test)

     u.class <- unique(colnames(out$post))
     n.class <- length(u.class)

     True.class <- dat.test$y
     Pred.class <- out$class

     post.prob <-0
     for(j in 1:n.class) {
         i <- which(True.class == u.class[j]) 
         post.prob <- post.prob + sum(out$post[i,j,drop=FALSE])
     }

     N <- length(True.class) 
     nMiss <- N- sum(True.class == Pred.class) 
     Er <- nMiss/nrow(dat.test)
     MiPP <- post.prob - nMiss
     sMiPP <- MiPP/N
     return(list(N.Miss=nMiss, ErrorRate=Er, MiPP=MiPP, sMiPP=sMiPP))
}


