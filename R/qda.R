#Computing Miss Error and MiPP after QDA
get.mipp.qda <- function(x.train, y.train, x.test, y.test){
     

     colnames(x.train) <- c(1:ncol(x.train))
     colnames(x.test)  <- c(1:ncol(x.test))

     fit <- qda(x.train, y.train)
     out <- predict(fit, x.test)

     u.class <- unique(colnames(out$post))
     n.class <- length(u.class)

     True.class <- y.test
     Pred.class <- out$class

     post.prob <-0
     for(j in 1:n.class) {
         i <- which(True.class == u.class[j]) 
         post.prob <- post.prob + sum(out$post[i,j])
     }

     N <- length(True.class) 
     nMiss <- N- sum(True.class == Pred.class) 
     Er <- nMiss/nrow(x.test)
     MiPP <- post.prob - nMiss
     sMiPP <- MiPP/N

     return(list(N.Miss=nMiss, ErrorRate=Er, MiPP=MiPP, sMiPP=sMiPP))
}


