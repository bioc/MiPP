##########################################################################
#
# Given a gene modle, MiPP evaluation function         
#   x: row=genes (only in the given gene model), columns=samples
#   y: class of samples
#
##########################################################################
mipp.eval <- function(x, y, probe.ID=NULL, n.fold=5, p.test=1/3, n.split.eval=100, rule="lda") {


     if(length(probe.ID)==0) probe.ID <- 1:nrow(x)
     nfold <- max(2, min(n.fold, nrow(x))) # 2 ~ N
     if(rule=="lda" | rule=="qda") require(MASS)
     if(rule=="svmlin" | rule=="svmrbf") require(e1071)
     cat("Please wait...")

                                              
     #Data manipulation
     n.sample <- ncol(x) 
     if(is.data.frame(x)==FALSE) x <- data.frame(x)
     colnames(x) <- 1:ncol(x)           
     rownames(x) <- 1:nrow(x)  
     x <- t(x) #convert (gene x sample) into (sample x gene)
     y <- factor(y) 
                                                
     n.gene <- ncol(x)
     u.y <- unique(y)
     n.y <- length(u.y)

     ###################################
     #Evaluate optimal models of splits 
     out.Er    <- rep(NA, n.split.eval)
     out.MiPP  <- rep(NA, n.split.eval)
     out.sMiPP <- rep(NA, n.split.eval)
     out2 <- data.frame(matrix(NA, 9, 2))
     colnames(out2) <- c("Item","Value") 
     out2[,1] <- c("mean ER","mean MiPP","mean sMiPP",
                         "5% ER","50% ER","95% ER", 
                         "5% sMiPP","50% sMiPP","95% sMiPP")
     for(j in 1:n.split.eval) { #Splits for evaluation
        i.test  <- c()
        for(i in 1:n.y) {
            part <- sample(which(y==u.y[i])) 
            n.part <- round(length(part)*p.test)
            i.test <- c(i.test, part[1:n.part])
        }

        y.train <- y[-i.test]
        y.test  <- y[ i.test]
        x.train <- x[-i.test,]
        x.test  <- x[ i.test,]
        if(is.data.frame(x.train)==FALSE) x.train <- data.frame(x.train)
        if(is.data.frame(x.test)==FALSE) x.test <- data.frame(x.test)

        if(is.data.frame(x.train)==FALSE) x.train <- data.frame(x.train)
        if(is.data.frame(x.test)==FALSE) x.test <- data.frame(x.test)

        tmp2 <- get.mipp(x.train, y.train, x.test,  y.test, rule=rule)
        out.Er[j]    <- tmp2$ErrorRate 
        out.MiPP[j]  <- tmp2$MiPP
        out.sMiPP[j] <- tmp2$sMiPP
     }
   
     out2[1,2] <- mean(out.Er)
     out2[2,2] <- mean(out.MiPP)
     out2[3,2] <- mean(out.sMiPP)
     out2[4:6,2] <- quantile(out.Er, probs=c(0.95, 0.50, 0.05))
     out2[7:9,2] <- quantile(out.sMiPP, probs=c(0.05, 0.50, 0.95))
     
     cat("Done.")

     return(out2) 

}

######END############################################################################
