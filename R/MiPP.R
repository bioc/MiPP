
##########################################################################
#
#        
#    MiPP(Misclassification Penalized Posterior)-based Classification
#
#                            by
#
#        HyungJun Cho, Mat Soukup, and Jae K. Lee
#
#                   Version 1.2.0 (2006-03-01)   
#
##########################################################################

.First.lib <- function(lib, pkg) {  
   invisible()
   if(.Platform$OS.type=="windows" && require(Biobase) && interactive() 
   && .Platform$GUI=="Rgui") { addVigs2WinMenu("MiPP") }
}


##########################################################################
#
# Main function         
#
##########################################################################
mipp <- function(x, y, x.test=NULL, y.test=NULL, probe.ID=NULL, rule="lda", 
                 method.cut="t.test", percent.cut = 0.01, 
                 model.sMiPP.margin=0.01, min.sMiPP=0.85, n.drops=2,
                 n.fold=5, p.test=1/3, n.split=20, n.split.eval=100){

     if(is.null(probe.ID)==TRUE) probe.ID <- 1:nrow(x)
     nfold <- max(2, min(n.fold, nrow(x))) # 2 ~ N
     if(length(unique(y)) < 2)  stop("The number of classes must be >=2")
     if((length(unique(y)) > 2) & (rule !="lda") & (rule !="qda")) { 
         stop("The rule should be 'lda' or 'qda' for multi-class problem.")
     }     
     if(rule=="lda" | rule=="qda") require(MASS)
     if(rule=="svmlin" | rule=="svmrbf") require(e1071)
     cat("Please wait")


     #####################################     
     #when there is an indepedent test set
     if(is.null(x.test) == FALSE) {
        cat(".")

        #Data manipulation
        n.train.sample <- ncol(x) 
        n.test.sample  <- ncol(x.test)
        if(is.data.frame(x)==FALSE) x <- data.frame(x)
        if(is.data.frame(x.test)==FALSE) x.test <- data.frame(x.test)
        colnames(x) <- 1:ncol(x)           
        rownames(x) <- 1:nrow(x)  
        colnames(x.test) <- 1:ncol(x.test)
        rownames(x.test) <- 1:nrow(x.test)  

        x <- t(x) #convert (gene x sample) into (sample x gene)
        x.test <- t(x.test)
        y <- factor(y) 
        y.test <- factor(y.test) 

        #pre-selection
        pre.model <- "ALL"
        ii <- 1:ncol(x)
        if(percent.cut < 1) {
           ii <- pre.select(x, y, percent.cut=percent.cut)
           pre.model <- ii          
        }
        if(length(ii) < 2) stop("There are too small number of candidate genes.")

        x.tr <- x[,ii]; y.tr <- y
        x.te <- x.test[,ii]; y.te <- y.test
        out <- mipp.rule(x.train=x.tr, y.train=y.tr, x.test=x.te, y.test=y.te, 
                         nfold=nfold, min.sMiPP=min.sMiPP, n.drops=n.drops, rule=rule) 
        out[,2] <- probe.ID[ii[out[,2]]]
     
        Select <- rep(" ", nrow(out))
        i <- min(which(out$sMiPP >= max(out$sMiPP))); Select[i] <- "*"
        j <- min(which(out$sMiPP >= out$sMiPP[i]-model.sMiPP.margin)); Select[j] <- "**"
        out <- cbind(out, Select)
 
        colnames(out) <- c("Order","Gene","Tr.ER","Tr.MiPP","Tr.sMiPP","Te.ER","Te.MiPP","Te.sMiPP", "Select")
        out$Tr.ER    <- round(out$Tr.ER, 4);    out$Te.ER    <-  round(out$Te.ER, 4)
        out$Tr.MiPP  <- round(out$Tr.MiPP, 2);  out$Te.MiPP  <-  round(out$Te.MiPP, 2)
        out$Tr.sMiPP <- round(out$Tr.sMiPP, 4); out$Te.sMiPP <-  round(out$Te.sMiPP, 4)

        cat(" Done. \n")
        return(list(rule=rule, n.fold=n.fold, n.train.sample=n.train.sample, n.test.sample=n.test.sample, pre.model=pre.model, model=out)) 

     }


     #####################################
     #when there is no indepedent test set
     if(is.null(x.test) == TRUE) {


        #Data manipulation
        n.sample <- ncol(x) 
        if(is.data.frame(x)==FALSE) x <- data.frame(x)
        colnames(x) <- 1:ncol(x)           
        rownames(x) <- 1:nrow(x)  
        x <- t(x) #convert (gene x sample) into (sample x gene)
        y <- factor(y) 

        #pre-selection
        pre.model <- "ALL"
        ii <- 1:ncol(x)
        if(percent.cut < 1) {
           ii <- pre.select(x, y, percent.cut=percent.cut)
           pre.model <- ii
        }
        if(length(ii) < 2) stop("There are too small number of candidate genes.")

        x.tr <- x[,ii] 
        y.tr <- y
        out <- cv.mipp.rule(x=x.tr, y=y.tr, nfold=nfold, p.test=p.test, n.split=n.split, n.split.eval=n.split.eval,
                               model.sMiPP.margin=model.sMiPP.margin, min.sMiPP=min.sMiPP, n.drops=n.drops, rule=rule)

        out$CV.out$Gene <- probe.ID[ii[out$CV.out$Gene]]

        for(i in 1:n.split) {
            k <- ncol(out$CVCV.out)-9 ###note
            k <- max(which(!is.na(out$CVCV.out[i,1:k])))
            kk <- as.numeric(out$CVCV.out[i,2:k])
            out$CVCV.out[i,2:k] <- probe.ID[ii[kk]]
        }

        rownames(out$CV.out) <- 1:nrow(out$CV.out)
        colnames(out$CV.out) <- c("Split","Order","Gene","Tr.ER","Tr.MiPP","Tr.sMiPP","Te.ER","Te.MiPP","Te.sMiPP", "Select")
        out$CV.out$Tr.ER    <- round(out$CV.out$Tr.ER, 4);    out$CV.out$Te.ER    <-  round(out$CV.out$Te.ER, 4)
        out$CV.out$Tr.MiPP  <- round(out$CV.out$Tr.MiPP, 2);  out$CV.out$Te.MiPP  <-  round(out$CV.out$Te.MiPP, 2)
        out$CV.out$Tr.sMiPP <- round(out$CV.out$Tr.sMiPP, 4); out$CV.out$Te.sMiPP <-  round(out$CV.out$Te.sMiPP, 4)

        #n <- ncol(out$CVCV.out)
        #out$CVCV.out[,(n-5):n]    <- round(out$CVCV.out[,(n-5):n], 4)
        #out$CVCV.out[,(n-4)]    <- round(out$CVCV.out[,(n-4)], 2)

        cat("Done. \n")

        return(list(rule=rule, n.fold=n.fold, n.sample=n.sample, n.split=n.split, n.split.eval=n.split.eval, 
                    sMiPP.margin=model.sMiPP.margin, p.test=p.test,
                    pre.model=pre.model, model=out$CV.out, model.eval=out$CVCV.out)) 

     }

}



##########################################################################
#
# MiPP-based selection with cross-validation         
#
##########################################################################
cv.mipp.rule <- function(x, y, nfold, p.test, n.split, n.split.eval, 
                         model.sMiPP.margin=0.01, min.sMiPP=0, n.drops=n.drops, rule="lda") {

    n.gene <- ncol(x)
    CV.out <- data.frame(matrix(NA, n.split, 10))
    colnames(CV.out) <- c("Split","Order","Gene","Train.ErrorRate","Train.MiPP","Train.sMiPP",
                          "ErrorRate","MiPP","sMiPP","Select")
    #colnames(CV.out) <- c("Split","Order","Gene","ErrorRate","MiPP","sMiPP","Select")


    u.y <- unique(y)
    n.y <- length(u.y)

    #################################
    #Select genes from n.split splits
    gene.list <- data.frame(matrix(NA, n.split, n.gene))
    rownames(gene.list) <- paste("S",1:n.split, sep="")
    colnames(gene.list) <- paste("G",1:n.gene, sep="")

    for(iter in 1:n.split) {

        cat(".")
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

        tmp <- mipp.rule(x.train=x.train,y.train=y.train,x.test=x.test,y.test=y.test,
                         nfold=nfold, min.sMiPP=min.sMiPP, n.drops=n.drops, rule=rule)

        Split <- rep(iter, nrow(tmp))
        Select <- rep(" ", nrow(tmp))
        i <- min(which(tmp$sMiPP >= max(tmp$sMiPP))); Select[i] <- "*"
        j <- min(which(tmp$sMiPP >= tmp$sMiPP[i]-model.sMiPP.margin)); Select[j] <- "**"
        gene.list[iter,1:j] <- tmp$Gene[1:j]
        tmp <- cbind(Split, tmp, Select)
        CV.out <- rbind(CV.out, tmp)

     }
     
     tmp <- apply(gene.list, 2, is.na)
     i <- which(apply(tmp, 2, sum) >= n.split)
     gene.list <- gene.list[,-i]
     CV.out <- CV.out[-c(1:n.split),]


     ###################################
     #Evaluate optimal models of splits 
     out.Er    <- matrix(NA, n.split, n.split.eval)
     out.MiPP  <- matrix(NA, n.split, n.split.eval)
     out.sMiPP <- matrix(NA, n.split, n.split.eval)
     out2 <- data.frame(matrix(NA, n.split, 9))
     rownames(out2) <- 1:n.split 
     colnames(out2) <- c("mean ER","mean MiPP","mean sMiPP",
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

        for(jj in 1:n.split) { #Split  

            k <- max(which(!is.na(gene.list[jj,])==TRUE))
            kk <- as.numeric(gene.list[jj,1:k])
            xx.train <- x.train[,kk]
            xx.test  <- x.test[,kk]
            if(is.data.frame(xx.train)==FALSE) xx.train <- data.frame(xx.train)
            if(is.data.frame(xx.test)==FALSE) xx.test <- data.frame(xx.test)

            tmp2 <- get.mipp(xx.train, y.train, xx.test,  y.test, rule=rule)
            out.Er[jj,j]    <- tmp2$ErrorRate 
            out.MiPP[jj,j]  <- tmp2$MiPP
            out.sMiPP[jj,j] <- tmp2$sMiPP

        }  
     }
   
     out2[,1] <- apply(out.Er, 1, mean)
     out2[,2] <- apply(out.MiPP, 1, mean)
     out2[,3] <- apply(out.sMiPP, 1, mean)
     out2[,4:6] <- t(apply(out.Er, 1, quantile, probs=c(0.95, 0.50, 0.05)))
     out2[,7:9] <- t(apply(out.sMiPP, 1, quantile, probs=c(0.05, 0.50, 0.95)))

     Split <- 1:n.split
     CVCV.out <- cbind(Split, gene.list, out2)
     return(list(genes=gene.list, CV.out=CV.out, CVCV.out=CVCV.out)) 

}



##########################################################################
#
# MiPP-based selection
#
##########################################################################
mipp.rule <- function(x.train, y.train, x.test=NULL, y.test=NULL, nfold=5, min.sMiPP=0, n.drops=2, rule="lda") {
       
     n.gene <- ncol(x.train)
     n.sample.train <- nrow(x.train)
     n.sample.test  <- nrow(x.test)
     colnames(x.train) <- 1:n.gene
     colnames(x.test)  <- 1:n.gene

     tmp <- round(n.sample.train/nfold*2)
     id <- rep((1:nfold),tmp)[1:n.sample.train] #CHECK
     i <- (1:n.sample.train)[sort.list(y.train)]
     id <- id[sort.list(i)]
        
     opt.genes <- c()
     opt.Er    <- c()
     opt.MiPP  <- c()
     opt.sMiPP <- c()

     opt.Er.train    <- c()
     opt.MiPP.train  <- c()
     opt.sMiPP.train <- c()
  

     ##################
     #Pick 1-gene model
     out <- matrix(0, nfold, n.gene)
     for(i in 1:nfold) {
        y.tr <- y.train[id!=i]
        y.te <- y.train[id==i]
        for(j in 1:n.gene) {
             x.tr <- data.frame(x.train[id!=i,j])
             x.te <- data.frame(x.train[id==i,j]) 
             out[i,j] <- get.mipp(x.tr, y.tr, x.te, y.te, rule=rule)$MiPP
        }
     }

     out.sum <- apply(out, 2, sum)    
     pick.gene <- min(which(out.sum >= max(out.sum)))
     pick.gene <- as.numeric(colnames(x.train)[pick.gene])
     opt.genes <- c(opt.genes, pick.gene)

     x.train.cand <- x.train[,-opt.genes]
     x.train.opt  <- data.frame(x.train[,opt.genes])
     colnames(x.train.opt) <- opt.genes

     #Evaluate by test set
     xx.train <- data.frame(x.train[,opt.genes])
     xx.test  <- data.frame(x.test[,opt.genes])
     tmp <- get.mipp(xx.train, y.train, xx.test, y.test, rule=rule)
     opt.Er    <-c(opt.Er, tmp$ErrorRate)
     opt.MiPP  <-c(opt.MiPP, tmp$MiPP)
     opt.sMiPP <-c(opt.sMiPP, tmp$sMiPP)

     #Evaluate by train set
     tmp <- get.mipp(xx.train, y.train, xx.train, y.train, rule=rule)
     opt.Er.train    <-c(opt.Er.train, tmp$ErrorRate)
     opt.MiPP.train  <-c(opt.MiPP.train, tmp$MiPP)
     opt.sMiPP.train <-c(opt.sMiPP.train, tmp$sMiPP)


     #########################
     #Pick k-gene model (k >1)
     i.stop <- 0 
     max.sMiPP <-  opt.sMiPP
     for(jj in 2:(n.gene-1)) {
        n.gene.cand <- n.gene-jj+1
        out <- matrix(0, nfold, n.gene.cand)
        for(i in 1:nfold) {
            y.tr <- y.train[id!=i]
            y.te <- y.train[id==i]
            for(j in 1:n.gene.cand) {
                x.tr <- data.frame(x.train.opt[id!=i,], x.train.cand[id!=i,j])
                x.te <- data.frame(x.train.opt[id==i,], x.train.cand[id==i,j])
                out[i,j] <- get.mipp(x.tr,y.tr, x.te, y.te, rule=rule)$MiPP
            }
        }

        out.sum <- apply(out, 2, sum)    
        pick.gene <- min(which(out.sum >= max(out.sum)))
        pick.gene <- as.numeric(colnames(x.train.cand)[pick.gene])
        opt.genes <- c(opt.genes, pick.gene)
        x.train.opt  <- x.train[, opt.genes]
        x.train.cand <- x.train[,-opt.genes]


        #Evaluate by test set
        xx.train <- x.train[,opt.genes]
        xx.test  <- x.test[,opt.genes]
        tmp <- get.mipp(xx.train, y.train, xx.test,  y.test, rule=rule)
        opt.Er    <-c(opt.Er, tmp$ErrorRate)
        opt.MiPP  <-c(opt.MiPP, tmp$MiPP)
        opt.sMiPP <-c(opt.sMiPP, tmp$sMiPP)

        #Evaluate by train set
        tmp <- get.mipp(xx.train, y.train, xx.train,  y.train, rule=rule)
        opt.Er.train    <-c(opt.Er.train, tmp$ErrorRate)
        opt.MiPP.train  <-c(opt.MiPP.train, tmp$MiPP)
        opt.sMiPP.train <-c(opt.sMiPP.train, tmp$sMiPP)

 
        #stopping rule
        if(max.sMiPP < tmp$sMiPP) {
           max.sMiPP <- tmp$sMiPP
           i.stop <- 0
        }
        else i.stop <- i.stop + 1 

        #stop if n.drops and at least min.sMiPP
        if((i.stop >= n.drops) & (max.sMiPP >= min.sMiPP)) break 

        #stop if min number of classes is less than 4
        n.min.genes <- 2
        if(rule=="qda") n.min.genes <- 4
        if(min(unique(table(y.train))) < (jj+n.min.genes)) break         

    }


    ##################
    #Output
    i <- 1:length(opt.genes)
    final.out <- data.frame(i, opt.genes, opt.Er.train, opt.MiPP.train, opt.sMiPP.train, opt.Er, opt.MiPP, opt.sMiPP)
    colnames(final.out) <- c("Order","Gene","Train.ErrorRate","Train.MiPP","Train.sMiPP","ErrorRate","MiPP","sMiPP")

    #final.out <- data.frame(i, opt.genes, opt.Er, opt.MiPP, opt.sMiPP)
    #colnames(final.out) <- c("Order","Gene","ErrorRate","MiPP","sMiPP")

    return(final.out)
}


######END############################################################################
