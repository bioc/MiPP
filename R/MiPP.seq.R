
##########################################################################
#
# Sequential gene selection         
#
##########################################################################
mipp.seq <- function(x, y, x.test=NULL, y.test=NULL, probe.ID=NULL, rule="lda", 
                 method.cut="t.test", percent.cut = 0.01, 
                 model.sMiPP.margin=0.01, min.sMiPP=0.85, n.drops=2,
                 n.fold=5, p.test=1/3, n.split=20, n.split.eval=100,
                 n.seq=3, cutoff.sMiPP=0.7, remove.gene.each.model="all"){

          
     p.ID <- 1:nrow(x)
     if(is.null(probe.ID) == TRUE) probe.ID <- 1:nrow(x)   
       
     #####################################     
     #when there is an indepedent test set
     if(is.null(x.test) == FALSE) {

       n.train.sample <- ncol(x) 
       n.test.sample  <- ncol(x.test)
       out2 <- data.frame(matrix(NA, 1, 9))
       colnames(out2) <- c("Order","Gene","Tr.ER","Tr.MiPP","Tr.sMiPP",
                              "Te.ER","Te.MiPP","Te.sMiPP","Select") 
     
       x.sub <- x
       x.test.sub <- x.test
       p.ID.sub <- p.ID
       Seq <- c()
       best.genes <- c()  
       for(iter in 1:n.seq) {
    
           cat("Seq "); cat(iter);  cat(" \n")
           out <- mipp(x=x.sub, y=y,  x.test = x.test.sub, y.test = y.test, probe.ID=p.ID.sub, rule=rule, 
                 method.cut=method.cut, percent.cut = percent.cut, 
                 model.sMiPP.margin=model.sMiPP.margin, min.sMiPP=min.sMiPP, n.drops=n.drops,
                 n.fold=n.fold)
                 
           #ADD SEQ                 
           out2 <- rbind(out2, out$model)
           Seq <- c(Seq, rep(iter, nrow(out$model)))         
     
           k <- which(out$model$Select=="**")
           nc <- ifelse(remove.gene.each.model=="first", 1, k)
           best.genes <- sort(unique(c(best.genes, out$model$Gene[1:nc])))  
           if(length(best.genes) < nrow(x)) { 
               x.sub <- x[-best.genes,]
               x.test.sub <- x.test[-best.genes,]
               p.ID.sub <- p.ID[-best.genes]
           }
           
           if((iter < n.seq) & (nrow(x)-length(best.genes) <= 2)) {
              cat("Two or less genes were left for the next Seq runs, so your run stopped after Seq", iter)
              break
           }    


       }           

       ###GENE ID
       out2 <- out2[-1,]
       out2$Gene <- probe.ID[out2$Gene]
       out2 <- cbind(Seq, out2)
       rownames(out2) <- 1:nrow(out2)
       if(length(best.genes) > 0) best.genes <- probe.ID[best.genes]
    
       return(list(rule=rule, n.fold=n.fold, n.train.sample=n.train.sample, n.test.sample=n.test.sample, model=out2, genes.selected=best.genes)) 

     }


     #####################################
     #when there is no indepedent test set
     if(is.null(x.test) == TRUE) {

       n.sample <- ncol(x) 
       CV.out2 <- data.frame(matrix(NA, 1, 10))
       colnames(CV.out2) <- c("Split", "Order","Gene","Tr.ER","Tr.MiPP","Tr.sMiPP",
                              "Te.ER","Te.MiPP","Te.sMiPP","Select") 
       CVCV.out2 <- data.frame(matrix(NA, n.seq*n.split, 1+n.sample+9))
       colnames(CVCV.out2) <- c("Split", paste("G", 1:n.sample, sep = ""), "mean ER","mean MiPP","mean sMiPP",
                         "5% ER","50% ER","95% ER", "5% sMiPP","50% sMiPP","95% sMiPP")

       x.sub <- x
       p.ID.sub <- p.ID
       Seq <- c()
       best.genes <- c()
       for(iter in 1:n.seq) {
       
           cat("Seq "); cat(iter);  cat(" \n")
           out <- mipp(x=x.sub, y=y, probe.ID=p.ID.sub, rule=rule, 
                 method.cut=method.cut, percent.cut = percent.cut, 
                 model.sMiPP.margin=model.sMiPP.margin, min.sMiPP=min.sMiPP, n.drops=n.drops,
                 n.fold=n.fold, p.test=p.test, n.split=n.split, n.split.eval=n.split.eval)
                 

           #ADD SEQ                 
           CV.out2 <- rbind(CV.out2, out$model)
           Seq <- c(Seq, rep(iter, nrow(out$model)))  

           nc <- ncol(out$model.eval)
           nc2 <-(iter-1)*n.split
           CVCV.out2[(nc2+1):(nc2+n.split),1:(nc-9)] <- out$model.eval[,1:(nc-9)]                
           CVCV.out2[(nc2+1):(nc2+n.split),(n.sample+2):(n.sample+10)] <- out$model.eval[,(nc-8):nc]               
         
           
           nc <- ifelse(remove.gene.each.model=="first", 2, (n.sample+1))
           k <- which(CVCV.out2[,(1+n.sample+7)] >= cutoff.sMiPP)
           if(length(k) > 0) {
              best.genes <- sort(unique(as.numeric(na.omit(as.vector(as.matrix(CVCV.out2[k,2:nc]))))))           
              if(length(best.genes) < nrow(x)) { 
                 x.sub <- x[-best.genes,]
                 p.ID.sub <- p.ID[-best.genes]
              }
           }
           
           if((iter < n.seq) & (nrow(x)-length(best.genes) <= 2)) {
              cat("Warning: Two or less genes were left for the next Seq runs, so your run stopped after Seq", iter)
              break
           }    

       }           

 

       ###GENE ID
       CV.out2$Gene <- probe.ID[CV.out2$Gene]
       kk <- as.numeric(as.vector(as.matrix(CVCV.out2[,2:(n.sample+1)])))           
       CVCV.out2[,2:(n.sample+1)] <- probe.ID[kk]
       
       #Remove missing columns and add seq
       CV.out2 <- CV.out2[-1,]
       CV.out2 <- cbind(Seq, CV.out2)
       rownames(CV.out2) <- 1:nrow(CV.out2)

       k <- max(apply(!is.na(CVCV.out2), 1, sum))-9        
       CVCV.out2 <- CVCV.out2[,-c((k+1):(n.sample+1))]
       
       Seq <- c() 
       for(i in 1:n.seq) Seq <- c(Seq, rep(i, n.split)) 
       CVCV.out2 <- cbind(Seq, CVCV.out2)
       if(length(best.genes) > 0) best.genes <- probe.ID[best.genes]

       return(list(rule=rule, n.fold=n.fold, n.sample=n.sample, n.split=n.split, n.split.eval=n.split.eval, 
                    sMiPP.margin=model.sMiPP.margin, p.test=p.test,
                    model=CV.out2, model.eval=CVCV.out2, genes.selected=best.genes)) 

     }

     cat(" Done. \n")

}



######END############################################################################
