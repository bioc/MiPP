#Choose rule
get.mipp <- function(x.train, y.train, x.test, y.test, rule){

     if(rule !="lda" & rule !="qda" & rule !="logistic" & rule !="svmlin" & rule !="svmrbf") 
     stop("No rule: ", rule)
     
     if(rule=="lda")      tmp <- get.mipp.lda(x.train, y.train, x.test,  y.test)
     if(rule=="qda")      tmp <- get.mipp.qda(x.train, y.train, x.test,  y.test)
     if(rule=="logistic") tmp <- get.mipp.logistic(x.train, y.train, x.test,  y.test)
     if(rule=="svmlin")   tmp <- get.mipp.svm.linear(x.train, y.train, x.test,  y.test)
     if(rule=="svmrbf")   tmp <- get.mipp.svm.rbf(x.train, y.train, x.test,  y.test)
     
     return(list(N.Miss=tmp$N.Miss, ErrorRate=tmp$ErrorRate, MiPP=tmp$MiPP, sMiPP=tmp$sMiPP))
}


