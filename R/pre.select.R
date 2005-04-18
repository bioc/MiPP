
pre.select <- function(x,y, percent.cut=0.01){

     tstat <- function(a, dep){    
              aaa <- t.test(a~dep, alternative = "two.sided", 
                            var.equal = FALSE)$statistic
              aaa <- abs(aaa)
              return(aaa)
     }
     
     c <- rep(0, length(y)); c[y==y[1]] <- 1 
     bbb <- apply(x, 2, tstat, dep=c)      
     q <- quantile(bbb, probs = (1-percent.cut))
     id <- which(bbb >=q)
     return(id) 
}

