
pre.select <- function(x, y, percent.cut=0.01){

     tstat <- function(a, dep){    
              aaa <- t.test(a~dep, alternative = "two.sided", 
                            var.equal = FALSE)$statistic
              aaa <- abs(aaa)
              return(aaa)
     }
     
     Fstat <- function(a, dep){    
              aaa <- anova(lm(a ~factor(dep)))$F[1]
              return(aaa)
     }
     
     u.y <- unique(y)
     n.y <- length(unique(y)) 
      
     if(n.y==2) {     
        cl <- rep(0, length(y)); cl[y==y[1]] <- 1 
        bbb <- apply(x, 2, tstat, dep=cl)
     }

     if(n.y > 2) {     
        cl <- rep(NA, length(y))
        for(i in 1:n.y) cl[y==u.y[i]] <- i 
        bbb <- apply(x, 2, Fstat, dep=cl)
     }
     
           
     q <- quantile(bbb, probs = (1-percent.cut))
     id <- which(bbb >=q)
     return(id) 
}

