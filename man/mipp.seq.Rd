\name{mipp.seq}
\alias{mipp.seq}
\title{MiPP-based Classification}
\description{
sequentially finds optimal sets of genes for classification
}
\usage{
mipp.seq(x, y, x.test = NULL, y.test = NULL, probe.ID = NULL, 
    rule = "lda", method.cut = "t.test", percent.cut = 0.01, 
    model.sMiPP.margin = 0.01, min.sMiPP = 0.85, n.drops = 2, 
    n.fold = 5, p.test = 1/3, n.split = 20, 
    n.split.eval = 100, n.seq=3, cutoff.sMiPP=0.7, remove.gene.each.model="all") 
}

\arguments{

\item{x}{data matrix}
\item{y}{class vector}
\item{x.test}{test data matrix if available}
\item{y.test}{test class vector if available}
\item{probe.ID}{probe set IDs; if NULL, row numbers are assigned.}
\item{rule}{classification rule: "lda","qda","logistic","svmlin","svmrbf"; 
the default is "lda".}
\item{method.cut}{method for pre-selection; t-test is available.}
\item{percent.cut}{proportion of pre-selected genes; the default is 0.01.}
\item{model.sMiPP.margin}{smallest set of genes s.t. sMiPP <= (max sMiPP-model.sMiPP.margin); the default is 0.01.}
\item{min.sMiPP}{Adding genes stops if max sMiPP is at least min.sMiPP; 
the default is 0.85.} 
\item{n.drops}{Adding genes stops if sMiPP decreases (n.drops) times, in 
addition to min.sMiPP criterion.; the default is 2.} 
\item{n.fold}{number of folds; default is 5.}
\item{p.test}{partition percent of train and test samples when test samples are not available; the default is 1/3 for test set.}
\item{n.split}{number of splits; the default is 20.}
\item{n.split.eval}{numbr of splits for evalutation; the default is 100.}
\item{n.seq}{Number of sequential gene selection}
\item{cutoff.sMiPP}{cutoff point to select gene models based on 5 percent sMiPP}
\item{remove.gene.each.model}{Re-run after removing all genes in the selected models if "all" 
      and the first gene for each of the selected models if "first"}
}

\value{
\item{model}{candiadate genes (for each split if no indep set is available}
\item{model.eval}{Optimal sets of genes for each split when no indep set is available}
}

\references{
Soukup M, Cho H, and Lee JK (2005). Robust classification modeling on microarray data 
using misclassification penalized posterior, Bioinformatics, 21 (Suppl): i423-i430.

Soukup M and Lee JK (2004). Developing optimal prediction models for cancer classification 
using gene expression data, Journal of Bioinformatics and Computational Biology, 1(4) 681-694
}

\author{
 Soukup M, Cho H, and Lee JK 
}


\examples{

##########
#Example 1: When an independent test set is available

data(leukemia)

#Normalize combined data
leukemia <- cbind(leuk1, leuk2)
leukemia <- mipp.preproc(leukemia, data.type="MAS4")

#Train set
x.train <- leukemia[,1:38]
y.train <- factor(c(rep("ALL",27),rep("AML",11)))

#Test set
x.test <- leukemia[,39:72]
y.test <- factor(c(rep("ALL",20),rep("AML",14)))


#Compute MiPP
out <- mipp.seq(x=x.train, y=y.train, x.test=x.test, y.test=y.test, n.fold=5, percent.cut=0.05, rule="lda", n.seq=3)

#Print candidate models
out$model

#Print the genes selected
out$genes.selected


##########
#Example 2: When an independent test set is not available

data(colon)

#Normalize data
x <- mipp.preproc(colon)
y <- factor(c("T", "N", "T", "N", "T", "N", "T", "N", "T", "N", 
       "T", "N", "T", "N", "T", "N", "T", "N", "T", "N",
       "T", "N", "T", "N", "T", "T", "T", "T", "T", "T", 
       "T", "T", "T", "T", "T", "T", "T", "T", "N", "T", 
       "T", "N", "N", "T", "T", "T", "T", "N", "T", "N", 
       "N", "T", "T", "N", "N", "T", "T", "T", "T", "N", 
       "T", "N"))


#Deleting comtaminated chips
x <- x[,-c(51,55,45,49,56)]
y <- y[ -c(51,55,45,49,56)]

#Compute MiPP
out <- mipp.seq(x=x, y=y, n.fold=5, p.test=1/3, n.split=5, n.split.eval=100, 
percent.cut= 0.1, rule="lda", n.seq=3)


#Print candidate models for each split
out$model

#Print optimal models and independent evaluation for each split
out$model.eval

#Print the genes selected
out$genes.selected

}

\keyword{models}


