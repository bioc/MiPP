\name{mipp.preproc}
\alias{mipp.preproc}
\title{
Preprocessing
}
\description{
Performs IQR normalization, thesholding, and log2-transformation 
}
\usage{
function (x, data.type = "MAS5")
}

\arguments{
\item{x}{data}
\item{data.type}{data type is MAS5, MAS4, or dChip}
}

\seealso{
\code{\link{mipp}}
}



\examples{

library(MiPP)

data(colon)
colon.nor <- mipp.preproc(colon)

}

\keyword{models}
