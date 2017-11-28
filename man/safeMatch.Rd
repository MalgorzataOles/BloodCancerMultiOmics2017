\name{safeMatch}
\alias{safeMatch}

\title{safe version of the match function that throws an error if there is no match}

\usage{ safeMatch (x, ...) }
\arguments{ 

\item{x}{ string to be matched, will be passed on as first argument to match  } 

\item{...}{ passed on to match }
}

\description{ While \code{\link{match}} returns an NA if no match is found, this function will throw an error }
 
\examples{ safeMatch("oranges",  c("apples", "oranges") ) }
