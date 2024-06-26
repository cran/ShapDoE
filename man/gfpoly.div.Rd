% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estsh.R
\name{gfpoly.div}
\alias{gfpoly.div}
\title{Polynomial division defined on GF(s) with a prime s}
\usage{
gfpoly.div(f1, f2, s)
}
\arguments{
\item{f1}{a vector represents the coefficients of the dividend polynomial. For example, if the dividend is x^5+2x^3+1, then f1=c(1,0,2,0,0,1).}

\item{f2}{a vector represents the coefficients of the dividend polynomial. For example, if the divisor is x^4+2, then f2=c(1,0,0,0,2).}

\item{s}{a prime, the order of the Galois filed (GF).}
}
\value{
a vector represents the coefficients of the resulting polynomial. For example, the result c(2,0,1,1) represents 2x^3+x+1.
}
\description{
Polynomial division defined on GF(s) with a prime s
}
\examples{
gfpoly.div(c(1,0,2,0,0,1),c(1,0,0,0,2),3)
}
