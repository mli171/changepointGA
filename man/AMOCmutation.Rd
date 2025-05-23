% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AMOCfunc.R
\name{AMOCmutation}
\alias{AMOCmutation}
\title{Jump mutation operator to produce offspring for AMOC problem}
\usage{
AMOCmutation(
  child,
  p.range = NULL,
  minDist,
  Pchangepoint = NULL,
  lmax = NULL,
  mmax = NULL,
  N = NULL
)
}
\arguments{
\item{child}{The child chromosome resulting from the \code{crossover} genetic
operator.}

\item{p.range}{The default value is \code{NULL}. If there is no requirement
on model order selection, such an auxiliary argument is needed for \code{GA}
and \code{IslandGA} functions.}

\item{minDist}{The minimum length between two adjacent changepoints in
\code{\link{AMOCselection}} operator, which is also the jump magnitude in the
\code{AMOCmutation} operator.}

\item{Pchangepoint}{An auxiliary argument is needed for \code{GA}
and \code{IslandGA} functions.}

\item{lmax}{An auxiliary argument is needed for \code{GA} and \code{IslandGA}
functions.}

\item{mmax}{An auxiliary argument is needed for \code{GA} and \code{IslandGA}
functions.}

\item{N}{An auxiliary argument is needed for \code{GA} and \code{IslandGA}
functions.}
}
\value{
The resulting child chromosome representation.
}
\description{
In a certain probability, the \code{mutation} genetic operator can be applied
to generate a new \code{child}. In this AMOC mutation operator, the new child
changepoint location can be down via a "jump" method. The child changepoint
location will jump \code{minDist} time units either to the left or right to
produce the changepoint location for the mutated child. The jump direction is
randomly decided with 0.5 probability.
}
