\name{offspring_uniformcrossover}
\alias{offspring_uniformcrossover}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
offspring_uniformcrossover
}
\description{
In uniform crossover, typically, each bit is chosen from either parent with equal probability. Other mixing ratios are sometimes used, resulting in offspring which inherit more genetic information from one parent than the other. In a uniform crossover, we don’t divide the chromosome into segments, rather we treat each gene separately. In this, we essentially flip a coin for each chromosome to decide whether or not it will be included in the off-spring.
}
\usage{
offspring_uniformcrossover(mom, dad, minDist, lmax, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mom}{Among two selected individuals, \code{mom} represents the selected chromosome representation with lower fitness function value.}
  \item{dad}{Among two selected individuals, \code{dad} represents the selected chromosome representation with larger fitness function value.}
  \item{minDist}{The minimum length between two adjacent changepoints.}
  \item{lmax}{The maximum possible length of the chromosome representation.}
  \item{n}{The length/sample size of the time series.}
}
\details{TBD
%%  ~~ If necessary, more details than the description above ~~
}
\references{TBD
%% ~put references to the literature/web site here ~
}
\author{Mo Li}
