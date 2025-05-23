% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{selection_linearrank}
\alias{selection_linearrank}
\title{The default parents selection genetic algorithm operator}
\usage{
selection_linearrank(pop, popFit)
}
\arguments{
\item{pop}{A matrix contains the chromosomes for all individuals. The number of
rows is equal to \code{lmax} and the number of columns is equal to the
\code{popsize}.}

\item{popFit}{A vector contains the objective function value (population fit)
being associated to each individual chromosome from above.}
}
\value{
A list contains the chromosomes for \code{dad} and \code{mom}.
}
\description{
The genetic algorithm require to select a pair of chromosomes, representing
\code{dad} and \code{mom}, for the \code{crossover} operator to
produce offspring (individual for next generation). The parents chromosomes
are randomly selectd from the initialized population by a linear ranking
method according to each individual's fittness in the input argument
\code{popFit}. By default, the dad has better fit/smaller fitness function
value/larger rank than \code{mom}.
}
