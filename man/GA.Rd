\name{GA}
\alias{changepointGA}
\title{GA}
\usage{
GA(fitness, n, GA_param, ga_operators, ...)
}
\description{
It runs the basic genetic algorithm (GA). There are 4 operator functions that need to be specified, \code{population}, \code{selection}, \code{crossover}, \code{mutation}.

Steady State
In steady state GA, we generate one or two off-springs in each iteration and they replace one or two individuals from the population. A steady state GA is also known as Incremental GA.



}
\arguments{
  \item{fitness}{A fitness function to be minimized, with certain changepoint configuration and chromosome representation.}
  \item{n}{The number of timepoints (sample size) of the data.}
  \item{GA_param}{A list contains the GA hyperparameters. See \code{\link{GA_param}} for details.}
  \item{ga_operators}{A list contains the GA operators, population, selection, crossover, mutation.See \code{\link{ga_operators}} for details.}
  \item{...}{Further arguments to be passed to \code{fitness} function.}
}
\details{\acronym{AOP} stands for
  \emph{autoregressive ordered probit} model, a dynamic generalizations of the
  ordered probit regression model or cumulative probit
  model, see Muller and Czado (2005) for details. Here, \acronym{AOP} models are fitted by the maximum pairwise
  likelihood method of order p (Varin and Vidoni, 2006, 2009).}

\value{
  \item{overbestfit}{The value of the fitness function value of the overall best-fitted chromosome for all generation.}
  \item{overbestchrom}{A vector contains the overall best-fitted chromosome for all generation.}
  \item{bestfit}{A vector that contains the fitness function value of the best chromosome for each generation.}
  \item{bestchrom}{A matrix that contains the best-fitted chromosome for each generation.}
  \item{count}{Number of generations that the GA runs.}
  \item{convg}{An integer code.}
  \item{}{0: The GA converge before the prespecified maximum number of iterations.}
  \item{}{1: Maximum number iterations reached before convergence.}
}

\examples{
}
