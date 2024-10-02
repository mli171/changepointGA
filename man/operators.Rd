\name{operators}
\alias{operators}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
operators
}
\description{
A list object contains the GA or IslandGA operators function names.
}
\arguments{
\item{population}{It could be a function or a matrix. The function should be designed for initializing a population. The default population initialization is random initialization with some imposed constraints. See \code{\link{random_population}} for example. The function returned object is a matrix, \code{pop}. The users can specified their own \code{population} function. It could also be a matrix object, which contain the user specified chromosome. By default, each column represents one individual chromosome.}
\item{selection}{A function can help select \code{mom} and \code{dad} from current generation population, where \code{dad} is set to have better fit (smaller fitness function values). The default for selection uses the linear rank selection method. See \code{\link{selection_linearrank}} for example. The function returned object is a list contain the chromosomes for \code{mom} and \code{dad}.}
\item{crossover}{A function can apply crossover to the chosen parents to produce child for next generation with specified probability. The default for crossover uses the uniform crossover method. See \code{\link{uniformcrossover}} for details in the default crossover operator. The function returned object is a vector contain the chromosomes for \code{child}.}
\item{mutation}{A function can apply mutation to the produced child with the specified probability \code{Pmutation} in \code{\link{GA_param}} and \code{\link{IslandGA_param}}. See \code{\link{mutation}} for details in the default mutation operator.}
}
\seealso{See Also as \code{\link{GA}}.}
\references{Lu, Q., Lund, R., & Lee, T. C. (2010). An MDL approach to the climate segmentation problem.
\emph{Ann. Appl. Stat.} 4 (1) 299 - 319.}
\examples{
operators = list(population = "random_population",
                 selection  = "selection_linearrank",
                 crossover  = "uniformcrossover",
                 mutation   = "mutation")
}
