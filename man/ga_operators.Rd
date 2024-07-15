\name{GA_operators}
\alias{GA_operators}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GA_operators
}
\description{
A list object contains the GA operators function names.
}
\arguments{
\item{population}{It could be a function or a matrix. The function should be designed for initializing a population. The default population initialization is random initialization with some imposed constraints. See \code{\link{random_population_cpp}} for example. The function returned object is a matrix, \code{pop}. The users can specified their own \code{population} function. It could also be a matrix object, which contain the user specified chromosome. By default, each column represents one individual chromosome.}
\item{selection}{A function can help select \code{mom} and \code{dad} from current generation population, where \code{dad} is set to have better fit (smaller fitness function values). The default for selection uses the linear rank selection method. See \code{\link{selection_linearrank_cpp}} for example. The function returned object is a list contain the chromosomes for \code{mom} and \code{dad}.}
\item{crossover}{A function can apply crossover to the chosen parents to produce child for next generation with specified probability. The default for crossover uses the uniform crossover method. See \code{\link{offspring_uniformcrossover_cpp}} for example. The function returned object is a vector contain the chromosomes for \code{child}.}
\item{mutation}{A function can apply mutation to the produced child with the specified probability \code{Pmutation} in \code{\link{GA_param}} and \code{\link{IslandGA_param}}. If order selection is not requested (\code{option = "cp"} in \code{GA_param} and \code{Island_GA}), the default mutation operator function uses \code{selectTau_cpp} to select a completely new individual with a new chromosome as the mutated child. For details, see \code{\link{selectTau_cpp}}. If order selection is needed (\code{option = "both"} in \code{GA_param} and \code{Island_GA}), we first decide whether to keep the produced child's model order with a probability of 0.5. If the child's model order is retained, the \code{selectTau_cpp} function is used to select a completely new individual with a new chromosome as the mutated child. If a new model order is selected from the candidate model order set, there is a 0.5 probability to either select a completely new individual with new changepoint locations or retain the original child's changepoint locations for the mutated child. Note that the current model orders in the child's chromosome are excluded from the set to avoid redundant objective function evaluation. Finally, the function returns a vector containing the modified chromosomes for mutated \code{child}.}
}
\details{
\itemize{Population Initialization. By default, the specified population function should include the following arguments as \code{population(popsize, N, minDist, Pchangepoint, mmax, lmax)}, where
\item{popsize} The size of population in each generation.
\item{N} The sample size of the time series.
\item{minDist} The required minimum distance between two adjacent changepoints.
\item{Pchangepoint} The probability of changepoints for every time series.
\item{mmax} The user specified maximum number of changepoints, by default, as \code{N/2 - 1}.
\item{lmax} The user specified maximum length of individual chromosome, by default, as \code{2+N/2 + 1}.
}
\itemize{Parent Selection. The genetic algorithm require to select a pair of chromosomes, representing dad and mom, for the \code{crossover} operator. The parents chromosomes are randomly selectd from initialized population by a linear ranking method. By default, the dad has better fit/smaller fitness function value/larger rank than mom. The usage of \code{selection(pop, popFit)} needs two arguments,
\item{pop} A matrix contains all individuals chromosomes from current population.
\item{popFit} A vector contains the fitness function values for the corresponding individuals chromosomes.
}
\itemize{Crossover. In a certial probability, the \code{crossover} genetic operator can be applied to the selected parents to produce the offspring for next generation. The usage of \code{crossover(mom, dad, minDist, lmax, N)} needs the following arguments,
\item{mom} A vector contains the selected chromosome representation with lower fitness function value.
\item{dad} A vector contains the selected chromosome representation with higher fitness function value.
\item{minDist} The required minimum distance between two adjacent changepoints.
\item{lmax} The user specified maximum length of individual chromosome, by default, as \code{2+N/2 + 1}.
\item{N} The sample size of the time series.
}
\itemize{mutation. In a certial probability, the \code{mutation} genetic operator can be applied to generate a new \code{child}. By default, the new child selection can be down by the similar individual selection method in population initialization. The usage of \code{mutation(minDist, Pchangepoint, lmax, mmax, N)} needs the following arguments,
\item{minDist} The required minimum distance between two adjacent changepoints.
\item{Pchangepoint} The probability of changepoints for every time series.
\item{mmax} The user specified maximum number of changepoints, by default, as \code{N/2 - 1}.
\item{lmax} The user specified maximum length of individual chromosome, by default, as \code{2+N/2 + 1}.
\item{N} The sample size of the time series.
}
}
\seealso{See Also as \code{\link{GA}}.}
\references{Lu, Q., Lund, R., & Lee, T. C. (2010). An MDL approach to the climate segmentation problem.
\emph{Ann. Appl. Stat.} 4 (1) 299 - 319.}
\examples{
GA_operators = list(population = "random_population_cpp",
                    selection  = "selection_linearrank",
                    crossover  = "offspring",
                    mutation   = "mutation")
}
