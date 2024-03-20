\name{ga_operators}
\alias{ga_operators}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ga_operators
}
\description{
A list object contains the GA operators function names.
}
\usage{
ga_operators = list(population = "random_population_cpp",
                    selection  = "selection_linearrank",
                    crossover  = "offspring",
                    mutation   = "mutation")
}
\arguments{
  \item{population}{A function can initialize a population. The default population initialization is random initialization with some imposed constraints. See \code{\link{random_population_cpp}} for more details The users can specified their own \code{population} function.}
  \item{selection}{A function can help select \code{mom} and \code{dad} from current generation population, where \code{dad} is set to have better fit (smaller fitness function values). The default for selection uses the linear rank selection method. See \code{\link{selection_linearrank}} for more details.}
  \item{crossover}{A function can apply crossover to the chosen parents to produce child for next generation with specified probability. The default for crossover uses the uniform crossover method. See \code{\link{offspring_uniformcrossover} for more details.}
  \item{mutation}{A function can apply mutation to the produced child with specified probability. The default for mutation uses the select a complete new individual with new chromosome as the mutated child. See \code{\link{selectTau_cpp}} for more details.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
