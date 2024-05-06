\name{random_population_cpp}
\alias{random_population_cpp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Capitalized Title (ideally limited to 65 characters)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
random_population_cpp(popsize, n, minDist, Pb, mmax, lmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{popsize}{
%%     ~~Describe \code{popsize} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{minDist}{
%%     ~~Describe \code{minDist} here~~
}
  \item{Pb}{
%%     ~~Describe \code{Pb} here~~
}
  \item{mmax}{
%%     ~~Describe \code{mmax} here~~
}
  \item{lmax}{
%%     ~~Describe \code{lmax} here~~
}
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
##--	or standard data sets, see data().

## The function is currently defined as
function (popsize, n, minDist, Pb, mmax, lmax) 
{
    .Call("_changepointGA_random_population_cpp", PACKAGE = "changepointGA", 
        popsize, n, minDist, Pb, mmax, lmax)
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
