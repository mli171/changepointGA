% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cptgaisl-class.R
\name{print.summary.cptgaisl}
\alias{print.summary.cptgaisl}
\title{Print Summary for a `cptgaisl` Object}
\usage{
\method{print}{summary.cptgaisl}(x, digits = getOption("digits"), max_display = 5, ...)
}
\arguments{
\item{x}{An object of class \code{cptgaisl}, typically produced by a GA-based changepoint detection routine.}

\item{digits}{Number of digits to print for probabilities and fitness. Default taken from \code{getOption("digits")}.}

\item{max_display}{Maximum number of suggested solutions to display if suggestions are provided.}

\item{...}{Additional arguments (currently not used).}
}
\value{
Invisibly returns \code{NULL}. Called for its side effect of printing to the console.
}
\description{
Displays key information about the settings and results from a changepoint detection
procedure using the Island model Genetic Algorithm (GA) stored in a `cptgaisl` object. This includes
the algorithm configuration, population settings, optimization mode, and final
solution such as the number and location of changepoints and model parameters (if applicable).
}
\details{
When the GA is run in \code{option = "cp"} mode, only changepoint locations are shown.
If \code{option = "both"}, the output includes the selected model hyperparameters along
with changepoint locations.

The function uses plain text output to print a formatted summary to the console. If
\code{x@suggestions} is provided, only up to \code{max_display} suggestions will be shown.
}
\seealso{
\code{\link[=summary,cptgaisl-method]{summary}}, \code{\link{cptgaisl-class}}, \code{\link{plot.cptgaisl}}
}
