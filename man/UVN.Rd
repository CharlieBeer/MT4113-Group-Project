\name{UVN}
\alias{UVN}

\title{
Univariate Newton Optimisation
}
\description{
Performs a univariate optimisation using Newton's method.
}
\usage{
UVN(f, inits, ...)
}

\arguments{
  \item{f}{The function to be optimised}
  \item{inits}{The initial value; a numeric scalar giving the starting point}
  \item{...}{Any other arguments passed}
}
\value{
  \item{estimate}{Optimised estimate}
  \item{feval}{Function evaluated at optimised estimate}
  \item{grad}{Gradient of function at optimised estimate}
  \item{tolerance}{Tolerance level reached through optimisation}
  \item{conv}{Whether or not the optimisation converged. 0 - converged, 1 - did not converge, 2 - max iterations reached}
  \item{niter}{Number of iterations ran}
}
\references{
Swallow, B. (2025, October 29th). Multivariate Optimisation [Lecture]. University of St Andrews. [https://moody.st-andrews.ac.uk/moodle/pluginfile.php/2128841/mod_resource/content/2/Chpater7_12.pdf]
}
\author{
Chenyu Lin
}

\examples{
f <- function(theta){
  x <- theta[1]
  return((x - 33)^3)
}

UVN(f = f, inits = 10, data = NULL, minimum = TRUE,
    tol = 1e-6, maxit = 100, method = "UVN",
    gradfn = NULL, hessfn = NULL, jacobfn = NULL)
}
