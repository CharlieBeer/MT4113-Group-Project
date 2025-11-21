\name{MVN}
\alias{MVN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multivariate Newton Optimisation
}
\description{
Performs a multivariate optimisation using Newton's method.
}
\usage{
MVN(f, inits, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{A function to be optimised}
  \item{inits}{The initial guess; a numeric vector of length the number of parameters f takes}
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
\examples{
f <- function(theta){
  x <- theta[1]
 y <- theta[2]
  return((x-3)^2 + (y+2)^2)
}

MVN(f = f, inits = c(4, -1), data = NULL, minimum = TRUE, tol = 0.000001, maxit = 100, method = "MVN", gradfn = NULL, hessfn = NULL, jacobfn = NULL)
}
\seealso{
  \code({\link(funoptim))}
}
