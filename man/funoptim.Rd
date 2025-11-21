\name{funoptim}
\alias{funoptim}
\title{
Unified Optimiser
}
\description{
Perfomrs a variety of optimisation methods based on the 'method' input.
}
\usage{
funoptim(f, inits=NULL, data=NULL, minimum=TRUE, tol=1e-8, maxit=100,
                     method=NULL, gradfn=NULL, hessfn=NULL, jacobfn=NULL))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{f}{A function to be optimised}
  \item{inits}{The initial guess; a numeric vector of length the number of parameters f takes}
  \item{data}{[optional] A relevant dataframe relating to f, with the same number of columns according to the number of parameters to be optimised}
  \item{minimum}{A logical value. TRUE or FALSE depending on whether f is to be minimised or maxiised respectively}
  \item{tol}{A numeric scalar. The difference from the true value at which the optimisation will be considered complete}
  \item{maxit}{A numeric scalar. The maximum numer of iterations before the optimisation terminates}
  \item{method}{A string. Takes values "GS", "BS", "UVN", "MVN", "GN" only, used to specify the method to be called upon}
  \item{gradfn}{[optional] A gradient function of f, uses finite differencing otherwise}
  \item{hessfn}{[optional] A hessian function of f, uses finite differencing otherwise}
  \item{jacobfn}{[optional] A jecobian function of f, uses finite differencing otherwise}
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
Swallow, B. (2025). Computer-intensive Statistics. [https://moody.st-andrews.ac.uk/moodle/pluginfile.php/2128841/mod_resource/content/2/Chpater7_12.pdf]
}
\examples{
# Multivariate Newton (MVN)

f <- function(theta){
  x <- theta[1]
  y <- theta[2]
  return((x-3)^2 + (y+2)^2)
}

funoptim(f = f, inits = c(4, -1))
}
