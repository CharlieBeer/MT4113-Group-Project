\name{GN}
\alias{GN}

\title{Gauss-Newton}
\description{
Performs Gauss-Newton optimisation on a function.
}
\usage{
GS(f, inits, ...)
}

\arguments{
  \item{f}{
    the function to be optimised - if no dataframe provided, f MUST return vector of residuals directly
    }
  \item{inits}{
  a vector of initial guesses for each parameter to be optimised.
  }
  \item{...}{
  any other arguments passed
  }
}

\value{
  \item{estimate}{
  optimised estimate
  }
  \item{feval}{
  function evaluated at optimised estimate
  }
  \item{grad}{
  gradient of function at optimised estimate
  }
  \item{tolerance}{
  tolerance level reached through optimisation
  }
  \item{conv}{
  whether or not the optimisation converged. 0 - converged, 1 - did not converge, 2 - max iterations reached
  }
  \item{niter}{
  number of iterations ran
  }

\references{
Swallow, B. (2025, October 29th). Multivariate Optimization [Lecture]. University of St Andrews. [https://moody.st-andrews.ac.uk/moodle/pluginfile.php/2128841/mod_resource/content/2/Chpater7_12.pdf]
}
\author{
Ilana Goldman
}

\examples{
test_f <- function(theta) {
  x <- theta[1]
  y <- theta[2]
  return(c(x - 3, y + 2))
}

test_inits <- c(0, 0)

GN(f = test_f, inits = test_inits, data = NULL)
