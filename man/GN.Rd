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
    the function to be optimised
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
Swallow, B. (2025, October 29th). Multivariate Optimisation [Lecture]. University of St Andrews. [https://moody.st-andrews.ac.uk/moodle/pluginfile.php/2128841/mod_resource/content/2/Chpater7_12.pdf]
}
\author{
Ilana Goldman
}

\examples{
test_f<-function(x){
  return (x^4-6*x^3+x^2-2*x^2-9)
}

test_inits <- c(1.0, 1.5, 0.2)

GN(test_f,test_inits)
}
