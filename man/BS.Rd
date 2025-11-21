\name{BS}
\alias{BS}

\title{
Bisection
}
\description{
Performs the Bisection method of optimization
}
\usage{
BS(f,inits,...)
}

\arguments{
  \item{f}{
    the function to be optimized
    }
  \item{inits}{
  the initial guess range, this should be two values of at which the function has different signed gradients.
  }
  \item{...}{
  any other arguments passed
  }
}

\value{
  \item{estimate}{
  Optimised estimate
  }
  \item{feval}{
  Function evaluated at optimised estimate
  }
  \item{grad}{
  Gradient of function at optimised estimate
  }
  \item{tolerance}{
  Tolerance level reached through optimisation
  }
  \item{conv}{
  Whether or not the optimisation converged. 0 - converged, 1 - did not converge, 2 - max iterations reached
  }
  \item{niter}{
  Number of iterations ran
  }

}
\references{
Swallow, B. (2025, October 29th). Multivariate Optimisation [Lecture]. University of St Andrews. [https://moody.st-andrews.ac.uk/moodle/pluginfile.php/2128841/mod_resource/content/2/Chpater7_12.pdf]
}
\author{
Charles Beer
}

\examples{
test_f<-function(x){
  return (x^5-3*x^4+12*x^3-5*x^2-3)
}

test_inits<-c(0.1,5)

BS(test_f,test_inits)
}
