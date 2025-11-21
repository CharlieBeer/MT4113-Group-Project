\name{GS}
\alias{GS}

\title{
Grid Search
}
\description{
Grid Search method
}
\usage{
GS(f,inits,...)
}

\arguments{
  \item{f}{
    the function to be optimized
    }
  \item{inits}{
  the initial guess, the function will look at +-10 of this first guess point
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
  return ((x-0.25)*(x-0.6)*(x-5)*(x-3))
}

test_inits<-c(0,10)

GS(test_f,test_inits,minimum=TRUE)
}


