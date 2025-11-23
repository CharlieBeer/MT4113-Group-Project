# -----------------------------------
# Parent Optimisation Function
# -----------------------------------
#' @name funoptim
#' @title Parent Optimisation Function
#' @description
#' Performs a variety of optimisation methods based on the 'method' input.
#' @param f Function to be optimised
#' @param inits Initial parameter values
#' @param data (optional) Dataframe.
#' @param minimum (optional) Search for minimum or maximum (takes TRUE if minimum search, FALSE if maximum search)
#' @param tol (optional) Tolerance level - defaults to 1e-8
#' @param maxit (optional) Maximum number of iterations run before stopping - defaults to 100
#' @param method optimisation method, takes one of "GS", "BS", "UVN", "MVN", "GN"
#' @param gradfn (optional) A gradient function
#' @param hessfn (optional) A hessian function
#' @param jacobfn (optional) A jacobian function

#' @return A list containing:
#' \item{estimate}{Optimised estimate}
#' \item{feval}{Function evaluated at optimised estimate}
#' \item{grad}{Gradient of function at optimised estimate}
#' \item{tolerance}{Tolerance level reached through optimisation}
#' \item{conv}{Whether or not the optimisation converged. 0 - converged, 1 - did not converge, 2 - max iterations reached}
#' \item{niter}{Number of iterations run}
#'
#' @references
#' Swallow, B. (2025). Univariate Optimization.
#' University of St Andrews.
#' \url{https://moody.st-andrews.ac.uk/moodle/pluginfile.php/2128841/mod_resource/content/2/Chpater7_12.pdf}
#' @author Charlie Beer, Holly Goldsmith, Ilana Goldman, Chenyu Lin
#' @examples
#' # Multivariate Newton (MVN)
#'f <- function(theta){
#'  x <- theta[1]
#'  y <- theta[2]
#'  return((x-3)^2 + (y+2)^2)
#'}
#'funoptim(f = f, inits = c(4, -1),method="MVN")
#' @export

funoptim <- function(f, inits=NULL, data=NULL, minimum=TRUE, tol=1e-8, maxit=100,
                     method=NULL, gradfn=NULL, hessfn=NULL, jacobfn=NULL) {
  # Error Checking:
  if(!is.function(f))
    stop("Error: 'f' must be a function.")
  if(!is.numeric(inits))
    stop("Error: 'inits' must be a numeric vector.")
  if(!identical(minimum, TRUE) && !identical(minimum, FALSE))
    stop("Error: 'minimum' must be either TRUE or FALSE")
  if(!method %in% c("GS", "BS", "UVN", "MVN", "GN"))
    stop("Error: method provided is invalid.")

  if (method=="GS"){
    return(GS(f, inits, data, minimum, tol, maxit,
              method, gradfn, hessfn, jacobfn))
  }
  if (method=="BS"){
    return(BS(f, inits, data, minimum, tol, maxit,
              method, gradfn, hessfn, jacobfn))
  }
  if (method=="UVN"){
    return(UVN(f, inits, data, minimum, tol, maxit,
              method, gradfn, hessfn, jacobfn))
  }
  if (method=="MVN"){
    return(MVN(f, inits, data, minimum, tol, maxit,
              method, gradfn, hessfn, jacobfn))
  }
  if (method=="GN"){
    return(GN(f, inits, data, minimum, tol, maxit,
              method, gradfn, hessfn, jacobfn))
  }

  stop("method not named correctly")

}



