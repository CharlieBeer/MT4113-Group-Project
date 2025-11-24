# --------------------
# Helper Functions
# -----------------------------------
# Tolerance Checker
# -----------------------------------
#' @name Norm_Ratio
#' @title Norm Ratio
#' @description
#' Computes the ratio of the Euclidean norm of 2 vectors
#'
#' @param delta step value for each parameter
#' @param theta estimate of the minimised values
#'
#' @returns the ratio of the Euclidean norm of the inputted vectors
#' @author Holly Goldsmith
#' @export

Norm_Ratio <- function(delta, theta){
  return((sqrt(sum(delta^2))) / (sqrt(sum(theta^2))))
}

# -----------------------------------
# Multivariate Newton's Method (MVN)
# -----------------------------------
#' @name MVN
#' @title Multivariate Newton Optimisation
#' @description
#' Performs multivariate optimisation using Newton's method
#'
#' @param f function to be minimised
#' @param inits vector of initial values to be minimised
#' @param data (optional) dataframe with appropriate number of columns for number of variables being optimised
#' @param minimum search for minimum or maximum (takes TRUE if minimum search, FALSE if maximum search)
#' @param tol tolerance level
#' @param maxit maximum number of iterations run before stopping
#' @param method identifier, takes "MVN" only to allow the parent function to call on it
#' @param gradfn (optional) gradient function of f, uses finite differencing otherwise
#' @param hessfn (optional) hessian function of f, uses finite differencing otherwise
#' @param jacobfn NULL variable, included to match parent function inputs
#'
#' @return A list containing:
#' \item{estimate}{Optimised estimate}
#' \item{feval}{Function evaluated at optimised estimate}
#' \item{grad}{Gradient of function at optimised estimate}
#' \item{tolerance}{Tolerance level reached through optimisation}
#' \item{conv}{Whether or not the optimisation converged. 0 - converged, 1 - did not converge, 2 - max iterations reached}
#' \item{niter}{Number of iterations run}
#'
#' @references
#' Swallow, B. (2025, October 29th). Multivariate Optimization (Lecture).
#' University of St Andrews.
#' \url{https://moody.st-andrews.ac.uk/moodle/pluginfile.php/2128841/mod_resource/content/2/Chpater7_12.pdf}
#' @author Holly Goldsmith
#'
#' @note
#' This function is one of 5 to be called with the parent function 'funoptim', also included in this package.
#' @examples
#' f <- function(theta){
#'x <- theta[1]
#'y <- theta[2]
#'return((x-3)^2 + (y+2)^2)
#'}
#'
#'MVN(f = f, inits = c(4, -1), data = NULL, minimum = TRUE, tol = 0.000001, maxit = 100, method = "MVN", gradfn = NULL, hessfn = NULL, jacobfn = NULL)
#' @export

MVN <- function(f, inits, data, minimum, tol, maxit,
                method = "MVN", gradfn, hessfn, jacobfn){

  # Local wrapper function that allows passing of the data argument to grad and hessian, if applicable
  Wrapper <- function(theta, ...){ # ... included when Wrapper was a global function, kept for defensive programming in case of later expansion
    if (is.null(data)){
      f(theta)
    } else {
      f(theta, data)
    }
  }

  # If searching for a maximum, work with the negative of the function, as Newton searches for minimums
  if (minimum == FALSE) {
    og_f <- f
    f <- function(theta, ...){
      return( -og_f(theta, ...))
    }}

  # If not provided with analytic functions, use numDeriv to estimate
  if (is.null(gradfn)){
    gradfn <- function(theta) {numDeriv::grad(Wrapper, theta)}
  }
  if (is.null(hessfn)){
    hessfn <- function(theta) {numDeriv::hessian(Wrapper, theta)}
  }

  theta <- inits # Set staring values
  # Set up some starting values for the loop
  niter <- 0
  loop <- TRUE

  while(loop){
    niter <- niter + 1
    # Compute gradient and hessian to update estimate
    g <- gradfn(theta)
    H <- hessfn(theta)

    delta <- solve(H, -g) # Compute step

    # Check stopping criteria
    if (Norm_Ratio(delta, theta) < tol | niter > maxit) loop <- FALSE
    theta <- theta + delta
  }
  # Check convergence
  # 0: tolerance reached, 1: tolerance not reached, 2: max. iterations reached
  if (Norm_Ratio(delta, theta) < tol) {
    conv <- 0
  } else if (niter > maxit) {
    conv <- 2
  } else {
    conv <- 1
  }

  # Flip back over for final evaluation if maximising
  if (minimum == FALSE) {
    f <- og_f}
  return(list(estimate = theta,
              feval = Wrapper(theta),
              grad = g,
              tolerance = Norm_Ratio(delta, theta),
              conv = conv,
              niter = niter))
}
