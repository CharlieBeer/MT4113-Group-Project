# -----------------------------------
# Gauss-Newton Method (GN)
# -----------------------------------
#' @name GN
#' @title Gauss-Newton Method
#' @description
#' Performs non-linear least squares optimisation using Gauss-Newton method
#'
#' @param f the function to be optimised
#' @param inits a vector of initial guesses for each parameter to be optimised.
#' @param data dataframe containing the data
#' @param minimum not relevant, included to match parent function inputs
#' @param tol tolerance level
#' @param maxit maximum number of iterations run before stopping
#' @param method identifier, takes "GN" only to allow the parent function to call on it
#' @param gradfn NULL variable, included to match parent function inputs
#' @param hessfn NULL variable, included to match parent function inputs
#' @param jacobfn (optional) Jacobian function of residuals, uses numerical approximation otherwise
#'
#' @return A list containing:
#' \item{estimate}{Optimised estimate}
#' \item{feval}{Function evaluated at optimised estimate}
#' \item{grad}{Gradient of function at optimised estimate}
#' \item{tolerance}{Tolerance level reached through optimisation}
#' \item{conv}{Whether or not the optimisation converged. 0 - converged, 1 - failed, 2 - max iterations reached}
#' \item{niter}{Number of iterations ran}
#' @references
#' Swallow, B. (2025). Non-Linear Least Squares.
#' University of St Andrews.
#' \url{https://moody.st-andrews.ac.uk/moodle/pluginfile.php/2128840/mod_resource/content/4/_book/non-linear-least-squares.html#gauss-newton-method}
#' @author Ilana Goldman
#' @examples
#' f <- function(theta, data){
#'   a <- theta[1]
#'   b <- theta[2]
#'   return(a * exp(b * data$x))
#' }
#'
#' data <- list(x = 1:10, y = 2.5 * exp(-1.5 * 1:10))
#' inits <- c(1, -1)
#' 
#' GN(f, inits, data)
#' @export

GN <- function(f, inits, data=NULL,minimum=TRUE, tol = 1e-10, maxit = 1000, method = "GN", gradfn = NULL, hessfn = NULL, jacobfn = NULL) {

  #checks: stop if there is no dataframe
  if (is.null(data)) {
    stop("Error: a dataframe is required for the Gauss-Newton method")
  }

  #residuals function
  resids <- function(theta, data) {
    pred <- f(theta, data)
    r <- data$y - pred
    return(r)
  }

  #compute jacobian if not provided
  if (is.null(jacobfn)) {
    jacobian_function <- function(theta, data) {
      numDeriv::jacobian(resids, theta, data = data)
    }
  } else {
    jacobian_function <- jacobfn
  }

  #compute Euclidean norm
  norm <- function(x) {
    return(sqrt(sum(x^2)))
  }

  #parameters for the loop
  theta_current <- inits
  iter <- 0
  conv <- 2  #default to max iterations reached - this can be updated if something else happens
  final_tolerance <- NA

  #begin loop
  for (iter in 1:maxit) {

    #compute Jacobian
    J <- jacobian_function(theta_current, data)
    #compute residuals
    residuals <- resids(theta_current, data)

    #compute gradient
    gradient <- t(J) %*% residuals

    #compute Hessian
    hessian <- t(J) %*% J

    #solve for the parameter update
    delta <- solve(hessian, -gradient)

    #update theta
    theta_new <- theta_current + delta

    #check the convergence - is it below the tolerance yet?
    increment <- norm(delta) / norm(theta_current)
    final_tolerance <- increment  #keep track of tolerance

    if (increment < tol) { #convergence has been reached
      conv <- 0  #update convergence code
      break
    }

    #update theta for the next iteration
    theta_current <- theta_new
  }

    #compute the final residuals and final value of the function
    final_residuals <- resids(theta_current,data)
    feval <- sum(final_residuals^2)

    #results
    result <- list(estimate = theta_current, feval = feval, grad = gradient,
                   tolerance = final_tolerance, conv = conv, niter = iter)

  return(result)
}
