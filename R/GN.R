# -----------------------------------
# Gauss-Newton Method (GN)
# -----------------------------------
#' @description
#' Performs non-linear least squares optimisation using Gauss-Newton method
#'
#' @param f function that returns predicted values given parameters and data
#' @param inits vector of initial values for parameters
#' @param data (optional, if f returns residuals) dataframe containing the data
#' @param tol tolerance level
#' @param maxit maximum number of iterations run before stopping
#' @param method identifier, takes "GN" only to allow the parent function to call on it
#' @param gradfn NULL variable, included to match parent function inputs
#' @param hessfn NULL variable, included to match parent function inputs
#' @param jacobfn (optional) Jacobian function of residuals, uses numerical approximation otherwise
#'
#' @returns a list containing the optimised estimate, the function evaluated at said estimate, the gradient of the function at the time, the tolerance level, whether the optimisation converged and the number of iterations ran.
#' @export

GN <- function(f, inits, data=NULL, tol = 1e-10, maxit = 1000, method = "GN", gradfn = NULL, hessfn = NULL, jacobfn = NULL) {

  #internal residuals function - modified to take data=NULL
  resids <- function(theta, data) {
    if (is.null(data)) {
      #when data=NULL return residuals directly
      return(f(theta))
    } else {
      #when data is provided, compute residuals as usual
      pred <- f(theta, data)
      r <- data$y - pred
      return(r)
    }
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
