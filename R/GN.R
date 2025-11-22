# -----------------------------------
# Gauss-Newton Method (GN)
# -----------------------------------
#' @description
#' Performs non-linear least squares optimisation using Gauss-Newton method
#'
#' @param f function that returns predicted values given parameters and data
#' @param inits vector of initial values for parameters
#' @param data dataframe containing the data
#' @param tol tolerance level
#' @param maxit maximum number of iterations run before stopping
#' @param method identifier, takes "GN" only to allow the parent function to call on it
#' @param gradfn NULL variable, included to match parent function inputs
#' @param hessfn NULL variable, included to match parent function inputs
#' @param jacobfn (optional) Jacobian function of residuals, uses numerical approximation otherwise
#'
#' @returns a list containing the optimised estimate, the function evaluated at said estimate, the gradient of the function at the time, the tolerance level, whether the optimisation converged and the number of iterations ran.
#' @export

GN <- function(f, inits, data=NULL, tol = 1e-10, maxit = 1000,
                       method = "GN", gradfn = NULL, hessfn = NULL, jacobfn = NULL) {

  #create a wrapper function to handle a case with data=NULL
  if (is.null(data)) {
    f_wrapper <- function(theta) f(theta)
    resids_wrapper <- function(theta) {
      preds <- f_wrapper(theta)
      return(preds)
    }
  } else {
    f_wrapper <- function(theta) f(theta, data)
    resids_wrapper <- function(theta) {
      preds <- f_wrapper(theta)
      r <- data$y - preds
      return(r)
    }
  }

  #compute jacobian if not provided
  if (is.null(jacobfn)) {
    jacobian_function <- function(theta) {
      numDeriv::jacobian(resids_wrapper, theta)
    }
  } else {
    if (is.null(data)) {
      jacobian_function <- function(theta) jacobfn(theta)
    } else {
      jacobian_function <- function(theta) jacobfn(theta, data)
    }
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
    J <- jacobian_function(theta_current)
    #compute residuals
    residuals <- resids_wrapper(theta_current)

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
    final_residuals <- resids_wrapper(theta_current)
    feval <- sum(final_residuals^2)

    #results
    result <- list(estimate = theta_current, feval = feval, grad = gradient, 
                   tolerance = final_tolerance, conv = conv, niter = iter)

  return(result)
}



