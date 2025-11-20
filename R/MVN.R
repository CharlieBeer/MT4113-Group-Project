# -----------------------------------
# Helper Function: Tolerance Checker
# -----------------------------------
#' @description
#' Computes the ratio of the Euclidean norm of 2 vectors
#'
#' @param delta step value for each parameter
#' @param theta estimate of the minimised values
#'
#' @returns the ratio of the Euclidean norm of the inputted vectors
#' @export

Norm_Ratio <- function(delta, theta){
  return((sqrt(sum(delta^2))) / (sqrt(sum(theta^2))))
}

# -----------------------------------
# Multivariate Newton's Method (MVN)
# -----------------------------------
#' @description
#' Performs multivariate optimisation using Newton's method
#'
#' @param f function to be minimised
#' @param inits vector of initial values to be minimised
#' @param data (optional) daatframe with appropriate number of columns for number of variables being optimised
#' @param minimum search for minimum or maximum (takes TRUE if minimum search, FALSE if maximum search)
#' @param tol tolerance level
#' @param maxit maximum number of iterations run before stopping
#' @param method identifier, takes "MVN" only to allow the parent function to call on it
#' @param gradfn (optional) gradient function of f, uses finite differencing otherwise
#' @param hessfn (optional) hessian function of f, uses finite differencing otherwise
#' @param jacobfn NULL variable, included to match parent function inputs
#'
#' @returns A list containing the optimised estimate, the function evaluated at said estimate, the gradient of the function at the time, the tolerance level, whether the optimisation converged and the number of iterations ran.
#' @export

MVN <- function(f, inits, data = NULL, minimum = TRUE, tol, maxit,
                method = "MVN", gradfn = NULL, hessfn = NULL, jacobfn = NULL){
  # If searching for a maximum, work with the negative of the function, as Newton searches for minimums
  if (minimum == FALSE) {
    f <- -f}

  # If not provided with analytic functions, use numDeriv to estimate
  if (is.null(gradfn)){
    gradfn <- function(theta, data) {grad(f, theta, data = data)}
  }
  if (is.null(hessfn)){
    hessfn <- function(theta, data) {hessian(f, theta, data = data)}
  }

  theta <- inits # Set staring values
  # Set up some starting values for the loop
  niter <- 0
  loop <- TRUE

  while(loop){
    niter <- niter + 1
    # Compute gradient and hessian to update estimate
    g <- gradfn(theta, data)
    H <- hessfn(theta, data)

    delta <- solve(H, -g) # Compute step

    # Check stopping criteria
    if (Norm_Ratio(delta, theta) < tol | niter > maxit) loop <- FALSE
    theta <- theta + delta
  }
  # Check convergence
  # 0: tolerance reached, 1: tolerance not reached, 2: max. iterations reached
  if (NormRatio(delta, theta) < tol) {
    conv <- 0
  } else if (niter > maxit) {
    conv <- 2
  } else {
    conv <- 1
  }

  # Flip back over for final evaluation if maximising
  if (minimum == FALSE) {
    f <- -f}
  return(list(estimate = theta,
              feval = f(theta, data),
              grad = g,
              tolerance = Norm_Ratio(delta, theta),
              conv = conv,
              niter = niter))
}
