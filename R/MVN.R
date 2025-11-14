# Need to input the grad and hessian functions from numDeriv

# Helper functions:

# Computes the ratio of the Euclidean norm of 2 vectors
# INPUTS:
# delta: step value for each parameter
# theta: estimate of minimised values
Norm_Ratio <- function(delta, theta){
  return((sqrt(sum(delta^2))) / (sqrt(sum(theta^2))))
}

# INPUTS:
# f: function to be minimised
# inits: vector of initial values to search through
# data: (optional) dataframe with appropriate number of columns for number of
# variables being optimised
# minimum: seach for minimum or maximum (takes values TRUE or FALSE)
# tol: tolerance level
# maxit: max number of iterations run before stopping
# method: MVN, to call the parent function
# gradfn:(optional) gradient function of f, uses finite differencing otherwise
# hessfn: (optional) hessian function of f, uses finite differencing otherwise
# jacobfn: NULL variable, included to match parent function inputs


MVN <- function(f, inits, data = NULL, minimum = TRUE, tol, maxit,
                method = "MVN", gradfn = NULL, hessfn = NULL, jacobfn = NULL){
  # If searching fro a maximum, work with the negative of the function, as Newton searches for minimums
  if (minimum == FALSE) {
    f <- -f}

  # If not provided with analytic fucntions, use numderiv to estimate
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
  # 0: tolerance reached, 1: tolerance reached, 2: max. iterationa reached
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
