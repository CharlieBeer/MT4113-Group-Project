library(numDeriv)

# the function of Univariate Newton Method
UVN <- function(f, inits, data, minimum, tol, maxit,
                method, gradfn, hessfn, jacobfn) {

  # setup for the target function
  f_tar <- function(theta){
    f(theta,data)
  }

  # find the maximum  and minimum
  if (minimum) {
    target <- function(theta) f_tar(theta)
  } else{
    target <- function(theta) -f_tar(theta)
  }

  # setup for the parameters
  par <- inits
  step <- tol + 1
  iter <- 0

  # calculation
  while (abs(step) > tol && iter < maxit) {

   g <- grad(target, par)
   H_mat <- hessian(target, par)
   h     <- H_mat[1, 1]

    if(abs(h) < 1e-10) {
      break
    }
    step <- g / h
    par <- par - step
    iter <- iter + 1
  }

  if (abs(step) <= tol) {
    conv <- 0
  } else if (iter >= maxit) {
    conv <- 1
  } else {
    conv <- 2
  }

  estimate <- par
  feval <- f(estimate, data)
  grad_final <- grad(f_tar, estimate)
  tolerance <- abs(step)
  niter <- iter

  result <- list(
    estimate  = estimate,
    feval     = feval,
    grad      = grad_final,
    tolerance = tolerance,
    conv      = conv,
    niter     = niter
  )

  return(result)

}
