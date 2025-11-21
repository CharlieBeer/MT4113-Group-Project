<<<<<<< HEAD
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
=======
# the function of Univariate Newton Method new
UVN <- function(f, inits, data = NULL, minimum = TRUE, tol, maxit,
                method = "UVN", gradfn = NULL, hessfn = NULL, jacobfn = NULL) {

  # change the sign for minimum
  if (minium == FALSE) {
    f <- -f
  }

  # set gradfn/grad and hessfn/hessian
  if (is.null(gradfn)) {
    gradfn <- function(theta, data) {
      grad(f, theta, data = data)
    }
  }
  if (is.null(hessfn)) {
    hessfn <- function(theta, data) {
      hessian(f, theta, data = data)
    }
  }

  # setup for the parameters
  par  <- inits
>>>>>>> origin/main
  step <- tol + 1
  iter <- 0

  # calculation
  while (abs(step) > tol && iter < maxit) {

<<<<<<< HEAD
   g <- grad(target, par)
   H_mat <- hessian(target, par)
   h     <- H_mat[1, 1]

    if(abs(h) < 1e-10) {
      break
    }
    step <- g / h
    par <- par - step
=======
    g     <- gradfn(par, data)
    H_mat <- hessfn(par, data)
    h     <- as.numeric(H_mat[1, 1])

    if(abs(h) < 1e-10) break

    step <- g / h
    par  <- par - step
>>>>>>> origin/main
    iter <- iter + 1
  }

  if (abs(step) <= tol) {
    conv <- 0
  } else if (iter >= maxit) {
    conv <- 1
  } else {
    conv <- 2
  }

<<<<<<< HEAD
  estimate <- par
  feval <- f(estimate, data)
  grad_final <- grad(f_tar, estimate)
  tolerance <- abs(step)
  niter <- iter
=======
  if (minimum == FALSE) {
    f <- -f
  }

  estimate <- par
  feval    <- f(estimate, data)
>>>>>>> origin/main

  result <- list(
    estimate  = estimate,
    feval     = feval,
<<<<<<< HEAD
    grad      = grad_final,
    tolerance = tolerance,
    conv      = conv,
    niter     = niter
  )

  return(result)

=======
    grad      = g,
    tolerance = abs(step),
    conv      = conv,
    niter     = iter
  )

  return(result)
>>>>>>> origin/main
}
