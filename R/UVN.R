# the function of Univariate Newton Method new
UVN <- function(f, inits, data = NULL, minimum = TRUE, tol, maxit,
                method = "UVN", gradfn = NULL, hessfn = NULL, jacobfn = NULL) {

  # change the sign for minimum
  if (minimum == FALSE) {
    f <- -f
  }

  # set gradfn/grad and hessfn/hessian
  if (is.null(gradfn)) {
    gradfn <- function(theta, data) {
      numDeriv::grad(f, theta, data = data)
    }
  }
  if (is.null(hessfn)) {
    hessfn <- function(theta, data) {
      hessian(f, theta, data = data)
    }
  }

  # setup for the parameters
  par  <- inits
  step <- tol + 1
  iter <- 0

  # calculation
  while (abs(step) > tol && iter < maxit) {

    g     <- gradfn(par, data)
    H_mat <- hessfn(par, data)
    h     <- as.numeric(H_mat[1, 1])

    if(abs(h) < 1e-10) break

    step <- g / h
    par  <- par - step
    iter <- iter + 1
  }

  if (abs(step) <= tol) {
    conv <- 0
  } else if (iter >= maxit) {
    conv <- 1
  } else {
    conv <- 2
  }

  if (minimum == FALSE) {
    f <- -f
  }

  estimate <- par
  feval    <- f(estimate, data)

  result <- list(
    estimate  = estimate,
    feval     = feval,
    grad      = g,
    tolerance = abs(step),
    conv      = conv,
    niter     = iter
  )

  return(result)
}
