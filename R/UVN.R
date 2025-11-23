# -----------------------------------
# Univariate Newton's Method (UVN)
# -----------------------------------
#' @description
#' Performs univariate optimisation using Newton's method
#'
#' @param f function to be minimised
#' @param inits initial value to be minimised
#' @param data (optional) single column dataframe of data in relation to the parameter to be optimised
#' @param minimum search for minimum or maximum (takes TRUE if minimum search, FALSE if maximum search)
#' @param tol tolerance level
#' @param maxit maximum number of iterations run before stopping
#' @param method identifier, takes "UVN" only to allow the parent function to call on it
#' @param gradfn (optional) gradient function of f, uses finite differencing otherwise
#' @param hessfn (optional) hessian function of f, uses finite differencing otherwise
#' @param jacobfn NULL variable, included to match parent function inputs
#'
#' @returns A list containing the optimised estimate, the function evaluated at said estimate, the gradient of the function at the time, the tolerance level, whether the optimisation converged and the number of iterations ran.
#' @export

UVN <- function(f, inits, data, minimum, tol, maxit,
                method = "UVN", gradfn, hessfn, jacobfn) {

  # Local wrapper function that allows passing of the data argument to grad and hessian, if applicable
  # ... included when Wrapper was a global function, kept for defensive programming in case of later expansion
  Wrapper <- function(theta, ...){
    if (is.null(data)){
      f(theta)
    } else {
      f(theta, data)
    }
  }

  # If searching for a maximum, work with the negative of the function, as Newton searches for minimums
  # Wrapper function follows -f
  if (minimum == FALSE) {
    og_f <- f
    f <- function(theta, ...){
      -og_f(theta, ...)
    }
    Wrapper <- function(theta, ...){
      if (is.null(data)){
        f(theta)
      } else {
        f(theta, data)
      }
    }
  }

  # If gradient/hessian is not provided, then use gradient/hessian from numDeriv
  if (is.null(gradfn)) {
    gradfn <- function(theta) {
      numDeriv::grad(Wrapper, theta)
    }
  }
  if (is.null(hessfn)) {
    hessfn <- function(theta) {
      numDeriv::hessian(Wrapper, theta)
    }
  }

  # Initialise parameter value, step size and iteration
  par  <- inits
  step <- tol + 1
  iter <- 0

  # Newton iterations
  while (abs(step) > tol && iter < maxit) {
    iter <- iter + 1

    g     <- gradfn(par)
    H     <- hessfn(par)
    h     <- as.numeric(H[1, 1])

    # Check stopping criteria
    if(abs(h) < 1e-10) break

    step <- g / h
    par  <- par - step
  }

  # Check convergence
  # 0: tolerance reached, 1: tolerance not reached, 2: max. iterations reached
  if (abs(step) <= tol) {
    conv <- 0
  } else if (iter >= maxit) {
    conv <- 2
  } else {
    conv <- 1
  }

  # Flip back over for final evaluation if maximising
  if (minimum == FALSE) {
   f <- og_f
  }

  result <- list(
    estimate  = par,
    feval     = Wrapper(par),
    grad      = g,
    tolerance = abs(step),
    conv      = conv,
    niter     = iter
  )

  return(result)
}

# =====================
# Test
# =====================
#library(tidyverse)
#MiningData <- read.csv("MiningData.csv") %>%
#  mutate(ratio = width/depth)
#alpha_0 <- median(MiningData$angle)
#beta_0 <- median(na.omit(-log(1 - (MiningData$angle)/alpha_0) / MiningData$ratio))
#RSS_UVN <- function(theta, data){
#  alpha <- theta
#  beta  <- beta_0
#  pred <- alpha * (1 - exp(-beta * data$ratio))
#  sum((data$angle - pred)^2)
#}
#UVN_test2 <- UVN(f = RSS_UVN, inits = alpha_0, data = MiningData %>% select(ratio, angle),
#                      minimum = TRUE, tol = 1e-6, maxit = 100, method = "UVN", gradfn = NULL,
#                      hessfn = NULL, jacobfn = NULL)
#UVN_test2
