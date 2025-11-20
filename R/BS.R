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
# Bisection
# -----------------------------------
#' @description
#' Performs Bisection
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




#Bisection Method
#The Bisection method works by taking two initial points, finding their derivatives and then if they
#are opposite signs we look at the midpoint, repeating until we find where the derivative is zero
#or within tolerance.
#this method requires a smooth function
#'@importFrom numDeriv grad
BS<-function(f, inits, data=NULL, minimum=TRUE, tol=1e-8, maxit=100,
             method=NULL, gradfn=NULL, hessfn=NULL, jacobfn=NULL){
  func<-function(x) f(x, data)
  right<-max(inits)
  left<-min(inits)
  if (right==left){stop("The input is not a range",call.=FALSE)}
  width<-right - left
  grad_left<- grad(func,left)
  grad_right<- grad(func,right)
  if (grad_left*grad_right>0){stop("the gradients of the range values have the same sign, Bisection will not work",call.=FALSE)}
  niter=0
  while (width>tol){
    niter=niter+1
    if (niter>maxit){
      estimate<-((left+right)/2)
      return (list(
        estimate=estimate,
        fval=func(estimate),
        grad=grad(func,estimate),
        tolerance=tol,
        conv=2,
        niter=niter
      ))

    }
    midpoint<-(left+right)/2
    grad_midpoint<-grad(func,midpoint)
    if (grad_midpoint*grad_left>0){
      left<-midpoint
      grad_left<-grad_midpoint
    }else{
      right<-midpoint
      grad_right<-grad_midpoint
    }
    width<-right-left
  }
  estimate=(left+right)/2
  return (list(
    estimate=estimate,
    fval=func(estimate),
    grad=grad(func,estimate),
    tolerance=tol,
    conv=0,
    niter=niter
      ))

}


#test
test_f<-function(x){
  return (x^5-3*x^4+12*x^3-5*x^2-3)
}

test_inits<-c(0.06,0.5)

BS(test_f,test_inits)
