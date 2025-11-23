# -----------------------------------
# Bisection
# -----------------------------------
#' @name BS
#' @title Bisection Method
#' @description
#'The Bisection method works by taking two initial points, finding their derivatives and then if they
#'are opposite signs we look at the midpoint, repeating until we find where the derivative is zero
#'or within tolerance.
#'This method requires a smooth function.
#' @param f Function to be optimised
#' @param inits A length-2 numeric vector giving the initial search interval.
#' The gradient at these two points must have opposite signs.
#' @param data (optional) Dataframe passed to the function if needed.
#' @param minimum (optional) Search for minimum or maximum (takes TRUE if minimum search, FALSE if maximum search)
#' @param tol (optional) Tolerance level
#' @param maxit (optional) Maximum number of iterations run before stopping
#' @param method (optional) NULL variable, included to match parent function inputs
#' @param gradfn (optional) NULL variable, included to match parent function inputs
#' @param hessfn (optional) NULL variable, included to match parent function inputs
#' @param jacobfn (optional) NULL variable, included to match parent function inputs
#'
#' @return A list containing:
#' \item{estimate}{Optimised estimate}
#' \item{feval}{Function evaluated at optimised estimate}
#' \item{grad}{Gradient of function at optimised estimate}
#' \item{tolerance}{Final interval width}
#' \item{conv}{Whether or not the optimisation converged. 0 - converged, 1 - failed, 2 - max iterations reached}
#' \item{niter}{Number of iterations}
#' @references
#' Swallow, B. (2025). Univariate Optimization.
#' University of St Andrews.
#' \url{https://moody.st-andrews.ac.uk/moodle/pluginfile.php/2128840/mod_resource/content/4/_book/univariate-optimization.html#bisection-method}
#' @author Charles Beer
#' @examples
#' test_f<-function(x){
#' return (x^5-3*x^4+12*x^3-5*x^2-3)
#' }
#' test_inits<-c(0.1,5)
#' BS(test_f,test_inits)
#'
#' @importFrom numDeriv grad
#'@export
BS<-function(f, inits, data=NULL, minimum=TRUE, tol=1e-8, maxit=100,
             method=BS, gradfn=NULL, hessfn=NULL, jacobfn=NULL){
  if (is.null(data)) {
    func<-function(x) f(x)
  } else {
    func<-function(x) f(x, data)
  }
  right<-max(inits)
  left<-min(inits)
  if (right==left){stop("The input is not a range",call.=FALSE)}
  width<-right - left
  grad_left<- numDeriv::grad(func,left)
  grad_right<- numDeriv::grad(func,right)
  if (grad_left*grad_right>0){stop("the gradients of the range values have the same sign, Bisection will not work",call.=FALSE)}
  niter=0
  while (width>tol){
    niter=niter+1
    if (niter>maxit){
      estimate<-((left+right)/2)
      return (list(
        estimate=estimate,
        fval=func(estimate),
        grad=numDeriv::grad(func,estimate),
        tolerance=tol,
        conv=2,
        niter=niter
      ))

    }
    midpoint<-(left+right)/2
    grad_midpoint<-numDeriv::grad(func,midpoint)
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
    feval=func(estimate),
    grad=numDeriv::grad(func,estimate),
    tolerance=width,
    conv=0,
    niter=niter
      ))

}


#test
test_f<-function(x){
  return (x^5-3*x^4+12*x^3-5*x^2-3)
}

test_inits<-c(0.1,5)

BS(test_f,test_inits)

