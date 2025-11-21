# -----------------------------------
# Grid Search
# -----------------------------------
#' @description
#' Performs Grid Search
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

GS<-function(f, inits, data=NULL, minimum=TRUE, tol=1e-8, maxit=100,
             method=NULL, gradfn=NULL, hessfn=NULL, jacobfn=NULL){

  if (is.null(data)) {
    func<-function(x) f(x)
  } else {
    func<-function(x) f(x, data)
  }

  initial<-inits[1]
  grid<-seq(initial-10,initial+10,length.out=100) #the grid that we will evaluate the function on

  fgrid<-sapply(grid,func) #the function evaluated at the grid points

  index<-if (minimum) which.min(fgrid) else which.max(fgrid) #the index where the optimal value is
  best_arg<-as.numeric(grid[index]) #the argument value that optimizes the function approximation
  best_val<-fgrid[index] #the value of the function at the min/max approximation

  return(list(
    estimate=best_arg,
    feval=best_val,
    grad=NA, #grid search doesn't find a gradient
    tolerance=tol,
    conv=0,
    niter=length(grid)
  ))



}

test_f<-function(x){
  return ((x-0.25)*(x-0.6)*(x-5)*(x-3))
}

test_inits<-c(0,10)

GS(test_f,test_inits,minimum=TRUE)

