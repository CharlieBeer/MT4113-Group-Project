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

  initial<-inits[1]
  grid<-seq(initial-10,initial+10,length.out=100) #the grid that we will evaluate the function on

  fgrid<-sapply(grid,function(x) f(x,data)) #the function evaluated at the grid points

  index<-if (minimum) which.min(fgrid) else which.max(fgrid) #the index where the optimal value is
  best_arg<-as.numeric(grid[index]) #the argument values that optimize the function approximation
  best_val<-fgrid[index] #the value of the function at the min/max approximation

  return(list(
    estimate=best_arg,
    fval=best_val,
    grad=NA,
    tolerance=tol,
    conv=0,
    niter=length(grid)
  ))



}
