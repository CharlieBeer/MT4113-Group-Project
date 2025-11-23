# -----------------------------------
# Grid Search
# -----------------------------------
#' @name GS
#' @title Grid search Method
#' @description
#'The grid search method works by evaluating a function on points of a grid and picking
#'the point which has either the maximum or minimum value depending on the given criteria
#' @param f Function to be optimised
#' @param inits An initial guess, the grid will be on +- 10 of this point
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
#' \item{tolerance}{Tolerance level reached through optimisation}
#' \item{conv}{Whether or not the optimisation converged. 0 - converged, 1 - failed, 2 - max iterations reached}
#' \item{niter}{Number of iterations}
#' @references
#' Swallow, B. (2025). Univariate Optimisation.
#' University of St Andrews.
#' \url{https://moody.st-andrews.ac.uk/moodle/pluginfile.php/2128840/mod_resource/content/4/_book/univariate-optimization.html#grid-search}
#' @author Charles Beer
#' @examples
#' test_f<-function(x){
#' return ((x-0.25)*(x-0.6)*(x-5)*(x-3))
#' }
#' test_inits<-c(0)
#' GS(test_f,test_inits,minimum=TRUE)
#'
#'@export

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

