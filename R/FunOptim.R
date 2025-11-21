# -----------------------------------
# Parent Optimisation Function
# -----------------------------------
#' @description
#' Unified function for calling the optimisation methods coded
#'
#' @param method optimisation method, takes one of "GS", "BS", "UVN", "MVN", "GN"
#' @returns A list containing the optimised estimate, the function evaluated at said estimate, the gradient of the function at the time, the tolerance level, whether the optimisation converged and the number of iterations ran.
#' @export

funoptim <- function(f, inits=NULL, data=NULL, minimum=TRUE, tol=1e-8, maxit=100,
                     method=NULL, gradfn=NULL, hessfn=NULL, jacobfn=NULL) {
  # Error Checking:
  if(!is.function(f))
    stop("Error: 'f' must be a function.")
  if(!is.numeric(inits))
    stop("Error: 'inits' must be a numeric vector.")
  if(!identical(minimum, TRUE) && !identical(minimum, FALSE))
    stop("Error: 'minimum' must be either TRUE or FALSE")
  if(!method %in% c("GS", "BS", "UVN", "MVN", "GN"))
    stop("Error: method provided is invalid.")

  if (method=="GS"){
    return(GS(f, inits, data, minimum, tol, maxit,
              method, gradfn, hessfn, jacobfn))
  }
  if (method=="BS"){
    return(BS(f, inits, data, minimum, tol, maxit,
              method, gradfn, hessfn, jacobfn))
  }
  if (method=="UVN"){
    return(UVN(f, inits, data, minimum, tol, maxit,
              method, gradfn, hessfn, jacobfn))
  }
  if (method=="MVN"){
    return(MVN(f, inits, data, minimum, tol, maxit,
              method, gradfn, hessfn, jacobfn))
  }
  if (method=="GN"){
    return(GN(f, inits, data, minimum, tol, maxit,
              method, gradfn, hessfn, jacobfn))
  }

  stop("method not named correctly")

}



