funoptim <- function(f, inits, data, minimum, tol, maxit,
                     method, gradfn, hessfn, jacobfn) {

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
