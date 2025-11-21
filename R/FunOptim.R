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

#test

test_data<-read.csv("Data/testdata.csv")

test_f_GS<-function(x){
  return ((x-0.25)*(x-0.6)*(x-5)*(x-3))
}
test_inits_GS<-c(0)
funoptim(test_f_GS,test_inits_GS,method="GS",minimum=TRUE)


test_f_BS<-function(x){
  return ((x-0.25)*(x-0.6)*(x-5)*(x-3))
}
test_inits_BS<-c(0,10)
funoptim(test_f_BS,test_inits_BS,method="BS",minimum=TRUE)

install.packages("EnvStats")
install.packages("LaplacesDemon")
library(EnvStats)  #zmlnorm functions
library(LaplacesDemon)  #invlogit function

dln <- function(p, data) {

  return(-sum(log(dzmlnorm(data, p[1], exp(p[2]), invlogit(p[3])))))

}

funoptim(f = dln, inits, data = test_data$x, minimum, tol,
         maxit, method="MVN", gradfn, hessfn, jacobfn)



f1 <- function(x) (x-2)^2
init <- c(5)

funoptim(f1, init, method="GS")
funoptim(f1, c(0,5), method="BS")
funoptim(f1, 5, method="UVN")
funoptim(function(t) f1(t[1]), c(5), method="MVN")
