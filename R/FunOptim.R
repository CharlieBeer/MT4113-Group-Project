funoptim <- function(f, inits=NULL, data=NULL, minimum=TRUE, tol=1e-8, maxit=100,
                     method=NULL, gradfn=NULL, hessfn=NULL, jacobfn=NULL) {

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
test_f<-function(x){
  return (x^2-3)
}

test_inits<-c(0.06,0.5)
test_data<-read.csv("Data/testdata.csv")
funoptim(test_f,test_inits,method="GS")
funoptim(test_f,test_inits,method="BS")
