#Bisection Method
#The Bisection method works by taking two initial points, finding their derivatives and then if they
#are opposite signs we look at the midpoint, repeating until we find where the derivative is zero
#or within tolerance.
#this method requires a smooth function
#'@importFrom numDeriv grad
BS<-function(f, inits, data=NULL, minimum=TRUE, tol=1e-8, maxit=100,
             method=NULL, gradfn=NULL, hessfn=NULL, jacobfn=NULL){
  right<-max(inits)
  left<-min(inits)
  if (right==left){stop("The input is not a range",call.=FALSE)}
  width<-right - left
  grad_left<- grad(f,left)
  grad_right<- grad(f,right)
  if (grad_left*grad_right>0){stop("the gradients of the range values have the same sign, Bisection will not work",call.=FALSE)}
  while (width>tol){
    midpoint<-(left+right)/2
    grad_midpoint<-grad(f,midpoint)
    if (grad_midpoint*grad_left>0){
      left<-midpoint
      grad_left<-grad_midpoint
    }else{
      right<-midpoint
      grad_right<-grad_midpoint
    }
    width<-right-left
  }
  return ((left+right)/2)

}


#test
test_f<-function(x){
  return (x^5-3*x^4+12*x^3-5*x^2-3)
}

test_inits<-c(0.06,0.5)

BS(test_f,test_inits)
