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
  width<-right - left
  grad_left<- grad(f,left)
  grad_right<- grad(f,right)





}
