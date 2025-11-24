#Test File

#BS
test_f<-function(x){
  return (x^5-3*x^4+12*x^3-5*x^2-3)
}
test_inits<-c(0.1,5)
funoptim(test_f,test_inits,method="BS")

#GN

#GS
test_f<-function(x){
  return ((x-0.25)*(x-0.6)*(x-5)*(x-3))
}
test_inits<-c(0)
funoptim(test_f,test_inits,minimum=TRUE,method="GS")
#MVN

#UVN
test_f1 <- funciton(theta){
  x <- theta[1]
  return((x-2)^2)
}
funoptim(test_f1, inits = 10, tol = 1e-6, maxit = 100, method="UVN")
