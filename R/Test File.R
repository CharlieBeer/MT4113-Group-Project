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
GS(test_f,test_inits,minimum=TRUE,method="GS")

#MVN
# No data:
f <- function(theta){
  x <- theta[1]
  y <- theta[2]
  return((x-3)^2 + (y+2)^2)
}

funoptim(f = f, inits = c(4, -1), method = "MVN")

# With data:


#UVN
