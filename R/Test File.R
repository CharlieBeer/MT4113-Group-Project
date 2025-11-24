#Test File

#BS
test_f<-function(x){
  return (x^5-3*x^4+12*x^3-5*x^2-3)
}
test_inits<-c(0.1,5)
funoptim(test_f,test_inits,method="BS")

#GN
test_f <- function(theta, data){
  a <- theta[1]
  b <- theta[2]
  return(a * exp(b * data$x))
}
test_data <- list(x = 1:10, y = 2.5 * exp(-1.5 * 1:10))
test_inits <- c(1, -1)

funoptim(test_f, test_inits, test_data, method="GN")

#GS
test_f<-function(x){
  return ((x-0.25)*(x-0.6)*(x-5)*(x-3))
}
test_inits<-c(0)
funoptim(test_f,test_inits,minimum=TRUE,method="GS")
#MVN
# No data:
f <- function(theta){
  x <- theta[1]
  y <- theta[2]
  return((x-3)^2 + (y+2)^2)
}

funoptim(f = f, inits = c(4, -1), method = "MVN")

# With data:
n<-50
x<-runif(n,0,2)
y<-2*exp(0.8*x)+rnorm(n,0,0.2)
test_data<-data.frame(x,y)

test_f<-function(theta,data){
  a<-theta[1]
  b<-theta[2]
  yhat<-a*exp(b*data$x)
  return(sum((data$y-yhat)^2))
}

test_inits<-c(1,1)

funoptim(test_f,test_inits,test_data,method="MVN")

#UVN
