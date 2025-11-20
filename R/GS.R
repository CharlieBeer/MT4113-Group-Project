GS<-function(f, inits, data=NULL, minimum=TRUE, tol=1e-8, maxit=100,
             method=NULL, gradfn=NULL, hessfn=NULL, jacobfn=NULL){

  initial<-inits[1]
  grid<-seq(initial-10,initial+10,length.out=100) #the grid that we will evaluate the function on

  fgrid<-sapply(grid,f) #the function evaluated at the grid points

  index<-if (minimum) which.min(fgrid) else which.max(fgrid) #the index where the optimal value is
  best_arg<-as.numeric(grid[index]) #the argument values that optimize the function approximation
  best_val<-fgrid[index] #the value of the function at the min/max approximation

  return(list(
    estimate=best_arg,
    fval=best_val,
    grad=NA,
    tolerance=tol,
    conv=0,
    niter=nrow(grid)
  ))



}
