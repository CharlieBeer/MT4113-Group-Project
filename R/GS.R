GS<-function(f, inits, data=NULL, minimum=TRUE, tol=1e-8, maxit=100,
             method=NULL, gradfn=NULL, hessfn=NULL, jacobfn=NULL){

  grids<-list() #we don't know how many parameters there are so we make an n-dimensional grid
  for (param in inits){
    grids[[length(grids)+1]]<-seq(param*0.5,param*2,length.out=100)
  }
  grid<-do.call(expand.grid,grids) #the grid that we will evaluate the function on

  fgrid<-apply(grid,1,function(x) f(x,data)) #the function evaluated at the grid points

  index<-if (minimum) which.min(fgrid) else which.max(fgrid) #the index where the optimal value is
  best_args<-as.numeric(grid[index,]) #the argument values that optimize the function approximation
  best_val<-fgrid[index] #the value of the function at the min/max approximation

  return(list(
    estimate=best_args,
    fval=best_val,
    grad=NA,
    tolerance=tol,
    conv=0,
    niter=nrow(grid)
  ))



}
