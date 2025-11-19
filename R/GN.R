GN <- function(f, inits, data, minimum = TRUE, tol = 1e-10, maxit = 1000, 
                       method = "GN", gradfn = NULL, hessfn = NULL, jacobfn = NULL) {
  
  #need advice on how to load Jacobian package desperately
  
  #compute initial residuals
  resids <- function(theta, data) {
    preds <- f(theta, data)
    r <- data$y - preds
    return(r)
  }
  
  #jacobian
  if (is.null(jacobfn)) {
    jacobian_function <- function(theta, data) {
      jacobian(resids, theta, data = data)
    }
  } else {
    jacobian_function <- jacobfn
  }
  
  #compute Euclidean norm
  norm <- function(x) {
    return(sqrt(sum(x^2)))
  }
  
  #parameters for the loop
  theta_current <- inits
  iter <- 0
  conv_code <- 2  #default to max iterations reached - this can be updated if something else happens
  final_tolerance <- NA
  
  #begin loop
  for (iter in 1:maxit) {
    
    #compute Jacobian
    J <- jacobian_function(theta_current, data)
    
    #compute residuals
    residuals <- resids(theta_current, data)
    
    #compute gradient
    gradient <- t(J) %*% residuals
    
    #compute Hessian
    hessian <- t(J) %*% J
    
    #solve for the parameter update
    delta <- solve(hessian, -gradient)
    
    #update theta
    theta_new <- theta_current + delta
    
    #check the convergence - is it below the tolerance yet?
    increment <- norm(delta) / norm(theta_current)
    final_tolerance <- increment  #keep track of tolerance
    
    if (increment < tol) { #convergence has been reached
      conv_code <- 0  #update convergence code
      break
    }
    
    #update theta for the next iteration
    theta_current <- theta_new
  }

  
  #compute the final residuals and final value of the function
  final_residuals <- resids(theta_current, data)
  feval <- sum(final_residuals^2)
  
  #results
  result <- list(estimate = theta_current, feval = feval, grad = gradient, tolerance = final_tolerance, conv = conv_code, niter = iter)
  
  return(result)
}
