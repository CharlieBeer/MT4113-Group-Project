# MT4113-Group-Project  
## Description of functions:  

### BS:
The Bisection method works by taking two initial points and finding their derivatives, if they
are opposite signs we look at the midpoint, repeating until we find where the derivative is either zero
or within tolerance.
This method requires a smooth function.

### GN: 
The Gauss-Newton method works by iteratively solving a linear system derived from the initial problem, and updating the vector of parameters in accordance, until the parameter estimates converge completely or are within tolerance. This function requires a dataframe.

### GS: 
The grid search method works by evaluating a function on points of a grid and picking
the point which has either the maximum or minimum value depending on the given criteria

### MVN:  
Performs multivariate optimisation using Newton's method.
The multivariate Newton method uses the gradient vector and Hessian matrix to expand Newton's method to multiple parameters.

### UVN:
The univariate Newton method starts from an initial point and iterates by first and second derivatives. 
The Newton iterations determine a local minimum or maximum at each step.
This process converges to either the minimum or maximum.

### FunOptim:
Can call any of the aforementioned optimisation methods in order to optimise an inputted function.
