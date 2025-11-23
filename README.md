# MT4113-Group-Project  
## Description of functions:  

### BS:
The Bisection method works by taking two initial points and finding their derivatives, if they
are opposite signs we look at the midpoint, repeating until we find where the derivative is either zero
or within tolerance.
This method requires a smooth function.

### GN: 
Performs non-linear least squares optimisation using Gauss-Newton method.

### GS: 
The grid search method works by evaluating a function on points of a grid and picking
the point which has either the maximum or minimum value depending on the given criteria

### MVN:  
Performs multivariate optimisation using Newton's method.

### UVN:  
Performs univariate optimisation using Newton's method

### FunOptim:
Can call any of the aforementioned optimisation methods in order to optimise an inputted function.
