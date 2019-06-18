#####################################################################
##
## ROUTINES RELATED TO GENERALIZED ADDITIVE MODELS
##
#####################################################################


bspline <- function(x, degree, knots) {
 # B-spline basis eval at vector of values x.
 # Normalized to sum to 1 at any x value.
 #
 # Input:
 #    x: vector of values at which to evaluate the B-spline basis
 #    degree: degree of the spline (0: piecewise constant, 1: linear etc.)
 #    knots: vector with positions of the knots
 #  Output: matrix[nx][nknots-degree-1] containing the B-spline basis
    ans= .Call("bsplineCI", as.double(x), as.integer(degree), as.double(knots))
    return(matrix(ans,nrow=length(x),byrow=TRUE))
}
