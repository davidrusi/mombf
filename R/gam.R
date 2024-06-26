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
    ans= bsplineCI(as.double(x), as.integer(degree), as.double(knots))
    #ans= .Call("bsplineCI", as.double(x), as.integer(degree), as.double(knots))
    return(matrix(ans,nrow=length(x),byrow=TRUE))
}



tensorbspline <- function(x, degree, knots, maineffects) {
    #Create tensor product B-splines
    # - x: matrix with covariate values
    # - degree: degree of each marginal B-spline. If a single number is given, it is used for all columns of x
    # - knots: list of length equal to ncol(x) indicating the knots for each marginal B-spline. Alternatively, a single vector with the common knots to be used for all columns of x
    # - maineffects: if TRUE answer is returned in terms of intercept + main effects + tensor, where tensor is orthogonalized with respect to the main effects. If maineffects==FALSE no main effects are returned, only a tensor that already accounts for the main effects
    # Output: list with the following elements
    # - main: intercept and main effects given by the univariate splines associated to each column of x
    # - tensor: tensor product B-splines
    # Note: the element tensor already includes terms for the main effects. If you wish to include 
    if (missing(degree)) stop("degree must be specified")
    if (missing(knots)) stop("knots must be specified")
    if (is.vector(x)) x= matrix(x, ncol=1)
    p= ncol(x)
    if (length(degree)==1) degree= rep(degree,p)
    if (!is.list(knots) & !is.vector(knots)) stop("knots should either be a list or a vector")
    if (is.list(knots) & (length(knots) != ncol(x))) stop("If knots is a list, it should have length equal to ncol(x)")
    des= vector('list',p)
    for (j in 1:p) {
        if (is.list(knots)) {
            des[[j]]= bspline(x[,j], degree=degree[j], knots=knots[[j]])
        } else {
            des[[j]]= bspline(x[,j], degree=degree[j], knots=knots)
        }
        des[[j]]= des[[j]][,apply(des[[j]],2,'sd')>0]
        colnames(des[[j]])= paste('x',j,'.',1:ncol(des[[j]]),sep='')
    }
    f= paste(paste('des[[',1:p,']]',sep=''), collapse=':')
    f= as.formula(paste('~ -1 +',f))
    tensor= model.matrix(f)
    colnames(tensor)= gsub('des\\[\\[.\\]\\]','',colnames(tensor))
    if (!maineffects) {
        main= NULL
    } else {
        #remove redundant columns from main effects
        main= cbind(1,do.call(cbind, des))
        colnames(main)[1]= 'intercept'
        myqr= qr(main)
        main= main[,myqr$pivot[1:myqr$rank]]
        #remove redundant columns from tensor
        tensor= residuals(lm(tensor ~ main))
        myqr= qr(tensor)
        tensor= tensor[,myqr$pivot[1:myqr$rank],drop=FALSE]
    }
    ans= list(main=main, tensor=tensor)
    return(ans)
}



#Obtain ICAR penalty matrix D taking the difference between each element and its neighbours
#
# If z is p x 1, then D z is also p x 1 and its j^th entry is z[j] - mean(neighbours of j)
#
# INPUT
# - neighbours: list where entry j returns the indexes of the neighbours of j
# - scale: if TRUE, the entries in D are divided by a constant such that tr(t(D) %*% D)= p
#
# OUTPUT: matrix D of dimension length(neighbours) x length(neighbours)
#         If scale=FALSE, d[j,j]= 1 and d[j,i]= -1/N[j] if i is a neighbour of j, where N[j] is the number of neighbours of j
#         If scale= TRUE, D is multiplied by p / tr(t(D) %*% D)
icar_dmatrix= function(neighbours, scale=TRUE) {
    D= diag(length(neighbours))
    N= sapply(neighbours, length)
    for (j in 1:nrow(D)) {
        if (N[j] > 0) {
            D[j, neighbours[[j]]] = -1/N[j]
        }
    }
    if (scale) {
        D= D * sqrt(ncol(D) / sum(colSums(D^2)))
    }
    return(D)
}
