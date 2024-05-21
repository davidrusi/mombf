
testfunction= function(A, oldcol, newcol) {
    A= Matrix::Matrix(A, sparse=TRUE)
    ans= testfunctionCI(A, as.integer(oldcol), as.integer(newcol))
    return(ans)
}

