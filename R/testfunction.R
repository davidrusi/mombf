
testfunction= function(x) {
    ans= testfunctionCI(as.double(x))
    #ans= .Call("testfunctionCI",as.double(x));
    return(ans);
}

