
testfunction= function(x) {
    ans= .Call("testfunctionCI",as.double(x));
    return(ans);
}

