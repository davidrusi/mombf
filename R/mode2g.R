###
### mode2g.R
###

mode2g <- function(prior.mode,prior='iMom',nu=1,dim=1) {
  if (!(prior %in% c('normalMom','tMom','iMom'))) stop("Currently only prior=='normalMom', 'tMom' or 'iMom' are implemented")
  if (prior=='normalMom') {
    return(prior.mode/2)
  } else if (prior=='tMom') {
    if (nu<3) stop('The tMom prior must have nu>2 degrees of freedom')
    return(prior.mode*(nu-2+dim)/(2*nu))    
  } else if (prior=='iMom') {
    return(prior.mode *(nu+dim)/2)
  }
}
