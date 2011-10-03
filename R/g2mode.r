g2mode <- function(g,prior='iMom',nu=1,dim=1) {
  if (!(prior %in% c('normalMom','tMom','iMom'))) stop("Currently only prior=='normalMom', 'tMom' or 'iMom' are implemented")
  if (prior=='normalMom') {
    return(2*g)
  } else if (prior=='tMom') {
    return(g*2*nu/(nu-2+dim))
  } else if (prior=='iMom') {
    return(2*g/(nu+dim))
  }
}
