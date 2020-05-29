.noAddInfo <- function(nInterval, eggDuration, longevity, a_init){

  minimize <- function(param, ...){
    d    <- numeric(nInterval)
    d[1] <- eggDuration
    for(i in 2:nInterval){
      d[i] <- d[1]*param*(param+1)^(seq_len(nInterval)[i]-2)
    }
    min <- abs(longevity - sum(d))
    return(min)
  }

  a_param   <- newuoa(par = a_init, fn = minimize)$par

  delta    <- numeric(nInterval)
  delta[1] <- eggDuration
  for(i in 2:nInterval){
    delta[i] <- delta[1]*a_param*(a_param+1)^(seq_len(nInterval)[i]-2)
  }

  output <- list(delta = delta, a = a_param)

  return(output)
}


.AddInfo <- function(nInterval, eggDuration, longevity, a_init, addInfo){

  minimize <- function(param, ...){
    d    <- c(eggDuration, addInfo)
    for(i in 2:nInterval){
      if(!is.na(d[i])) next
      d[i] <- d[1]*param*(param+1)^(seq_len(nInterval)[i]-2)
    }
    min <- abs(longevity - sum(d))
    return(min)
  }

  a_param   <- newuoa(par = a_init, fn = minimize)$par

  delta    <- c(eggDuration, addInfo)
  for(i in 2:nInterval){
    if(!is.na(delta[i])) next
    delta[i] <- delta[1]*a_param*(a_param+1)^(seq_len(nInterval)[i]-2)
  }

  output <- list(delta = delta, a = a_param)

  return(output)
}
