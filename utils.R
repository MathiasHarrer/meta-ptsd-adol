hedgesG = function(d, N){
  return(d*(1-(3/(4*N-9))))
}

hdi_lo = function(x){
  require(HDInterval)
  return(hdi(x)[1,])
}