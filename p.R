p = function(N) {
  # dennie allen "Sat Aug  9 21:40:36 2014"
  # calc partitons of N by 10's using the Hardy-Ramanujan formula
  # inspired by Ken Ono, Emory U, 'New Theories Reveal the Nature of Numbers'
  df = data.frame(n = seq(0, N, by = 10))
  P = data.frame()
  df <- within(df, P <- sapply(df,function(n){floor(1/(4*n*sqrt(3)) * 		exp(pi*sqrt(2*n/3)))}))
  colnames(df) <- c("n", "P(n)")
  return(data.matrix(df))  # finally! data.matrix() displays colnames.
}
p(100)
