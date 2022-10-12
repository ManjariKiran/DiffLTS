dist_diag <- function(buffer, loop, res, exclude){
# p <- pairdist(loop)
  bin <- p/res
  res1 <- res/1000
  start <- (bin)*res1
  m=(buffer*2)+1
  M <- matrix(data=NA,nrow=m,ncol=m)
for(k in 1:m){
  rem = k-1
  for(i in k:m){
    j = i-rem
      M[i,j] <- start
    }
    start <- start - res1
}
  start <- (bin)*res1
for(k in 1:m){
    rem = k-1
    for(j in k:m){
      i = j-rem
      M[i,j] <- start
    }
     start <- start + res1 
  }
#  M[is.na(M)] <- 1
  M[M < exclude] <- NA
  return(M)
  }
dist_diag(5, l,res, 0)
