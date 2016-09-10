# reconstruct a covariance of X matrix, from a QR decomposition R.

nc=2
nr=4
x =matrix(nrow=nr, ncol=nc, data=c(1,0,04,2,6,-1,0))

y =c(2,-2,2,-4)

qq=qr(x)
R=qr.R(qq)
Q = qr.Q(qq,complete=F) # n x p, ready to by multiplied with R (R is really R1, in that it omits the lower all-zeros portion of R)
Qfull=qr.Q(qq,complete=T) # n x n


x2 = qr.Q(qq) %*% R  # matches x

(t(R) %*% R)/3

m=(R %*% t(R))

Q %*% R

t(R) %*% t(Q) %*% Q %*% R

t(x) %*% x


u=rep(0,nc)
for (i in 1:nc) {
  u[i] = mean(x[,i])
}

# nr = nrow(X)
R_to_cov=function(R, u, nr) {

  umat=matrix(ncol=nc, nrow=nc,data=rep(0,nc*nc))

  for (i in 1:nc) {
    for (j in 1:nc) {
      umat[i,j] = u[i]*u[j]*nr
    }
  }

  #((t(x) %*% x) - umat)/(nr-1)

  (t(R) %*% R - umat)/(nr-1)
}



co2 = (t(R) %*% R - umat)/(nr-1)

co2/cov(x)



mycov=function(x) {
  nc=ncol(x)
  nr=nrow(x)
  # compute mean in u
  u=rep(0,nc)
  for (i in 1:nc) {
   u[i] = mean(x[,i])
  }
  # subtract mean
  xn = x
  for (i in 1:nc) {
   xn[,i] = xn[,i] - u[i]
  }

  # compute covariance matrix X'X
  co = t(xn) %*% xn / (nr-1)
  co
}
