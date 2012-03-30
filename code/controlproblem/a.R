require(hdf5, lib.loc="/local/R/lib")
#require(hdf5)
n <- 64
L <- 100
m <- n+1
xf <- seq(0,L,l=m)
yf <- seq(0,L,l=m)

xc <- 0.5 * (xf[-1] + xf[-m])
yc <- 0.5 * (yf[-1] + yf[-m])

xcm <-matrix( rep(xc,n), ncol = n)
ycm <-matrix( rep(yc,each=n), ncol = n)
wx <- pi * xcm / L
wy <- pi * ycm / L

if (FALSE){
  topg <- t ( -500 *  sin(wx)*sin(wy))
  thck <- t ( 600 - 400 * sin(wx)*sin(wy) )
  thck <- ifelse(thck < 220,0,thck)
  btrc <- t ( 450 - 200 * (sin(wx) + sin(2*wx))*sin(wy))
} else {
  topg <- 0.5* t ( -500 *  cos(2*wx)*cos(2*wy))
  thck <- 0.5* t ( 600 - 400 * cos(2*wx)*cos(2*wy) )
  thck <- ifelse(thck < 120,0,thck)
  btrc <- t ( 450 - 250 * (cos(2*wx) + sin(2*wx))*cos(2*wy))

}





vely <- t( matrix(0 + rnorm(n*n,sd=0),n,n))
velx <- t( matrix(100 + rnorm(n*n,sd=0),n,n))
velx <- matrix(0,n,n)
vely <- matrix(0,n,n)
divuh <- matrix(0,n,n)

sg <- topg+thck
sf <- (1-910/1028) * thck
s <- ifelse(sg > sf, sg, sf)
b <- s - thck

C <- ifelse(sg > sf, btrc, 0)


hdf5save("a.hdf5","thck","topg","vely","velx", "btrc","divuh")
filled.contour(C,col=rainbow(64))
#contour(C,add=TRUE)
