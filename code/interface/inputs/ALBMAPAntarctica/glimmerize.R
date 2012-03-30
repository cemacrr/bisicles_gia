require(ncdf, lib.loc="/local/R/lib")
nc <- open.ncdf("Antarctica-5km.nc")
x <- get.var.ncdf(nc,"x0")
y <- get.var.ncdf(nc,"y0")
beta <- get.var.ncdf(nc,"beta")
topg <- get.var.ncdf(nc,"topg")
thck <- get.var.ncdf(nc,"thk")
close.ncdf(nc)

nx <- length(x)
ny <- length(y)

dx <- x[2]-x[1]
x1 <- c(x[1]-dx,x,x[nx]+dx)
x0 <- 0.5*(x1[1:(nx+1)] + x1[1:(nx+1)+1])
y1 <- c(y[1]-dx,y,y[ny]+dx)
y0 <- 0.5*(y1[1:(nx+1)] + y1[1:(nx+1)+1])

nctocc <- function(a)
  {
    n <- dim(a)[1]
    m <- dim(a)[2]

    sqn <- 2:n
    sqm <- 2:m
    
    0.25 * ( a[sqn,sqm] + a[sqn-1,sqm] + a[sqn-1,sqm-1] + a[sqn,sqm-1] )
   
  }

cctonc <- function(cc,pad)
  {
    row <- rep(pad,nx+2)
    column <- rep(pad,ny)
    rbind(row,
          cbind( column ,  cc , column),
          row)

  }


betanc <- cctonc(beta,20)
betacc <- nctocc(betanc)

#betacc <- 0*betanc + ifelse(x0 < 4500 | x0 > 6395500,1,0) 

thcknc <- cctonc(thck,0)
topgnc <- cctonc(topg,-1e+4)

ncxdim <- dim.def.ncdf("x1","m",x1)
ncydim <- dim.def.ncdf("y1","m",y1)
ccxdim <- dim.def.ncdf("x0","m",x0)
ccydim <- dim.def.ncdf("y0","m",y0)
timedim <- dim.def.ncdf("time","a",0)
#varx1 <- var.def.ncdf("x1","m",ncxdim,-9999,prec="double")
#vary1 <- var.def.ncdf("y1","m",ncydim,-9999,prec="double")

varbeta <- var.def.ncdf("beta","G",list(ccxdim,ccydim,timedim),
                        -9999,prec="double")
varthck <- var.def.ncdf("thk","m",list(ncxdim,ncydim,timedim),
                        -9999,prec="double")
vartopg <- var.def.ncdf("topg","m",list(ncxdim,ncydim,timedim),
                        -9999,prec="double")


ncs <- create.ncdf("Antarctica-glimmer-5km.nc", list(varbeta,varthck,vartopg))

put.var.ncdf(ncs,varbeta,betacc)
put.var.ncdf(ncs,varthck,thcknc)
put.var.ncdf(ncs,vartopg,topgnc)

close.ncdf(ncs)
