require(ncdf)
#require(ncdf, lib.loc="/local/R/lib")


yflip <- function(a)
  {
    (a[1:dim(a)[1],rev(1:dim(a)[2])])
  }

nc <- open.ncdf("pigv5.1km.cellcentered.nc")
velx <- -yflip(get.var.ncdf(nc,"us"))
vely <- -yflip(get.var.ncdf(nc,"vs"))
velc <- yflip(get.var.ncdf(nc,"ucoef"))
topg <- yflip(get.var.ncdf(nc,"topg"))
thck <- yflip(get.var.ncdf(nc,"thk"))
btrc <- yflip(get.var.ncdf(nc,"beta"))
ucoef <- yflip(get.var.ncdf(nc,"ucoef"))
divuh <- 0*btrc + 10
close.ncdf(nc)

#require(hdf5, lib.loc="/local/R/lib")

mcoarsen <- function(a,n)
  {
    ifelse(coarsen(a,n) < 1, 0, 1)
  }


coarsen <- function(a,n)
  {
    if (n ==0)
      {
        a
      } else {
        nx <- dim(a)[1]
        ny <- dim(a)[2]
        
        sqx <- seq(1,(nx-1),by=2)
        sqy <- seq(1,(ny-1),by=2)
        
        b <- (a[sqx,sqy] + a[sqx+1,sqy] + a[sqx,sqy+1] + a[sqx+1,sqy+1])/4

        coarsen(b,n-1)
        
      }
  }

sav <- function(ncrse,fnam)
  {

    #btrc <- ifelse(btrc < 500, 100, 1e+4)
    
    velx <- coarsen(velx,ncrse)
    vely <- coarsen(vely,ncrse)
    topg <- coarsen(topg,ncrse)
    thck <- coarsen(thck,ncrse)

    
    
    btrc <- coarsen(btrc,ncrse)
    divuh <- coarsen(divuh,ncrse)
    velc <- mcoarsen( velc, ncrse)

    rhoi <- 918
    rhoo <- 1028
    sg <- (topg + thck)
    sf <- (1-rhoi/rhoo)*thck
    s <- ifelse(sg > sf, sg, sf)
    f <- ifelse(sg > sf, 0, 1)

    velc <- ifelse(f == 0,velc,0)

    m <- ifelse(thck > 0, f, 2)
    image(m)
    contour(velc,add=TRUE)
    
    #vely <- ifelse(velc > 0, vely, 10000)
    #velx <- ifelse(velc > 0, velx, 10000)
    #hdf5save(fnam,"thck","topg","vely","velx", "velc", "btrc","divuh")

    
    

 
  }

#sav(2,"pigv5.4km.hdf5")
#sav(1,"pigv5.2km.hdf5")
sav(0,"pigv5.1km.hdf5")

#hdf5save("pigv5.1km.hdf5","thck","topg","vely","velx", "btrc","divuh")
  
#rhoi <- 918
#g <- 9.81
#rhoo <- 1028

#sg <- (topg + thck)
#sf <- (1-rhoi/rhoo)*thck
#s <- ifelse(sg > sf, sg, sf)
#f <- ifelse(sg > sf, 0, 1)
#m <- ifelse(thck > 0, f, 2)

#image(m)

#any(m==1 & thck==0)

#p <- hdf5load("pig.fwd.a.hdf5",load=FALSE)
#ucol <- colorRampPalette(c("white","blue","yellow"),space="Lab")(64)
#image(log(1+sqrt(p$vely^2+p$velx^2)),col=ucol)
#../flatToAMR2d.Linux.64.g++.gfortran.DEBUG.OPT.ex pigv5.4km.hdf5 pigv5.4km.2d.hdf5 btrc divuh thck topg velc velx vely
