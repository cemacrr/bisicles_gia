#require(hdf5, lib.loc="/local/R/lib")
require(hdf5)
p <- hdf5load("fwd.a.hdf5",load=FALSE)
sg <- p$topg+p$thck
sf <- (1-910/1028) * p$thck
s <- ifelse(sg > sf, sg, sf)
grounded <- ifelse(sg > sf, 1, 0)
sav <- function(fnam,sd)
  {
    n <- length(p$thck)
    nf <- function()
      {
        (1 + rnorm(n,sd=sd))
      }
    
    thck <- p$thck * nf()
    topg <- p$topg * nf()

    nsg <- topg+thck
    nsf <- (1-910/1028) * thck
    ngrounded <- ifelse(nsg > nsf, 1, 0)

    thck <- ifelse(ngrounded==grounded,thck,p$thck)
    topg <- ifelse(ngrounded==grounded,topg,p$topg) 

    nsg <- topg+thck
    nsf <- (1-910/1028) * thck
    ngrounded <- ifelse(nsg > nsf, 1, 0)
    par(mfrow=c(2,2))
    
   
    
    velx <- p$velx * nf()
    vely <- p$vely * nf()
    btrc <- p$btrc * 0 + 300
    divuh <- p$divuh * 0.0
    hdf5save(fnam,"thck","topg","velx","vely", "btrc","divuh")
  }

#sav("ctrl.a.noise0.hdf5")
#sav("ctrl.a.noise0.01.hdf5",0.01)
sav("ctrl.a.pgl.noise0.03.hdf5",0.00)


#f <- function(x){cos(2*pi*x)+sin(2*pi*x)}
#curve(f,-1,2)
#abline(v=0)
#abline(v=1)
