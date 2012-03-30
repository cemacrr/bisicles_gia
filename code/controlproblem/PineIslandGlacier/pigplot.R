require(hdf5, lib.loc="/local/R/lib")
#require(hdf5)
source("read.chombo.R")

nx <- 64
ny <- 96
xf <- seq(0,256,l=nx+1)
yf <- seq(0,256,l=ny+1)
xc <- 0.5*(xf[1:nx]+xf[(1:nx)+1])
yc <- 0.5*(yf[1:ny]+yf[(1:ny)+1])

gradsq <- function(z)
  {
    n <- dim(z)

    dx <- (xf[2]-xf[1])*1e+3
    dy <- (yf[2]-yf[1])*1e+3

    maskx <- ifelse ( z[2:n[1],] > 0 & z[1:(n[1]-1),] > 0 , 1, 0)
    masky <- ifelse ( z[,2:n[2]] > 0 & z[,1:(n[2]-1)] > 0 , 1, 0)
    gx <- maskx * (z[2:n[1],] - z[1:(n[1]-1),])/dx
    gy <- masky * (z[,2:n[2]] - z[,1:(n[2]-1)])/dy
    sum ( gx^2) + sum(gy^2)
    
  }


L2norm <- function(x)
  {
    sqrt ( sum ( x^2) /length(x) ) 
  }

maxnorm <- function(x)
  {
    max(abs(x))
  }


#pdf("~/Desktop/ctrl.pig.a.pdf",width=8,height=12)

sq  <-  c(11)
for (nfmax in sq)
  {
    #cL2err <- NULL
    vxL2err <- NULL
    vyL2err <- NULL
    #cmaxerr <- NULL
    vxmaxerr <- NULL
    divL2err <- NULL
    gs <- NULL
    nf <- seq(1,nfmax,by=1)
    for (iter in nf)
      {
        p <- read.chombo(paste("out",formatC(iter,flag="0",width=6),".2d.hdf5",sep=""),0)

        C <- p$level.0$boxes[[1]]$v3
        velx <- p$level.0$boxes[[1]]$v5
        vely <- p$level.0$boxes[[1]]$v6
        velox <- p$level.0$boxes[[1]]$v7
        veloy <- p$level.0$boxes[[1]]$v8 
        divuh <-  p$level.0$boxes[[1]]$v8
        divuho <-  p$level.0$boxes[[1]]$v10
        
        mask <- ifelse(C * abs(velox) > 0,1,0)
        #cL2err <- c(cL2err,L2norm(cpc-cbtrc))
        vyL2err <- c(vyL2err,L2norm( mask*(vely-veloy)))
        vxL2err <- c(vxL2err,L2norm( mask*(velx-velox)))
        divL2err <- c(divL2err,L2norm(divuh-divuho))
        
        gs <- c(gs,gradsq(C))
      }
    
    modvelobs <- sqrt(velox^2 + veloy^2)
    modvelcomp <- sqrt(velx^2 + vely^2)

    nx <- dim(velox)[1]
    ny <- dim(velox)[2]
    xf <- seq(0,256,l=nx+1)
    yf <- seq(0,384,l=ny+1) 

    xc <- 0.5 * ( xf[1:(nx)] + xf[2:(nx+1)])
    yc <- 0.5 * ( yf[1:(ny)] + yf[2:(ny+1)])
    par(mfrow=c(3,2))


    plot(range(nf),c(1e-2,2),type='n',log='y',
         xlab="function/gradient evaluations",
         ylab="error norms",
         main=paste("Error decay")
         )
    lines(nf,vxL2err/vyL2err[1],type='l',col="blue")
    lines(nf,vyL2err/vyL2err[1],type='l',col="orange")
    lines(nf,gs/gs[1],type='l',col="pink")
    #lines(nf,vxmaxerr/max(vxmaxerr),type='l',col="purple",lty=2,lwd=2)
    

    lines(nf,divL2err/divL2err[1],type='l',col="green")
 
    ucol <- colorRampPalette(c("white","blue","yellow"),space="Lab")(64)
    zl <- c(0,5000)

    pc <- ifelse(1 + C > 1+max(zl),1+max(zl),1 + C)

    ccol <- colorRampPalette(c("white","yellow","blue"),space="Lab")(64)
    
    image(xf,yf,log(pc),col=ccol,zlim=log(1+zl),
          #ylim=c(56,156),xlim=c(56,156),
          main="basal friction coefficient C(x,y)",
          sub=paste("after ",nfmax,"function/gradient evaluations"),
          xlab="x (km)", ylab = "y (km)")
    contour(xc,yc,C,add=TRUE,lev=c(0,1,2,4,8))

    zl <- c(1,4000)

    #vf <- function(a){ifelse(a > max(zl),max(zl),a)}
    vf <- function(a){log(a+1,10)}
    
    image(xf,yf,vf(modvelobs),col=ucol,zlim=vf(zl))
    contour(xc,yc,C,add=TRUE,lev=c(.1))
    image(xf,yf,vf(modvelcomp),col=ucol,zlim=vf(zl))
    contour(xc,yc,C,add=TRUE,lev=c(.1))
    
    verrlev <- seq(0,1000,by=200)
    errzl <- range(verrlev)
    errv <- mask*abs(modvelcomp-modvelobs)
    image(xf,yf,errv,col=ucol,zlim=errzl)
    contour(xc,yc,errv,add=TRUE,lev=verrlev,zlim=errzl)

    dcol <- colorRampPalette(c("red","white","blue"))(64)
    image(xf,yf,divuh,col=dcol,zlim=c(-150,150))
    contour(xc,yc,mask*divuh,lev=c(-100,-10,10,100),add=TRUE)
    
  }

#dev.off()
