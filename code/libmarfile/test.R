source("R/amrfile.R")


ub <- c(1e-3,5000)

uf <- function(u)
  {
    u <- ifelse(u > ub[2],ub[2],u)
    u <- ifelse(u < ub[1],ub[1],u)
    log(u,10)
  }

ubound <- log(ub,10)
s <- seq(0,1,l=5)
lc <- rgb(s,0,0,0.15)
#cg <- colorRampPalette(c("#000000","#3333FF","#EEEEFF"),space="rgb")(128)
#cg <- colorRampPalette(c("blue","white","red"),space="rgb")(128)
#cg <- c("#000000",colorRampPalette(c("blue","white","orange"),space="rgb")(128))
#cg <- c("#000000",colorRampPalette(c("#54A9FF","white","orange"),space="rgb")(128))
cg <- colorRampPalette(c("#003366","#54A9FF","#EEEEFF"),space="rgb")(128)
#cg <- colorRampPalette(c("#003366","#54A9FF","#EEEEFF","#FFFF00","#FF0000"),space="rgb")(128)

maskcol <- rgb(c(0.5,1), c(0.5,1), c(0.5,1), c(0.2,0.0))

uplot <- function(amrID, x0, y0, mag = 1.0, maxlev = 100, mesh=FALSE){
  #add a color map of the velocity field and a contour marking
  #the grounding line, to a plot, with the bottom left at x0, y0
  maxlev <- min(amr.query.nlevel(amrID) - 1,maxlev)
  for (lev in 0:maxlev){
    maxfab <- amr.query.nfab(amrID,lev) - 1
    for (i in 0:maxfab)
      {
        u <- amr.read.fab(amrID,lev,i,1)
        v <- amr.read.fab(amrID,lev,i,2)
        modu <- sqrt(u$v^2 + v$v^2)
        bf <-  amr.read.fab(amrID,lev,i,7)
        s <-  amr.read.fab(amrID,lev,i,4)
        thk <-  amr.read.fab(amrID,lev,i,0)
        sab <- s$v - (1-918/1028)*thk$v
        x <- mag*u$x + x0
        y <- mag*u$y + y0
        image(x,y,uf(modu),add=TRUE,zlim=ubound,col=cg,useRaster=TRUE)
                                        #image(x,y,bf$v,add=TRUE,zlim=c(0,1),useRaster=TRUE,col=maskcol)
        #contour(x,y,sab,levels=c(1),add=TRUE,drawlabels=FALSE,lwd=3)
        if (mesh & lev > 0)
          {
            dx <- x[2]-x[1]
            xm <- min(x)+dx/2
            xp <- max(x)-dx/2
            dy <- y[2]-y[1]
            ym <- min(y)+dy/2
            yp <- max(y)-dy/2
            rect(xm,ym,xp,yp,border=lc[lev+1])
          }
        contour(x,y,sab,levels=c(1),add=TRUE,drawlabels=FALSE,lwd=1,col="black")
      }
  }
}

par(xaxs="i",yaxs="i",las=1,mar=c(5,5,1,1))

waismap <- function(risfris,pigthw)
  {

par(xaxs="i",yaxs="i",las=1,mar=c(1,1,1,1))
plot(c(0,5*640e+3),c(0,5*640e+3),type='n',axes=FALSE,xlab="x (km)",ylab="y (km)")
box()
sc <- seq(0,5*640,by=1e3)
#axis(1,at=sc*1e3,lab=sc)
#axis(2,at=sc*1e3,lab=sc)


risfrisID <- amr.load(risfris)
uplot(risfrisID,0,0,mesh=FALSE)
amr.free(risfrisID)

#length scale on the main map
xos <- 150e+3
yos <- 150e+3
xl <- 1000e+3
yl <- 40e+3
xt <- c(200,600)*1e3
dxt <- 0.5*(xt[2]-xt[1])
rect(xos,yos,xos+xl,yos+yl,col="yellow")
rect(xos+xt,yos,xos+xt+dxt,yos+yl,col="black")
text(xos+xl+dxt,yos+yl/2,"1000 km",col="yellow")



pthwID <- amr.load(pigthw)

#X2.5 inset
mag <- 2.5
xo <- 5*640e+3 - mag*(512e+3)
yo <- 5*640e+3 - mag*(768e+3)
xm <- -10e3 + xo; xp <- xm +  mag*(512 + 20) * 1e+3
ym <- -10e3 + yo; yp <- ym +  mag* (768 + 20) * 1e+3
polygon(c(xm,xm,xp,xp),c(ym,yp,yp,ym),col="#54A9FF",lwd=2)
uplot(pthwID,xo,yo,mag=mag,mesh=FALSE)


#length scale on X2.5  inset
xos <- mag*150e+3 + xo
yos <- mag*720e+3 + yo
xl <- mag*200e+3 
yl <- 40e+3
xt <- c(50,150)*1e3*mag
dxt <- 0.5*(xt[2]-xt[1])
rect(xos,yos,xos+xl,yos+yl,col="yellow")
rect(xos+xt,yos,xos+xt+dxt,yos+yl,col="black")
text(xos+xl+200e+3,yos+yl/2,"200 km",col="yellow")


#X1 inset
xo <- 90e+3
yo <- 900e+3
xm <-  xo -5e3; xp <- xm + (512 + 10) * 1e+3
ym <-  yo -5e3; yp <- ym + (768 + 10) * 1e+3
polygon(c(xm,xm,xp,xp),c(ym,yp,yp,ym),col="#54A9FF",lwd=2)
uplot(pthwID,xo,yo,mag=1.0)

amr.free(pthwID)

#speed scale
xo <- 5*640e+3 - 220e+3
yo <- 100e+3
yl <- 750e+3
xl <- 150e+3

x <- seq(xo,xo+xl,by=10e+3)
y <-  seq(yo,yo+yl,by=10e+3)
#u <- uf(seq(1e-2,ubound[2],l=length(y)))
u <- seq(uf(ub[1]),uf(ub[2]),l=length(y))

a<- 400e+3
b <- 40e+3
rect(xo-a,yo-2*b,xo+xl+b,yo+yl+7*b,col="white")

image(x,y,matrix(rep(u,each=length(x)),length(x),length(y)),
      col=cg,add=TRUE,useRaster=TRUE)

sq <- c(1,26,51,76,101)

text(xo-a/10,y[sq],format(10^u[sq],digits=2),pos=2)
segments(xo,y[sq],xo+xl,y[sq])

text(xo-a/4,yo+yl+5*b,expression(abs(u) * phantom(0)* ( m*a^{-1} )))


}


pdf("~/Desktop/test.pdf",width=8.0,height=8.0,pointsize=11)

waismap("plot.risfris.5km.l1l2.CONTROL_RUN_EXTRA_ACC.melt.2lev.000293.2d.hdf5",
        "plot.pigthwaites.1km.l1l2.t263.shelf-steady.000000.2d.hdf5")


waismap("plot.risfris.5km.l1l2.P16_A1B_FES_HA3-SMB_A1B_RAC_HA3_EXTRA_ACC.melt.2lev.002296.2d.hdf5",
        "plot.pigthwaites.1km.l1l2.t263.gllocalized-melt-soften.000738.2d.hdf5")

dev.off()
