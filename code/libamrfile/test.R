source("R/amrfile.R")


ub <- c(1e-3,5000)
ub <- c(1,5000)
uf <- function(u)
  {
    u <- ifelse(u > ub[2],ub[2],u)
    u <- ifelse(u < ub[1],ub[1],u)
    log(u,10)
  }

ubound <- log(ub,10)


cg <- colorRampPalette(c("#003366","#54A9FF","#FFFF00","#FF0000"),space="rgb")(128)

uplot <- function(amrID){
  #add a color map of the velocity field and a contour marking
  #the grounding line
  maxlev <- amr.query.nlevel(amrID) - 1
  for (lev in 0:maxlev){
    maxfab <- amr.query.nfab(amrID,lev) - 1
    for (i in 0:maxfab)
      {
        u <- amr.read.fab(amrID,lev,i,1,ng=1)
        v <- amr.read.fab(amrID,lev,i,2,ng=1)
        modu <- sqrt(u$v^2 + v$v^2)
        s <-  amr.read.fab(amrID,lev,i,4,ng=1)
        thk <-  amr.read.fab(amrID,lev,i,0,ng=1)
        sab <- s$v - (1-918/1028)*thk$v
        x <- u$x
        y <- u$y 
        modu <- ifelse(thk$v > 0.0, modu, NA)
        image(x,y,uf(modu),add=TRUE,zlim=ubound,col=cg,useRaster=TRUE)
        contour(x,y,sab,levels=c(1),add=TRUE,drawlabels=FALSE,lwd=1,col="black")
       
      }
  }
}

umodify <- function(amrID){
  #test modiication : double the x-velocity
  maxlev <- amr.query.nlevel(amrID) - 1
  for (lev in 0:maxlev){
    maxfab <- amr.query.nfab(amrID,lev) - 1
    for (i in 0:maxfab)
      {
        ucomp = 1
        u <- amr.read.fab(amrID,lev,i,ucomp,ng=1)
        u$v <- u$v * 2.0
        s <- amr.write.fab(amrID,lev,i,u$v,ucomp,ng=1)
      }
  }
}



#speed scale
uscale <- function()
  {
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
    
    #sq <- c(1,26,51,76,101)
    sq <- seq(1,length(y),l=6)
    text(xo-a/10,y[sq],paste(c("<",rep("",4),">"),format(10^u[sq],digits=1)),pos=2)
    segments(xo,y[sq],xo+xl,y[sq])
    
    text(xo-a/4,yo+yl+5*b,expression(abs(u) * phantom(0)* ( m*a^{-1} )))
  }



par(xaxs="i",yaxs="i",las=1,mar=c(1,1,1,1))
plot(c(0,768e+3),c(0,768e+3),type='n',axes=FALSE)
box()
pthwID <- amr.load("plot.amundsen.2d.hdf5")
#uplot(pthwID)

umodify(pthwID)


amr.write(pthwID, "plot.rmod.amundsen.2d.hdf5")

uplot(pthwID)

amr.free(pthwID)


