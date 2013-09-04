require(libamrfile)


cg <- colorRampPalette(c("#003366","#54A9FF","#FFFF00","#FF0000"),space="rgb")(128)


uplot <- function(fab){
  #test plot : color plot of log(|u|,10)
  print(names(fab))
  modu <-sqrt(fab$v^2 + fab$v2^2)
  modu <- ifelse(modu > 5000,5000,modu)
  image(fab$x,fab$y,log(modu,10),add=TRUE,zlim=log(c(1,5000),10),col=cg,useRaster=TRUE)
  
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


par(xaxs="i",yaxs="i",las=1,mar=c(1,1,1,1))
plot(c(0,768e+3),c(0,768e+3),type='n',axes=FALSE)
box()
pthwID <- amr.load("plot.amundsen.2d.hdf5")
umodify(pthwID)
amr.apply(pthwID, umodf, comp=c(1,2))
amr.free(pthwID)


