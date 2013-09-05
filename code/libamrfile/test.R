require(libamrfile)


cg <- colorRampPalette(c("#003366","#54A9FF","#FFFF00","#FF0000"),space="rgb")(128)


uplot <- function(fab){
  #test plot : color plot of log(|u|,10)
  print(names(fab))
  modu <-sqrt(fab$v[,,1]^2 + fab$v[,,2]^2)
  modu <- ifelse(modu > 5000,5000,modu)
  image(fab$x,fab$y,log(modu,10),add=TRUE,zlim=log(c(1,5000),10),col=cg,useRaster=TRUE)
  
}


par(xaxs="i",yaxs="i",las=1,mar=c(1,1,1,1))
plot(c(0,768e+3),c(0,768e+3),type='n',axes=FALSE)
box()
pthwID <- amr.load("plot.amundsen.2d.hdf5")


names <- amr.query.compnames(pthwID)
ncomp <- length(names)
comp <- 0:(ncomp-1)
uxc <- comp[names=="xVel"]
uyc <- comp[names=="yVel"]

amr.apply(pthwID, uplot, comp=c(uxc,uyc),maxlev=2)
amr.free(pthwID)


