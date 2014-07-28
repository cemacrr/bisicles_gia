
amr.load <- function(f)
  {

    #test the file is present first...
    if (is.na(file.info(f)$size)){
      stop(c("file ",f," not found"))
      rc <- -1
    } else {
      r <- .C("amr_read_file_R",
              status=integer(1),
              amrID=integer(1),
              file=as.character(f))
      
      if (r$status == 0)
        {
          rc <- r$amrID
        }
      else
        {
          stop(c("file ",f," not opened"))
          rc <- -1
        }
    }
    rc
  }

amr.create <- function(nx, ny, dx, ncomp, nghost)
  {
    rc <- -1
    
    r <- .C("amr_create_coarse",
            status=integer(1),
            amrID=integer(1),
            nx=as.integer(nx),
            ny=as.integer(ny),
            dx=as.double(dx),
            ncomp=as.integer(ncomp),
            nghost=as.integer(nghost))
    if (r$status == 0)
      {
        rc <- r$amrID
      }
    else
      {
        stop(c("amr hierarchy not created"))
        rc <- -1
      }
    rc
    
  }

amr.write <- function(amrID, f)
  {
    r <- .C("amr_write_file_R",
            status=integer(1),
            amrID=as.integer(amrID),
            file=as.character(f))
    
    r$status
  }

amr.free <- function(amrID)
  {
    r <- .C("amr_free",
            status=integer(1),
            amrID=as.integer(amrID))
  }

amr.free.all <- function(amrID)
  {
    r <- .C("amr_free_all")
  }
amr.query.nlevel <- function(amrID)
  {
     r <- .C("amr_query_n_level",
            status=integer(1),
            nlevel=integer(1),
            amrID=as.integer(amrID))

     if (r$status == 0){
       r$nlevel
     } else {
       -1
     }
  }

amr.query.ncomp <- function(amrID, level)
{
  r <- .C("amr_query_n_comp",
          status=integer(1),
          ncomp=integer(1),
          amrID=as.integer(amrID))
   if (r$status == 0){
     r$ncomp
   } else {
     -1
   }

}


amr.query.compnames <- function(amrID)
{
  n <- amr.query.ncomp(amrID)
  t <- "spamspamspamspamspamspamspamspam"
  blen <- nchar(t)
  compnames <- NULL
  
  for (i in 0:(n-1))
    {
  
      r <- .C("amr_query_comp_name_R",
              status=integer(1),
              buf=as.character(t),
              amrID=as.integer(amrID),
              comp=as.integer(i),
              buflen=as.integer(blen))
      if (r$status == 0)
        {
          compnames <- c(compnames,r$buf)
        }
      else
        {
          compnames <- c(compnames,"error in name")
        }
    }

  compnames
}


amr.set.compnames <- function(amrID, names)
  {
    n <- amr.query.ncomp(amrID)
    if (n != length(names))
      {
        stop ("length(names) != number of components")
      }

    rc <- 0
    
    for (i in 1:n)
      {
      
        r <- .C("amr_set_comp_name_R",
                status=integer(1),
                name=as.character(names[i]),
                amrID=as.integer(amrID),
                comp=as.integer(i-1))
        rc <- max(rc,r$status)
      }
    rc
  }

amr.query.nfab <- function(amrID, level)
  {
    
    r <- .C("amr_query_n_fab",
            status=integer(1),
            nfab=integer(1),
            amrID=as.integer(amrID),
            level=as.integer(level))

     if (r$status == 0){
       r$nfab
     } else {
       -1
     }
    
}

amr.read.fab <- function(amrID, level, fab, comp, ng=0)
  {
    r <- .C("amr_query_fab_dimensions_2d",
           status=integer(1),nx=integer(1),ny=integer(1),
           ncomp=integer(1),amrID=as.integer(amrID),
           level=as.integer(level),fab=as.integer(fab))

    if (r$status == 0) {
  
      s <- .C("amr_read_fab_data_2d",
              status=integer(1),
              v=matrix(0,r$nx+2*ng,r$ny+2*ng),
              x=numeric(r$nx+2*ng),
              y=numeric(r$ny+2*ng),
              amrID=as.integer(amrID),
              level=as.integer(level),
              fab=as.integer(fab),
              comp=as.integer(comp),
              nghost=as.integer(ng))
      if (s$status == 0){
        s
      } else {
        -2
      }
    } else {
      -1
    }
  }


amr.read.box <- function(amrID, level, lo, hi, comp)
  {
  
    nx <- as.integer(hi[1]-lo[1]+1)
    ny <- as.integer(hi[2]-lo[2]+1)
    
    s <- .C("amr_read_box_data_2d",
            status=integer(1),
            v=matrix(0,nx,ny),
            x=numeric(nx),
            y=numeric(ny),
            amrID=as.integer(amrID),
            level=as.integer(level),
            lo=as.integer(lo),
            hi=as.integer(hi),
            comp=as.integer(comp))
    if (s$status == 0){
      s
    } else {
      -1
    }
  }



amr.write.fab <- function(amrID, level, fab, fabdata, comp, ng=0)
  {
    s <- .C("amr_write_fab_data_2d",
            status=integer(1),
            v=as.double(fabdata),
            nx=as.integer(dim(fabdata)[1]-2*ng),
            ny=as.integer(dim(fabdata)[2]-2*ng),
            amrID=as.integer(amrID),
            level=as.integer(level),
            fab=as.integer(fab),
            comp=as.integer(comp),
            nghost=as.integer(ng))
    
    if (s$status == 0){
      s
    } else {
      -2
    }
  }


amr.apply <- function(amrID, f, minlev = 0, maxlev = -1, comp = 0, nghost = 1, ...) 
{
#apply f(fab,...) to every fab in the amr hierarchy
  if (maxlev < 0){
    maxlev <- amr.query.nlevel(amrID) - 1
  }

  if (maxlev > amr.query.nlevel(amrID) - 1){
    stop("maxlev > amr.query.nlevel(amrID) - 1)")
  }
  
  if (minlev < 0){
    stop("minlev < 0")
  }

  for (lev in 0:maxlev){
    nfab <- amr.query.nfab(amrID,lev) - 1
    for (ifab in 0:nfab)
      {
        fab <- amr.read.fab(amrID,lev,ifab,comp[1],ng=nghost)
        nx <- length(fab$x)
        ny <- length(fab$y)
        if (length(comp) > 1)
          {
            for (icomp in comp[-1])
              {
                t <- amr.read.fab(amrID,lev,ifab,icomp,ng=nghost)
                fab$v <- cbind(fab$v, t$v)
              }
          }
        dim(fab$v) <- c(nx,ny,length(comp))
        fab$amrID <- amrID
        fab$level <- lev
        fab$ID <-  ifab
        fab$comp <- comp
        fab$nghost <- nghost
        r <- f(fab,...)

      }
  }
}
