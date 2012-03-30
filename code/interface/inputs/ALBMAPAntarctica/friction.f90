module ncdump
 
 use netcdf
 !use mgrelax
 implicit none
  
 contains
 subroutine nccheck(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop 2
    end if
  end subroutine nccheck

  subroutine ncload(x,y,thck,topg,umod,stemp,n,m)
  
    integer n,m,i,j,nc_id,var_id
    real(kind=8), dimension(1:n) :: x
    real(kind=8), dimension(1:m) :: y
    real(kind=8), dimension(1:n,1:m,1) :: thck, topg, umod, stemp

    call nccheck( nf90_open("ALBMAP_i2s_4BISICLES.nc", NF90_NOWRITE, nc_id) )

    !x
    !call nccheck( nf90_inq_varid(nc_id, "x", var_id) )
    !call nccheck( nf90_get_var(nc_id, var_id , x) )
    !y
    !call nccheck( nf90_inq_varid(nc_id, "y", var_id) )
    !call nccheck( nf90_get_var(nc_id, var_id , y) )

    x(1) = 2500.0d0
    do i = 2,n
       x(i)=x(i-1)+5000.0d0
    end do

    y(1) = 2500.0d0
    do j = 2,m
       y(j)=y(j-1)+5000.0d0
    end do

    ! thickness
    call nccheck( nf90_inq_varid(nc_id, "thick", var_id) )
    call nccheck( nf90_get_var(nc_id, var_id , thck) )
    
    !topgraphy
    call nccheck( nf90_inq_varid(nc_id, "topg", var_id) )
    call nccheck( nf90_get_var(nc_id, var_id , topg) )

    !speed
    call nccheck( nf90_inq_varid(nc_id, "speed", var_id) )
    call nccheck( nf90_get_var(nc_id, var_id , umod) )

    !surface temperature
    call nccheck( nf90_inq_varid(nc_id, "temp", var_id) )
    call nccheck( nf90_get_var(nc_id, var_id , stemp) )

    call nccheck( nf90_close(nc_id) )
    return
  end subroutine ncload

  subroutine ncsave(x,y,thck,topg,umod,beta,usrf,lsrf,stemp,mask,n,m)

    integer n,m, nc_id, time_id, x0_id, y0_id, x1_id, y1_id, &
         thck_id, topg_id, beta_id, umod_id, usrf_id, mask_id, lsrf_id, &
         stemp_id
    integer node_dim_id(3),cell_dim_id(3)
    real(kind=8), dimension(1:n) :: x
    real(kind=8), dimension(1:m) :: y
    real(kind=8), dimension(1:(n+1)) :: xn
    real(kind=8), dimension(1:(m+1)) :: yn
    real(kind=8), dimension(1:n,1:m,1) :: thck, topg, umod, beta, usrf, lsrf, stemp, mask
    real(kind=8), dimension(1:(n+1),1:(m+1),1) :: thckn, topgn
    real(kind=8) :: time

    call nccheck( nf90_create("Antarctica-5km.nc", NF90_CLOBBER, nc_id) )

    !node centered dimensions
    call nccheck( nf90_def_dim(nc_id, "time", 1, node_dim_id(3)) )
    call nccheck( nf90_def_dim(nc_id, "x1", n+1, node_dim_id(1)) )
    call nccheck( nf90_def_dim(nc_id, "y1", m+1, node_dim_id(2)) )

    !cell centred dimensions
    cell_dim_id(3) = node_dim_id(3)
    call nccheck( nf90_def_dim(nc_id, "x0", n, cell_dim_id(1)) )
    call nccheck( nf90_def_dim(nc_id, "y0", m, cell_dim_id(2)) )

    !dimension array defintions
    call nccheck( nf90_def_var(nc_id, "time", nf90_real8, node_dim_id(3), time_id))
    call nccheck( nf90_put_att(nc_id, time_id, "units", "a") )
    call nccheck( nf90_def_var(nc_id, "x0", nf90_real8, cell_dim_id(1), x0_id))
    call nccheck( nf90_put_att(nc_id, x0_id, "units", "m") )
    call nccheck( nf90_def_var(nc_id, "y0", nf90_real8, cell_dim_id(2), y0_id))
    call nccheck( nf90_put_att(nc_id, y0_id, "units", "m") )
    call nccheck( nf90_def_var(nc_id, "x1", nf90_real8, node_dim_id(1), x1_id))
    call nccheck( nf90_put_att(nc_id, x1_id, "units", "m") )
    call nccheck( nf90_def_var(nc_id, "y1", nf90_real8, node_dim_id(2), y1_id))
    call nccheck( nf90_put_att(nc_id, y1_id, "units", "m") )



    !cell centered array definitions
    call nccheck( nf90_def_var(nc_id, "thk", nf90_real8, cell_dim_id, thck_id) )
    call nccheck( nf90_def_var(nc_id, "topg", nf90_real8, cell_dim_id, topg_id) )
    call nccheck( nf90_def_var(nc_id, "umod", nf90_real8, cell_dim_id, umod_id) )
    call nccheck( nf90_def_var(nc_id, "beta", nf90_real8, cell_dim_id, beta_id) )
    call nccheck( nf90_def_var(nc_id, "usrf", nf90_real8, cell_dim_id, usrf_id) )
    call nccheck( nf90_def_var(nc_id, "lsrf", nf90_real8, cell_dim_id, lsrf_id) )
    call nccheck( nf90_def_var(nc_id, "stemp", nf90_real8, cell_dim_id, stemp_id) )
    call nccheck( nf90_def_var(nc_id, "mask", nf90_real8, cell_dim_id, mask_id) )
    !node centered array definitions


    call nccheck( nf90_enddef(nc_id) )
    
    !dimension array data
    time = 0.0d0
    call nccheck( nf90_put_var(nc_id, time_id, time) )
    call nccheck( nf90_put_var(nc_id, x0_id, x) )
    call nccheck( nf90_put_var(nc_id, y0_id, y) )
    !node centred x and y
    xn(1:n) = x(1:n) - (x(2)-x(1))
    xn(n+1) = x(n) + (x(2)-x(1)) 
    call nccheck( nf90_put_var(nc_id, x1_id, xn) )
    yn(1:m) = y(1:m) - (y(2)-y(1))
    yn(m+1) = y(m) + (y(2)-y(1)) 
    call nccheck( nf90_put_var(nc_id, y1_id, yn) )

    !cell centered array data 
    call nccheck( nf90_put_var(nc_id, thck_id, thck) )
    call nccheck( nf90_put_var(nc_id, topg_id, topg) )
    call nccheck( nf90_put_var(nc_id, umod_id, umod) )
    call nccheck( nf90_put_var(nc_id, beta_id, beta) )
    call nccheck( nf90_put_var(nc_id, usrf_id, usrf) )
    call nccheck( nf90_put_var(nc_id, lsrf_id, lsrf) )
    call nccheck( nf90_put_var(nc_id, mask_id, mask) )
    call nccheck( nf90_put_var(nc_id, stemp_id, stemp) )
    !node centered array data

    call nccheck( nf90_close(nc_id) )
    return
  end subroutine ncsave
end module ncdump

subroutine computeA(flwa,stemp,thck,n,ewn,nsn)
  implicit none
  integer :: ewn,nsn
!  real(kind=8), parameter :: Cgas = 0.16612d0, Q = 7.88d4, R=8.31d0, K=1.17d0, B0=2.207, &
!       tR = 273.39d0, enhance=1.0d0
  real(kind=8) :: Cgas, Q, R, K , B0, tR, enhance, rhog, icepmeltf
  real(kind=8), dimension(1:ewn,1:nsn) :: flwa, stemp, thck
  real(kind=8) :: n, factor, theta, tmp, pmp, pfactor
  integer i, j
  
  Cgas = 0.16612d0
  Q = 7.88d4
  R=8.31d0
  K=1.17d0
  B0=2.207
  tR = 273.39d0
  enhance=1.0d0
  factor = enhance* ( (1.0d0/B0)**n )
  pfactor = 9.81d0 * 918d0 * 9.7456d-8 
  do i = 1,ewn
     do j = 1,nsn
        
        theta = (stemp(i,j) + 273.16) * (1.0d0 + 0.01 * thck(i,j) * pfactor)
        pmp = 273.0
        theta = (min(pmp,theta))

        tmp = 3.0*Cgas/(tR - theta)**K - Q/R / theta
        flwa(i,j) = factor * exp(tmp)
     end do
  end do
  return
end subroutine computeA

program t
  use ncdump
  use mgrelax
  implicit none
  integer, parameter :: ewn = 1280, nsn = 1280
  real (kind=8), parameter :: rhoi = 918.0d0, rhoo = 1028.0d0, grav = 9.81d0
  real (kind=8), parameter :: cm = 1.0d-2, glen_n = 3.0,  lambda = 8.0e+3
  real (kind=8), parameter :: minbeta=20, maxbeta=1.0e+6
  real (kind=8), parameter :: depth_thresh = 120.d0, thck_thresh = 120.0d0
  integer, parameter :: typ_grounded = 0, typ_iceshelf = 1, typ_opensea = 2


  real(kind=8), dimension(1:ewn) :: x
  real(kind=8), dimension(1:nsn) :: y
  real(kind=8), dimension(1:ewn,1:nsn) :: thck, topg, umod, beta, usrf, lsrf, &
       dsx, dsy, umodsia, mask, stemp, flwa

  real(kind=8), dimension(0:ewn+1,0:nsn+1) :: umg, vmg, dumg, rmg
  integer, dimension(1:ewn,1:nsn) :: typ
  real(kind=8) :: time, dx, mu, resNorm
  integer iter, niter, i, j

  call ncload(x,y,thck,topg,umod,stemp,ewn,nsn)

  dx = x(2)-x(1)

  !decide usrf from flotation criterion
  usrf = max(topg + thck, thck * (1.0 - rhoi/rhoo))
  lsrf = usrf - thck

  
  typ = typ_grounded ! grounded ice
  where (lsrf - topg .gt. cm)
     typ = typ_iceshelf ! ice shelf
  end where

  where (topg(1024:1280,1:nsn).eq.0.0d0)
     topg(1024:1280,1:nsn) = -1.0d+4
  end where
  where (topg(1:ewn,1024:1280).eq.0.0d0)
     topg(1:ewn,1024:1280) = - 1.0d+4
  end where

  where ((lsrf - topg .gt. depth_thresh) .and. (thck .lt. thck_thresh))
     typ = typ_opensea ! open sea
     thck = 0.0
     usrf = 0.0
     lsrf = 0.0
  end where

  !compute surface gradients
  dsx = 0
  dsy = 0
  
  where ((typ (2:ewn-1,1:nsn) .eq. typ_grounded) &
       .and. (typ (3:ewn,1:nsn) .eq. typ_grounded) &
       .and. (typ (1:ewn-2,1:nsn) .eq. typ_grounded))
     dsx(2:ewn-1,1:nsn) = (usrf(3:ewn, 1:nsn) - usrf(1:ewn-2, 1:nsn)) &
          /(2.0d0 * dx)
  end where


  where ((typ (1:ewn,2:nsn-1) .eq. typ_grounded) &
       .and. (typ (1:ewn,3:nsn) .eq. typ_grounded) &
       .and. (typ (1:ewn,1:nsn-2) .eq. typ_grounded))
     dsy(1:ewn,2:nsn-1) = (usrf(1:ewn,3:nsn) - usrf(1:ewn,1:nsn-2)) &
          /(2.0d0 * dx)
  end where

  !A in Glen's law
  call computeA(flwa,stemp,thck,glen_n,ewn,nsn)
  
  !basal traction coefficient
  where (typ.eq.typ_grounded)
     umodsia = (2.0d0*flwa *thck**(glen_n+1)) / (glen_n+1.0d0) * (rhoi* grav)**glen_n  & 
          * sqrt(dsx**2 + dsy**2)**glen_n
     beta = (1 + rhoi* grav * thck * sqrt(dsx**2 + dsy**2)) / ( max(1.0d-6,umod - umodsia))
  end where

  beta = min(beta,maxbeta)
  beta = max(beta,minbeta)

  !smooth beta
  if (1 .eq. 1) then
     !smooth beta by solving umg -  lambda**2 * div(grad(umg) = vmg
     mu = 1.0*(lambda**2) / (dx**2)
     rmg  = 0.0d0
     vmg = log(minbeta)
     vmg(1:ewn,1:nsn) = log(beta)
     !vmg = beta
     umg = vmg

     call resid(umg,vmg,rmg,mu,ewn,nsn)
     resNorm = sum(abs(rmg) )
     write(*,*) 'initial res= ', resNorm
     niter = 10
     do iter = 1, niter
        dumg = 0.0
        call vcycle(dumg,rmg,mu,ewn,nsn,10,4)
        umg = umg + dumg
        call resid(umg,vmg,rmg,mu,ewn,nsn)
        resNorm = sum(abs(rmg) )
        write (*,*) 'MG iter', iter, " norm = ", resNorm
    
     end do
     beta = exp(umg(1:ewn,1:nsn))
     beta = min(beta,maxbeta)
     beta = max(beta,minbeta)
  end if

  mask = typ
  call ncsave(x,y,thck,topg,umod,beta,usrf,lsrf,flwa,mask,ewn,nsn)
 
end program t
