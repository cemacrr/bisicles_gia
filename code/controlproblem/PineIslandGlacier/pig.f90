module ncdump
 
 use netcdf
 use mgrelax
 implicit none
  
 contains
 subroutine nccheck(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then
       print *, trim(nf90_strerror(status))
       stop 2
    end if
  end subroutine nccheck

  subroutine yflip(a,n,m)
    integer :: i,n,m
    real(kind=8), dimension(1:n,1:m) :: a,b
    
    do i = 1,m
       b(1:n,m-i+1) = a(1:n,i)
    end do
    a = b
  end subroutine yflip

  subroutine ncdump2D(x,y, thck, topg, lsrf, usrf, us, vs, ucoef, divuh, &
    divuhc, umod, beta, mask, n, m)
    
    integer n, m, nc_id, time_id,  x_id, y_id, &
         topg_id, lsrf_id,usrf_id, thck_id,  &
         us_id, vs_id, umod_id, beta_id, mask_id, ucoef_id, divuhc_id, divuh_id
    integer :: node_dim_ids(2), cell_dim_ids(2)
    real(kind=8), dimension(1:n) :: x
    real(kind=8), dimension(1:m) :: y
    !cell centred arrays
    real(kind=8), dimension(1:n,1:m) :: thck, topg, lsrf, usrf,  & 
         us , vs, umod, mask, beta, ucoef, divuh, divuhc

   

    call nccheck( nf90_create("pigv5.1km.cellcentered.nc", NF90_CLOBBER, nc_id) )

    !call nccheck( nf90_def_dim(nc_id, "time", 1, cell_dim_ids(3)) )
    call nccheck( nf90_def_dim(nc_id, "x", n, cell_dim_ids(1)) )
    call nccheck( nf90_def_dim(nc_id, "y", m, cell_dim_ids(2)) )


    !call nccheck( nf90_def_var(nc_id, "time", nf90_real8, cell_dim_ids(1), time_id))
    !call nccheck( nf90_put_att(nc_id, time_id, "units", "a") )

    call nccheck( nf90_def_var(nc_id, "x", nf90_real8, cell_dim_ids(1), x_id))
    call nccheck( nf90_def_var(nc_id, "y", nf90_real8, cell_dim_ids(2), y_id))
    call nccheck( nf90_put_att(nc_id, x_id, "units", "m") )
    call nccheck( nf90_put_att(nc_id, y_id, "units", "m") )

   
   
    call nccheck( nf90_def_var(nc_id, "thck", nf90_real8, cell_dim_ids, thck_id) )
    
    call nccheck( nf90_def_var(nc_id, "topg", nf90_real8, cell_dim_ids, topg_id) )
    
    call nccheck( nf90_def_var(nc_id, "lsrf", nf90_real8, cell_dim_ids, lsrf_id) )
    call nccheck( nf90_def_var(nc_id, "usrf", nf90_real8, cell_dim_ids, usrf_id) )
    call nccheck( nf90_def_var(nc_id, "velx", nf90_real8, cell_dim_ids, us_id) )
    call nccheck( nf90_def_var(nc_id, "vely", nf90_real8, cell_dim_ids, vs_id) )
    call nccheck( nf90_def_var(nc_id, "umod", nf90_real8, cell_dim_ids, umod_id) )
    call nccheck( nf90_def_var(nc_id, "btrc", nf90_real8, cell_dim_ids, beta_id) )
    call nccheck( nf90_def_var(nc_id, "mask", nf90_real8, cell_dim_ids, mask_id) )
    call nccheck( nf90_def_var(nc_id, "velc", nf90_real8, cell_dim_ids, ucoef_id) )
    call nccheck( nf90_def_var(nc_id, "divuh", nf90_real8, cell_dim_ids, divuh_id) )
    call nccheck( nf90_def_var(nc_id, "divuhc", nf90_real8, cell_dim_ids, divuhc_id) )
    call nccheck( nf90_enddef(nc_id) )
    
    !call nccheck( nf90_put_var(nc_id, time_id, time) )
    call nccheck( nf90_put_var(nc_id, x_id, x) )
    call nccheck( nf90_put_var(nc_id, y_id, y) )
 
    call yflip(thck,n,m)
    call nccheck( nf90_put_var(nc_id, thck_id, thck) )
    call yflip(topg,n,m)
    call nccheck( nf90_put_var(nc_id, topg_id, topg) )
    call yflip(lsrf,n,m)
    call nccheck( nf90_put_var(nc_id, lsrf_id, lsrf) )
    call yflip(usrf,n,m)
    call nccheck( nf90_put_var(nc_id, usrf_id, usrf) )
    call yflip(us,n,m)
    call nccheck( nf90_put_var(nc_id, us_id, us) )
    call yflip(vs,n,m)
    vs = -vs
    call nccheck( nf90_put_var(nc_id, vs_id, vs) )
    call yflip(umod,n,m)
    call nccheck( nf90_put_var(nc_id, umod_id, umod) )
    call yflip(beta,n,m)
    call nccheck( nf90_put_var(nc_id, beta_id, beta) )
    call yflip(mask,n,m)
    call nccheck( nf90_put_var(nc_id, mask_id, mask) )
    call yflip(ucoef,n,m)
    call nccheck( nf90_put_var(nc_id, ucoef_id, ucoef) )
    call yflip(divuh,n,m)
    call nccheck( nf90_put_var(nc_id, divuh_id, divuh) )
    call yflip(divuhc,n,m)
    call nccheck( nf90_put_var(nc_id, divuhc_id, divuhc) )

    call nccheck( nf90_close(nc_id) )
   
  end subroutine ncdump2D

end module ncdump

subroutine applydrop(drop,usrf,lsrf,topg,thck,ewn,nsn)
  implicit none
  integer, intent(in) :: ewn,nsn
  real(kind=8),parameter :: rhoi = 918.0d0, rhoo = 1028.0d0, grav = 9.81d0, cm = 1.0e-5
  real (kind = 8), dimension(1:ewn,1:nsn), intent(inout) :: topg, lsrf, usrf,  thck
  real (kind = 8), intent(in) :: drop
  
  usrf = usrf - drop
  topg = topg - drop
  lsrf = lsrf - drop
  ! some fiddling to ensure the division into sheet/shelf is
  ! as the data says, and usrf is correct also. topography may change
  ! near the grounding line.
  usrf = max(usrf,0.0)
  lsrf = max(lsrf,topg)
  thck = max(usrf - lsrf,0.0)
  where (lsrf .gt. topg + 1)
     thck = usrf / (1.0 - rhoi/rhoo) 
     lsrf = usrf - thck
     topg = min(topg,lsrf - cm)
  elsewhere
     thck = usrf - lsrf
     where (topg+thck.le.(1.0 - rhoi/rhoo)*thck)
        thck = - rhoo * topg / rhoi + cm
     end where
  end where
  !now use the bisicles expressions, 
  usrf = max(topg + thck, thck * (1.0 - rhoi/rhoo))
  lsrf = usrf - thck
 
end subroutine applydrop

subroutine testdrop(drop,res,usrf,lsrf,topg,thck,ewn,nsn)
  implicit none
  integer, intent(in) :: ewn,nsn
  real (kind = 8), intent(in) :: drop
  real (kind = 8), intent(out) :: res
  real (kind=8),parameter :: rhoi = 918.0d0, rhoo = 1028.0d0, grav = 9.81d0, cm = 1.0e-2
  real (kind = 8), dimension(1:ewn,1:nsn), intent(in) :: topg, lsrf, usrf,  thck
  
  real (kind = 8), dimension(1:ewn,1:nsn) :: ttopg, tlsrf, tusrf,  tthck
  tusrf = usrf
  tlsrf = lsrf
  ttopg = topg
  tthck = thck
  call applydrop(drop,tusrf,tlsrf,ttopg,tthck,ewn,nsn)
  
  res = sum(abs(ttopg - (topg-drop)), mask = ((lsrf.gt.topg) .and. (lsrf-topg.lt.5)))

end subroutine testdrop


program t
   use ncdump
  implicit none
  
  integer, parameter :: ewn = 292, nsn = 465, unin=9
  integer, parameter :: nx = 256, ny = 384
  integer, parameter :: wb = 35, sb = nsn - ny - 32, eb = wb + nx - 1, nb = sb + ny -1
  integer, parameter :: m = 5, niter = 9
  real(kind=8),parameter :: scyr = 31556926.0d0, maxsrf=1250.0d0 
  real(kind=8),parameter :: rhoi = 918.0d0, rhoo = 1028.0d0, grav = 9.81d0
  real(kind=8),parameter :: dx = 1.0d+3, dy = dx,  tol = 1.0d-3
  real(kind=8),parameter :: minbeta = 5.0, maxbeta = 1.0e+4, seabeta = 100.0, shelfbeta = 1.0e-5,cm = 1.0e-2
  real (kind = 8), dimension(1:ewn,1:nsn) :: topg, lsrf, usrf,  thck, us, vs, mask
  real (kind = 8), dimension(1:ewn,1:nsn) :: ttopg, tlsrf, tusrf,  tthck
  real (kind = 8), dimension(wb:eb,sb:nb) :: topgb, lsrfb, usrfb,  thckb, usb, vsb , maskb, ucoef,divuh,divuhc
  real (kind = 8), dimension(wb:eb,sb:nb) ::    umodb, betab,betar, dsx, dsy, umodsia
  real (kind = 8), dimension(wb-1:eb+1,sb-1:nb+1) :: umg,vmg,dumg,rmg
  integer , dimension(1:ewn,1:nsn) ::  typ
  real (kind = 8), dimension(1:ewn) :: x
  real (kind = 8), dimension(1:nsn) :: y
  integer, parameter :: nre = 4
  real (kind = 8), dimension(-nre:nre,-nre:nre) :: ra
  real (kind = 8) :: mu, lambda, flwa, glen_n, res0,res1,res,drop0,drop1,drop
  integer :: ew, ns, iter, radius, i,j,k

  open(unin,file='topg_5pg.dat',status='old')
  read(unin,*) ((topg(ew,ns), ew=1,ewn), ns=1,nsn)
  close(unin)

  open(unin,file='lsrf_5pg.dat',status='old')
  read(unin,*) ((lsrf(ew,ns), ew=1,ewn), ns=1,nsn)
  close(unin)
  
  open(unin,file='usrf_5pg.dat',status='old')
  read(unin,*) ((usrf(ew,ns), ew=1,ewn), ns=1,nsn)
  close(unin)

  open(unin,file='uvel_5pg.dat',status='old')
  read(unin,*) ((us(ew,ns), ew=1,ewn), ns=1,nsn)
  close(unin)

  open(unin,file='vvel_5pg.dat',status='old')
  read(unin,*) ((vs(ew,ns), ew=1,ewn), ns=1,nsn)
  close(unin)

  typ = 0 ! 0 for grounded ice

  x(1) = 0.0
  do ew = 2,ewn
     x(ew) = x(ew-1) + dx
  end do
  y(1) = 0.0
  do ns = 2,nsn
     y(ns) = y(ns-1) + dy
  end do
  
  !raw data - need to set up typ,mask,thck
  typ = 0 ! grounded ice
  thck = usrf - lsrf
  where (lsrf - topg .gt. cm)
     typ = 1 ! ice shelf
  end where

  where ((lsrf - topg .gt. 120.0) .and. (thck .lt. 120.0))
     typ = 2 ! open sea
     thck = 0.0
     usrf = 0.0
     lsrf = 0.0
  end where
  
  !double copy of typ
  mask = lsrf - topg

  thckb = thck(wb:eb,sb:nb)
  topgb = topg(wb:eb,sb:nb)


  

  lsrfb = lsrf(wb:eb,sb:nb)
 
  usrfb = usrf(wb:eb,sb:nb)
  maskb = mask(wb:eb,sb:nb)
  usb = us(wb:eb,sb:nb)
  vsb = -vs(wb:eb,sb:nb)
  
!!$  rhostar = rhoi
!!$  where (lsrf .gt.  topg .and. thck .gt. cm)
!!$     rhostar = (1.0-usrf/(thck))*rhoo
!!$  end where
!!$
!!$  write(*,*) sum (rhostar,mask =(lsrf .gt.  topg .and. thck .gt. cm) ) &
!!$       /   count(mask =(lsrf .gt.  topg .and. thck .gt. cm) ), &
!!$       maxval (rhostar,mask =(lsrf .gt.  topg .and. thck .gt. cm) ), &
!!$       minval (rhostar,mask =(lsrf .gt.  topg .and. thck .gt. cm) )

  drop =  15.2! Value from Anne Le Brocq

  call applydrop(drop,usrf,lsrf,topg,thck,ewn,nsn)
  
  !now use the bisicles expressions, 
  usrf = max(topg + thck, thck * (1.0 - rhoi/rhoo))
  lsrf = usrf - thck
 
  thckb = thck(wb:eb,sb:nb)
  topgb = topg(wb:eb,sb:nb)
  lsrfb = lsrf(wb:eb,sb:nb)
  usrfb = usrf(wb:eb,sb:nb)
  maskb = mask(wb:eb,sb:nb)
  usb = -us(wb:eb,sb:nb)
  vsb = -vs(wb:eb,sb:nb)


  !attempt to get rid of the missing values
  if (1.eq.0) then
     do ew = wb,eb
        do ns = sb,nb
           if ( abs(usb(ew,ns)) .gt. 9000.0) then
              usb(ew,ns) = usb(ew,ns-1)
              vsb(ew,ns) = vsb(ew,ns-1)
           end if
        end do
     end do
  end if

  !assume same 1/ (s * sigma) everywhere
  ucoef = 1.0d0;
  divuhc = 1.0d0;

  ! at the moment, we don't know what div(uh) should be, but
  ! we don't mind if it is large (and positive) in shelf regions 
  divuh = 0.0d0; 
  where (typ(wb:eb,sb:nb) .ne. 0)
      divuhc = 0.0d0
  end where

  where (usb(wb:eb,sb:nb) .gt. 9000.0)
     ucoef = 0.0d0
     usb = 0
     vsb = 0
  end where

  where (typ(wb:eb,sb:nb) .eq. 2)
     usb = 0.0
     vsb = 0.0
     ucoef = 0.0d0
  end where

  !set uceof = 0 close to open sea
  radius = 4
  do ew = wb+radius,eb-radius
     do  ns = sb+radius,nb-radius
        if (typ(ew,ns).eq.2) then
           do i = -radius,radius
              do j = -radius,radius
                 ucoef(ew+i,ns+j) = 0.0;
              end do
           end do
        end if
     end do
  end do


  !expand the region where ucoef = 0 by a few cells
  do iter = 1,5
     do ew = wb+1,eb-1
        do ns = sb+1,nb-1
           if (ucoef(ew,ns).gt.0.1) then
           ucoef(ew,ns) = 0.25 * ( &
                ucoef(ew-1,ns) + ucoef(ew+1,ns) &
                +ucoef(ew,ns-1) + ucoef(ew,ns+1))
           end if
        end do
     end do
!!$     do ns = sb+1,nb-1
!!$        do ew = wb+1,eb-1
!!$           if (ucoef(ew,ns).gt.0.1) then
!!$           ucoef(ew,ns) = 0.25 * ( &
!!$                ucoef(ew-1,ns) + ucoef(ew+1,ns) &
!!$                +ucoef(ew,ns-1) + ucoef(ew,ns+1))
!!$           end if
!!$        end do
!!$     end do
  end do
  where (ucoef.lt.0.95)
     ucoef = 0.0d0
  end where

  umodb = sqrt(usb**2 + vsb**2)

 !compute grad s
 dsx = 0
 dsy = 0
 
 where ((typ (wb+1:eb-1,sb:nb) .eq. 0) &
      .and. (typ (wb+2:eb,sb:nb) .eq. 0) &
      .and. (typ (wb:eb-2,sb:nb) .eq. 0))
    dsx(wb+1:eb-1,sb:nb) = (usrfb(wb+2:eb, sb:nb) - usrfb(wb:eb-2, sb:nb)) &
         /(2.0d0 * dx)
 end where

 where ((typ (wb:eb,sb+1:nb-1) .eq. 0) &
      .and. (typ (wb:eb,sb+2:nb) .eq. 0) &
      .and. (typ (wb:eb,sb:nb-2) .eq. 0))

    dsy(wb:eb,sb+1:nb-1) = (usrfb(wb:eb, sb+2:nb) - usrfb(wb:eb, sb:nb-2)) &
         /(2.0d0 * dy)

 end where

 !attempt ot compute a basal stress coeffcient
 flwa = 5.0d-17
 glen_n = 3.0
 umodsia = 0.0
 betar = 1.0e+6
 where (typ(wb:eb,sb:nb).eq.0)
    umodsia = (2.0d0*flwa *thckb**(glen_n+1)) / (glen_n+1.0d0) * (rhoi* grav)**glen_n  & 
         * (sqrt(dsx**2 + dsy**2))**glen_n
    betar = (1 + rhoi* grav * thckb * sqrt(dsx**2 + dsy**2)) / ( max(1.0d-4,umodb - umodsia))
 end where
 betar = min(betar,maxbeta)
 betar = max(betar,minbeta)
 !where (typ(wb:eb,sb:nb) .gt. 0)
 !   betar = minbeta
 !end where

 !fill in the basal traction where we don't know it
 if (1.eq.1) then
    do ns = sb+100,nb
       do ew = wb,eb
          if (typ(ew,ns) .gt. 0) then
             call random_number(ra)
             ra = ra + 1.0e-2 * ( 1.0e+2 + umodb(ew-nre:ew-nre,(ns-1-2*nre):(ns-1) ))
             ra = ra / sum(ra)
             betar(ew,ns) =  sum(ra * betar(ew-nre:ew-nre,(ns-1-2*nre):(ns-1) )) 
          end if
       end do
    end do
 else
    where (typ(wb:eb,sb:nb).gt.0)
       betar=shelfbeta
    end where
 end if

 where (ucoef.lt.1.0e-2)
    betab=shelfbeta
 end where

 if (1 .eq. 1) then
     !smooth beta by solving umg -  lambda**2 * div(grad(umg) = vmg,
    ! vmg = log(betar), u = log(betab)
    lambda = 4.0e+3
    mu = (lambda**2) / (dx**2)
    rmg  = 0.0d0
    vmg = log(maxbeta);
    vmg(wb:eb,sb:nb) = log(betar)
    umg = vmg

    call resid(umg,vmg,rmg,mu,nx,ny)
    do iter = 1, niter
       dumg = 0.0
       call vcycle(dumg,rmg,mu,nx,ny,10,4)
       umg = umg + dumg
       call resid(umg,vmg,rmg,mu,nx,ny)
    
    end do
    betab = exp(umg(wb:eb,sb:nb))
    betab = min(betab,maxbeta)
    betab = max(betab,minbeta)
 else
    betab = betar
 end if
 


 where (typ(wb:eb,sb:nb).eq.1)
    betab=shelfbeta
 end where
 
 where (typ(wb:eb,sb:nb).eq.2)
    betab=seabeta
 end where
 
! limit the maximum surface elevationf
 where (usrf .gt. maxsrf)
    usrf = maxsrf
    thck = usrf - lsrf
 end where
 
   
 maskb = lsrfb - topgb
 call ncdump2D(x(wb),y(sb),thckb,topgb,lsrfb,usrfb,usb,vsb,ucoef,divuh,divuhc,&
      umodb,betab,maskb,nx,ny)

  

end program t 
