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

  subroutine ncdump2D(x,y, thck, topg, lsrf, usrf, us, vs, umod, beta, temp, mask, n, m)
    
    integer n, m, nc_id, time_id, x1_id, y1_id, x0_id, y0_id, &
         topg_id, lsrf_id,usrf_id, thck_id, us_id, vs_id, umod_id, beta_id, temp_id,mask_id
    integer :: node_dim_ids(3), cell_dim_ids(3)
    real(kind=8), dimension(1:n) :: x
    real(kind=8), dimension(1:m) :: y
    ! DFM (12-18-11) Note that all quantities assumed cell-centered here.
    !cell centred arrays
    real(kind=8), dimension(1:n,1:m,1) :: thck, topg, lsrf, usrf,  & 
         us , temp, vs, umod, mask
    ! cell centred arrays
    real(kind=8), dimension(1:n,1:m,1) :: beta
    
    real(kind=8), dimension(1:n-1) :: xs
    real(kind=8), dimension(1:m-1) :: ys
    real(kind=8) :: time

    time = 0
    xs(1:n-1) = 0.5 * (x(1:n-1) + x(2:n))
    ys(1:m-1) = 0.5 * (y(1:m-1) + y(2:m))

    call nccheck( nf90_create("Antarctica-5km.test.nc", NF90_CLOBBER, nc_id) )

    call nccheck( nf90_def_dim(nc_id, "time", 1, node_dim_ids(3)) )
    call nccheck( nf90_def_dim(nc_id, "x1", n, node_dim_ids(1)) )
    call nccheck( nf90_def_dim(nc_id, "y1", m, node_dim_ids(2)) )
    
    cell_dim_ids(3) = node_dim_ids(3)
    call nccheck( nf90_def_dim(nc_id, "x0", n-1, cell_dim_ids(1)) )
    call nccheck( nf90_def_dim(nc_id, "y0", m-1, cell_dim_ids(2)) )

    call nccheck( nf90_def_var(nc_id, "time", nf90_real8, node_dim_ids(1), time_id))
    call nccheck( nf90_put_att(nc_id, time_id, "units", "a") )

    call nccheck( nf90_def_var(nc_id, "x1", nf90_real8, node_dim_ids(1), x1_id))
    call nccheck( nf90_def_var(nc_id, "y1", nf90_real8, node_dim_ids(2), y1_id))
    call nccheck( nf90_put_att(nc_id, x1_id, "units", "m") )
    call nccheck( nf90_put_att(nc_id, y1_id, "units", "m") )

    
    call nccheck( nf90_def_var(nc_id, "x0", nf90_real8, cell_dim_ids(1), x0_id))
    call nccheck( nf90_def_var(nc_id, "y0", nf90_real8, cell_dim_ids(2), y0_id))
    call nccheck( nf90_put_att(nc_id, x0_id, "units", "m") )
    call nccheck( nf90_put_att(nc_id, y0_id, "units", "m") )


    call nccheck( nf90_def_var(nc_id, "thk", nf90_real8, node_dim_ids, thck_id) )
    call nccheck( nf90_def_var(nc_id, "topg", nf90_real8, node_dim_ids, topg_id) )
    call nccheck( nf90_def_var(nc_id, "lsrf", nf90_real8, node_dim_ids, lsrf_id) )
    call nccheck( nf90_def_var(nc_id, "usrf", nf90_real8, node_dim_ids, usrf_id) )
    call nccheck( nf90_def_var(nc_id, "temp", nf90_real8, node_dim_ids, temp_id) )
    call nccheck( nf90_def_var(nc_id, "us", nf90_real8, node_dim_ids, us_id) )
    call nccheck( nf90_def_var(nc_id, "vs", nf90_real8, node_dim_ids, vs_id) )
    call nccheck( nf90_def_var(nc_id, "umod", nf90_real8, node_dim_ids, umod_id) )
!    call nccheck( nf90_def_var(nc_id, "beta", nf90_real8, cell_dim_ids, beta_id) )
    call nccheck( nf90_def_var(nc_id, "beta", nf90_real8, node_dim_ids, beta_id) )
    call nccheck( nf90_def_var(nc_id, "mask", nf90_real8, node_dim_ids, mask_id) )
    call nccheck( nf90_enddef(nc_id) )
    
    call nccheck( nf90_put_var(nc_id, time_id, time) )

    call nccheck( nf90_put_var(nc_id, x1_id, x) )
    call nccheck( nf90_put_var(nc_id, y1_id, y) )

    call nccheck( nf90_put_var(nc_id, x0_id, xs) )
    call nccheck( nf90_put_var(nc_id, y0_id, ys) )

    call nccheck( nf90_put_var(nc_id, thck_id, thck) )
    call nccheck( nf90_put_var(nc_id, topg_id, topg) )
    call nccheck( nf90_put_var(nc_id, lsrf_id, lsrf) )
    call nccheck( nf90_put_var(nc_id, usrf_id, usrf) )
    call nccheck( nf90_put_var(nc_id, temp_id, temp) )
    call nccheck( nf90_put_var(nc_id, us_id, us) )
    call nccheck( nf90_put_var(nc_id, vs_id, vs) )
    call nccheck( nf90_put_var(nc_id, umod_id, umod) )
    call nccheck( nf90_put_var(nc_id, beta_id, beta) )
    call nccheck( nf90_put_var(nc_id, mask_id, mask) )
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
  real (kind=8),parameter :: rhoi = 910.0d0, rhoo = 1028.0d0, grav = 9.81d0, cm = 1.0e-2
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
  
  integer, parameter :: ewn = 1200, nsn = 1200, unin=9
  integer, parameter :: nx = 1200, ny = 1200
  integer, parameter :: wb = 1, sb = nsn - ny + 1, eb = wb + nx - 1, nb = sb + ny -1
  integer, parameter :: m = 5, niter = 9
  real(kind=8),parameter :: scyr = 31556926.0d0, maxsrf=1250.0d0 
  real(kind=8),parameter :: rhoi = 910.0d0, rhoo = 1028.0d0, grav = 9.81d0
  real(kind=8),parameter :: dx = 5.0d+3, dy = dx,  tol = 1.0d-3
  real(kind=8),parameter :: minbeta = 5.0, maxbeta = 1.0e+8, cm = 1.0e-2
  real (kind = 8), dimension(1:ewn,1:nsn) :: topg, lsrf, usrf,  thck, us, vs, temp, mask, rhostar
  real (kind = 8), dimension(1:ewn,1:nsn) :: ttopg, tlsrf, tusrf,  tthck
  real (kind = 8), dimension(wb:eb,sb:nb) :: topgb, lsrfb, usrfb,  thckb, usb, vsb , tempb, maskb
  real (kind = 8), dimension(wb:eb,sb:nb) ::    umodb, betab,betar, dsx, dsy, umodsia
!  real (kind = 8), dimension(wb:eb-1,sb:nb-1) ::    betac
  real (kind = 8), dimension(wb:eb,sb:nb) ::    betac
  real (kind = 8), dimension(wb-1:eb+1,sb-1:nb+1) :: umg,vmg,dumg,rmg
  integer , dimension(1:ewn,1:nsn) ::  typ
  real (kind = 8), dimension(1:ewn) :: x
  real (kind = 8), dimension(1:nsn) :: y
  integer, parameter :: nre = 4
  real (kind = 8), dimension(-nre:nre,-nre:nre) :: ra
  real (kind = 8) :: mu, lambda, flwa, glen_n, res0,res1,res,drop0,drop1,drop
  real (kind = 8) :: resNorm
  integer :: ew, ns, iter
  
!  print *, wb, eb, sb, nb

  thck = 1.0

  write(*,*) 'reading thickness...'
  open(unin,file='thk.data',status='old')
  read(unin,*) ((thck(ew,ns), ew=1,ewn), ns=1,nsn)
  close(unin)

  write(*,*) 'reading topography...'
  open(unin,file='topg.data',status='old')
  read(unin,*) ((topg(ew,ns), ew=1,ewn), ns=1,nsn)
  close(unin)

  write(*,*) 'reading temperature...'
  open(unin,file='temp.data',status='old')
  read(unin,*) ((temp(ew,ns), ew=1,ewn), ns=1,nsn)
  close(unin)


  write (*,*) 'reading velocities...'
  open(unin,file='balvelmag.data',status='old')
  read(unin,*) ((umodb(ew,ns), ew=1,ewn), ns=1,nsn)
  close(unin)
!!!!  
!!!!  open(unin,file='usrf_5pg.data',status='old')
!!!!  read(unin,*) ((usrf(ew,ns), ew=1,ewn), ns=1,nsn)
!!!!  close(unin)
!!!!
!!!!  open(unin,file='uvel_5pg.data',status='old')
!!!!  read(unin,*) ((us(ew,ns), ew=1,ewn), ns=1,nsn)
!!!!  close(unin)
!!!!
!!!!  open(unin,file='vvel_5pg.data',status='old')
!!!!  read(unin,*) ((vs(ew,ns), ew=1,ewn), ns=1,nsn)
!!!!  close(unin)

  write(*,*) '... done reading'
  typ = 0 ! 0 for grounded ice

  x(1) = 0.0
!  print *, thck
  do ew = 2,ewn
     x(ew) = x(ew-1) + dx
  end do
  y(1) = 0.0
  do ns = 2,nsn
     y(ns) = y(ns-1) + dy
  end do
  
  !now use the bisicles expressions, 
  usrf = max(topg + thck, thck * (1.0 - rhoi/rhoo))
  lsrf = usrf - thck

!  call ncdump2D(x(wb),y(sb),thck,topg,lsrf,usrf,usb,vsb,umodb,betac,temp,maskb,nx,ny)



  !raw data - need to set up typ,mask,thck
  typ = 0 ! grounded ice
!  thck = usrf - lsrf
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
!  mask = typ


  thckb = thck(wb:eb,sb:nb)
  topgb = topg(wb:eb,sb:nb)
  tempb = temp(wb:eb,sb:nb)
  lsrfb = lsrf(wb:eb,sb:nb)
  usrfb = usrf(wb:eb,sb:nb)
  maskb = mask(wb:eb,sb:nb)
  usb = -us(wb:eb,sb:nb)
  vsb = -vs(wb:eb,sb:nb)
  call ncdump2D(x(wb),y(sb),thckb,topgb,lsrfb,usrfb,usb,vsb,umodb,betar,tempb,maskb,nx,ny)
  !stop
  !call ncdump2D(x,y,thck,topg,lsrf,usrf,usb,vsb,thck,thck,tempb,mask,ewn,nsn)
  !stop
 
  
  rhostar = rhoi
  where (lsrf .gt.  topg .and. thck .gt. cm)
     rhostar = (1.0-usrf/(thck))*rhoo
  end where

  write(*,*) sum (rhostar,mask =(lsrf .gt.  topg .and. thck .gt. cm) ) &
       /   count(mask =(lsrf .gt.  topg .and. thck .gt. cm) ), &
       maxval (rhostar,mask =(lsrf .gt.  topg .and. thck .gt. cm) ), &
       minval (rhostar,mask =(lsrf .gt.  topg .and. thck .gt. cm) )
  
  drop0 = 0.0
  drop =  15.2! Value from Anne Le Brocq
  drop1 = 2.0 * drop

  if (1.eq.0) then
     call testdrop(drop0,res0,usrf,lsrf,topg,thck,ewn,nsn)
     call testdrop(drop1,res1,usrf,lsrf,topg,thck,ewn,nsn)
     
     do iter = 1,20
        drop = (drop0 + drop1)/2
        call testdrop(drop,res,usrf,lsrf,topg,thck,ewn,nsn)
        write(*,*) "drop[",iter,"] <- ", drop, "; res[",iter,"] <- ", res
        
        if (res .gt. res0 .and. res .gt. res1) then
           write(*,*) 'failed to minimise res'
           stop
        else if (res0 .lt. res1) then
           if (res .gt. res0) then
              res1 = res
              drop1 = drop
           else
              res0 = res
              drop0 = drop
           end if
        else 
           if (res .gt. res1) then
              res0 = res
              drop0 = drop
           else
              res1 = res
              drop1 = drop
           end if
        end if
     end do
  end if
   drop =  15.2! Value from Anne Le Brocq

! DFM (11/3/11) -- don't worry about drop for the moment
!  call applydrop(drop,usrf,lsrf,topg,thck,ewn,nsn)
  
  !now use the bisicles expressions, 
  usrf = max(topg + thck, thck * (1.0 - rhoi/rhoo))
  lsrf = usrf - thck
 
  !do ns = 1,nsn
  !   do ew = 1,ewn
  !      if (lsrf(ew,ns) - topg(ew,ns) .gt. 0.0001)  then
  !         thck(ew,ns) = 0.0
  !         usrf(ew,ns) = thck(ew,ns) * (1028.0d0 - 910.0d0)/1028.0d0
  !         lsrf(ew,ns) = usrf(ew,ns)  - thck(ew,ns) 
  !      end if
  !   end do
  !end do

  !throw away any floating ice less than 120 m thick in more then
  !120 m of ocean
  !do ew = 1,ewn
  !   do ns = 1,nsn
  !      if ((lsrf(ew,ns) - topg(ew,ns) .gt. 120.0) &
  !           .and. (thck(ew,ns) .lt. 120.0)) then 
  !         typ(ew,ns) = 2 ! 2 for open sea
  !         thck(ew,ns) = 0.0
  !         usrf(ew,ns) = thck(ew,ns) * (1028.0d0 - 910.0d0)/1028.0d0
  !         lsrf(ew,ns) = usrf(ew,ns)  - thck(ew,ns) 
  !     end if
  !  end do
  !end do
 
  thckb = thck(wb:eb,sb:nb)
  topgb = topg(wb:eb,sb:nb)
  lsrfb = lsrf(wb:eb,sb:nb)
  usrfb = usrf(wb:eb,sb:nb)
  maskb = mask(wb:eb,sb:nb)
  usb = -us(wb:eb,sb:nb)
  vsb = -vs(wb:eb,sb:nb)


  !attempt to get rid of the missing values
  do ew = wb,eb
     do ns = sb,nb
        if ( abs(usb(ew,ns)) .gt. 9000.0) then
           usb(ew,ns) = usb(ew,ns-1)
           vsb(ew,ns) = vsb(ew,ns-1)
        end if
     end do
  end do
  where (typ(wb:eb,sb:nb) .eq. 2)
     usb = 0.0
     vsb = 0.0
     umodb = 0.0
  end where

!  umodb = sqrt(usb**2 + vsb**2)

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

 usb = dsx
 vsb = dsy

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

! DFM -- I really don't understand what's going on here...
 ! do ns = sb+100,nb
 !    do ew = wb,eb
 !    if (typ(ew,ns) .gt. 0) then
 !       call random_number(ra)
 !       ra = ra + 1.0e-2 * ( 1.0e+2 + umodb(ew-nre:ew-nre,(ns-1-2*nre):(ns-1) ))
 !       ra = ra / sum(ra)

 !       betar(ew,ns) =  sum(ra * betar(ew-nre:ew-nre,(ns-1-2*nre):(ns-1) )) 
 !    end if
 !    end do
 ! end do

 betab = betar
 call ncdump2D(x(wb),y(sb),thckb,topgb,lsrfb,usrfb,usb,vsb,umodb,betab,temp,maskb,nx,ny)
! stop
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
!    write(*,*) niter 
    resNorm = sum(abs(rmg) )
    write(*,*) 'initial res= ', resNorm
    do iter = 1, niter
       dumg = 0.0
       call vcycle(dumg,rmg,mu,nx,ny,10,4)
       umg = umg + dumg
       call resid(umg,vmg,rmg,mu,nx,ny)
       resNorm = sum(abs(rmg) )
       write (*,*) 'MG iter', iter, " norm = ", resNorm
    
    end do
    betab = exp(umg(wb:eb,sb:nb))
    betab = min(betab,maxbeta)
    betab = max(betab,minbeta)
 else
    betab = betar
 end if
 !interpolate to cell centers
 

 betac = betab
! betac(wb:eb-1,sb:nb-1) = 0.25d0 * ( betab(wb:eb-1,sb:nb-1) & 
!      + betab(wb+1:eb,sb:nb-1) &
!      + betab(wb+1:eb,sb+1:nb)  &
!      + betab(wb:eb-1,sb+1:nb)) 

 ! limit the maximum surface elevationf
 where (usrf .gt. maxsrf)
    usrf = maxsrf
    thck = usrf - lsrf
 end where
 
   
 if (1.eq.0) then
 !set the outer m ring to be uniform in the normal direction
 ! and ensure the east and west sides are grounded
 ! left(west),right(east),bottom(north)
 do ns = sb,nb
    do ew = wb,wb+m
       thckb(ew,ns) = thckb(wb+m+1,ns)
       usrfb(ew,ns) = usrfb(wb+m+1,ns)
       lsrfb(ew,ns) = lsrfb(wb+m+1,ns)
       topgb(ew,ns) = topgb(wb+m+1,ns)

       thckb(ew,ns) = max(thckb(ew,ns), -topgb(ew,ns) * rhoo/rhoi + 10.0)
       lsrfb(ew,ns) = topg(ew,ns)
       usrfb(ew,ns) = topg(ew,ns) +  thckb(ew,ns)
    end do
 end do
       
 do ns = sb,nb
    do ew = eb-m,eb
       thckb(ew,ns) = thckb(eb-m-1,ns)
       usrfb(ew,ns) = usrfb(eb-m-1,ns)
       lsrfb(ew,ns) = lsrfb(eb-m-1,ns)
       topgb(ew,ns) = topgb(eb-m-1,ns)
    end do
 end do


 do ns = sb,sb + m
    do ew = wb,eb
       thckb(ew,ns) = thckb(ew,sb+m+1)
       usrfb(ew,ns) = usrfb(ew,sb+m+1)
       lsrfb(ew,ns) = lsrfb(ew,sb+m+1)
       topgb(ew,ns) = topgb(ew,sb+m+1)
    end do
 end do

 do ns = nb - m,nb
    do ew = wb,eb
       thckb(ew,ns) = thckb(ew,nb-m-1)
       usrfb(ew,ns) = usrfb(ew,nb-m-1)
       lsrfb(ew,ns) = lsrfb(ew,nb-m-1)
       topgb(ew,ns) = topgb(ew,nb-m-1)
    end do
 end do
end if
if (1.eq.0) then
  !cut off an annoying chunk near wb,bn
  thckb(wb:wb+5,nb-9:nb) = 0.0
  lsrfb(wb:wb+5,nb-9:nb) = 0.0
  usrfb(wb:wb+5,nb-9:nb) = 0.0
  typ(wb:wb+5,nb-9:nb) = 2
end if

  maskb = lsrfb - topgb
  
  !do ns = 100,ns-20
  !   where(thckb(wb:eb,ns-1) .gt. 100.0)
  !      thckb(wb:eb,ns) = max(thckb(wb:eb,ns),1.0)
  !   end where
  !end do

  !thckb = max(thckb,10.0)
  

 call ncdump2D(x(wb),y(sb),thckb,topgb,lsrfb,usrfb,usb,vsb,umodb,betac,temp,maskb,nx,ny)

  

end program t 
