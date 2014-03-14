!tarmac.f90 : A collection of subroutines which are entirely ignorant
!of the horizontal mesh. 
! Copyright (C) 2012 BISICLES contributors
!
! Please refer to Copyright.txt, in Chombo's root directory.
!
module glunnamed
  implicit none

  !physical constants, all SI units
  real(kind=8) scyr,  rhoi , rhoo,  grav , shci, lhci, coni , pmlt, trpt
  ! seconds in a year, density of ice, density of water, acceleration due to gravity 
  ! specific heat capacity of ice, latent heat of fusion of water, thermal conductivity of ice,  
  ! factor for dependence of melting point on pressure, triple point of water (K)

  integer, parameter :: groundedmaskval = 1, floatingmaskval = 2, openseamaskval = 4, openlandmaskval = 8
  real(kind=8), parameter :: temp_eps = 1.0e-3
  
  contains

    real(kind=8) function fo_upwind(u,sminus,splus)
      !first-order upwind scheme.
      real(kind=8), intent(in) :: u,sminus,splus
      if (u .ge. 0.0) then
         fo_upwind = u * sminus
      else
         fo_upwind = u * splus 
      end if
      return
    end function fo_upwind


    subroutine tdmasolve(n,x,a,b,c,r)
      ! solve the n x n tridiagonal system
      ! overwriting the coefficients a,b,c and r
      ! diagonal elements must not be zero
      ! |b1 c1 0 0 0  ...    0| |x1|   |r1|
      ! |a2 b2 c2 0 0 0  ... 0| |x2| = |r2|
      ! |0  a3 b3 c3 0 ......0| |x3|   |r3|
      ! etc      
      integer, intent(in) :: n
      real(kind=8), dimension(1:n), intent(inout) :: x,b,r
      real(kind=8), dimension(1:n), intent(in) :: a,c
      !locals
      integer :: i
      real(kind=8) :: t

      do i = 2,n
         t = a(i)/b(i-1)
         b(i) = b(i) - t*c(i-1)
         r(i) = r(i) - t*r(i-1)
      end do

      x(n) = r(n)/b(n)
      
      do i = n-1, 1,-1
         x(i) = (r(i)-c(i)*x(i+1))/b(i)
      end do
     

    end subroutine tdmasolve

    subroutine sigma_advect(rhs, temp, stemp, btemp, husig ,  & 
         fsig, dt, mask, n)
      !increment rhs with interlayer advection
      !TODO : replace first order upwind fluxes with something nice
      integer :: n
      real(kind=8), dimension(1:n), intent(inout) :: rhs
      real(kind=8), dimension(1:n), intent(in) :: temp ! temperature
      real(kind=8), dimension(1:n+1), intent(in) :: husig ! thickness * u^{\sigma}
      real(kind=8), dimension(1:n+1), intent(in) :: fsig ! sigma at layer faces
      real(kind=8), intent(in) :: stemp, btemp, dt ! surface & basal temperature, time step
      integer, intent(in) :: mask
      !locals
      integer l
      real(kind=8) flux
      
      !flux across the surface : Dirichlett condition
      flux = dt * (fo_upwind(husig(1) ,stemp,temp(1)) )
      rhs(1) = rhs(1) + flux / (fsig(2)-fsig(1))

      do l = 2, n
         flux = dt * (fo_upwind(husig(l) ,temp(l-1),temp(l)) )
         rhs(l-1) = rhs(l-1) - flux / (fsig(l)-fsig(l-1))
         rhs(l) = rhs(l) + flux / (fsig(l+1)-fsig(l))
      end do

      !flux across the base
      flux = dt * (fo_upwind(husig(n+1) ,temp(n),btemp) )
!!$      if (mask.eq.floatingmaskval) then
!!$         ! Dirichlett condition
!!$         flux = dt * (fo_upwind(husig(n+1) ,temp(n),btemp) )
!!$      else
!!$         ! no advection flux
!!$         flux =  0.0
!!$      end if
      rhs(n) = rhs(n) - flux / (fsig(n+1)-fsig(n))

    end subroutine sigma_advect

  
    subroutine fo_diffusive_advance(temp,stemp,sflux,sdiric,btemp,bflux,rhs,chi, &
         thckold,thcknew,fsig,dt,mask,n)
      ! solve the equation 
      ! H_nT_n - dt * 1/H * \chi * d^2 T_n / dsigma^2 = H_To + rhs
      ! with either a Dirichlett boundary condition T = stemp 
      ! or a flux boundary condition - chi dT/dz = sflux at sigma = 0
      ! and a flux boundary condition - chi dT/dz = bflux at sigma = 1
      ! (first order - backward Euler integration)
      integer :: n
      real(kind=8), intent(inout) :: stemp,btemp,sflux,bflux
      real(kind=8), intent(in) :: chi,thckold,thcknew,dt
      real(kind=8), dimension(1:n), intent(in) :: rhs
      real(kind=8), dimension(1:n), intent(inout) :: temp
      real(kind=8), dimension(1:n+1), intent(in) :: fsig ! sigma at layer faces
      integer, intent(in) :: mask
      logical, intent(in) :: sdiric
      !locals
      real(kind=8), dimension(1:n) :: a,c,b,r ! TDMA coefficients
      real(kind=8), dimension(1:n) :: cdsig
      real(kind=8), dimension(1:n+1) :: fdsig
      real(kind=8) :: chistar, tpmp, eps, bmb, dtempdz
      logical :: bdiric
      integer :: i
      
      eps = 1.0e-3 


      do i = 1,n
         cdsig(i) = fsig(i+1)-fsig(i)
      end do
      do i = 2,n
         fdsig(i) = 0.5 * cdsig(i) + 0.5 * cdsig(i-1)
      end do
      fdsig(1) = cdsig(1)
      fdsig(n+1) = cdsig(n)
     
      chistar = 2.0 * dt *chi / (thcknew+thckold) 
      !interior layers 2 to n-1
      do i = 2, n-1
         a(i) = -chistar / (cdsig(i)*fdsig(i+1))
         c(i) = -chistar / (cdsig(i)*fdsig(i))
         b(i) = thcknew - a(i) - c(i) 
         r(i) = thckold * temp(i) + rhs(i)
      end do

      

      ! top layer
      ! hold temperature eps K below the pressure melting point
      ! top layer : Dirichlett boundary condition temp = stemp
      tpmp = trpt -  eps
      if (stemp .gt. tpmp) stemp = tpmp
      
      if (sdiric) then
         !Dirichlett boundary condition
         a(1) = 0.0d0
         b(1) =  thcknew + 3. * chistar / (cdsig(1)*fdsig(1))
         c(1) = - chistar / (cdsig(1)*fdsig(1))
         r(1) = thckold * temp(1) +  rhs(1) + &
              2. * stemp * chistar / (cdsig(1)*fdsig(1))
      else
         !Flux boundary condition
         a(1) = 0.0d0
         c(1) =  - chistar / (cdsig(1)*fdsig(1))
         b(1) = thcknew - c(1)
         r(1) = thckold * temp(1) +  rhs(1)  + dt * sflux/cdsig(1) 
      end if

      ! bottom layer
      ! hold temperature eps K below the pressure melting point
      tpmp = trpt - pmlt * rhoi * grav * thcknew - eps
      bdiric = ((mask.eq.floatingmaskval) .or. (btemp .ge. tpmp))

      if (bdiric) then
         ! ice at pressure melting point 
         tpmp = trpt - pmlt * rhoi * grav * thcknew - eps
         btemp = tpmp
         ! Dirichlett condition, temp = btemp
         a(n) =  - chistar / (cdsig(n)*fdsig(n+1))
         b(n) = thcknew + 3. * chistar / (cdsig(n)*fdsig(n+1))
         c(n) = 0.0
         r(n) =  thckold * temp(n) +  rhs(n) + &
              2. * btemp * chistar / (cdsig(n)*fdsig(n))
      else
         !ice above pressure melting point, 
         !heat flux condition
         a(n) =  - chistar / (cdsig(n)*fdsig(n+1))
         b(n) = thcknew - a(n)
         c(n) = 0.0
         r(n) = thckold * temp(n) +  rhs(n)  + dt * bflux/cdsig(n) 
      end if
      
      call tdmasolve(n,temp,a,b,c,r)

   
      if (sdiric) then
         sflux = 2.0d0 * chistar * (temp(1) - stemp) /  cdsig(1)
      else
         stemp = 1.5 * temp(1) - 0.5 * temp(2)
         tpmp  = trpt - eps
         if (stemp .gt. tpmp) stemp = tpmp
      end if

      if (bdiric) then
         bflux = - 2.0d0 * chistar * (temp(n) - btemp) /  cdsig(n)
      else
         btemp = 1.5 * temp(n) - 0.5 * temp(n-1)
         tpmp = trpt - pmlt * rhoi * grav * thcknew  - eps
         if (btemp .gt. tpmp) then 
            btemp = tpmp
            !compute basal melt rate
            dtempdz = temp(n) + pmlt * rhoi * grav * thcknew * 0.5*(fsig(i) + fsig(i+1))
            dtempdz = 2.0 / thcknew * (dtempdz - trpt) / cdsig(n)
            bmb = (coni * scyr * dtempdz + bflux * (shci*rhoi))/(rhoi * lhci)
         end if
      end if

      !This is a kludge, and should be replaced with something conservative
      do i = 1, n
         tpmp = trpt - pmlt * rhoi * grav * thcknew * 0.5*(fsig(i) + fsig(i+1)) - eps
         if (temp(i) .gt. tpmp) temp(i) = tpmp 
      end do

      return

    end subroutine fo_diffusive_advance

    subroutine set_constants( ascyr , arhoi , arhoo, agrav , ashci, alhci, aconi , apmlt, atrpt)
      
      implicit none
      real(kind=8) ascyr,  arhoi , arhoo,  agrav , ashci, alhci, aconi , apmlt, atrpt
      scyr = ascyr
      rhoi = arhoi 
      rhoo = arhoo
      grav = agrav 
      shci = ashci
      lhci = alhci
      coni = aconi 
      pmlt = apmlt
      trpt = atrpt
    end subroutine set_constants

end module glunnamed


subroutine tarmac_set_constants( ascyr , arhoi , arhoo, agrav , ashci, alhci, aconi , apmlt, atrpt)
  use glunnamed
  implicit none
  real(kind=8) ascyr,  arhoi , arhoo,  agrav , ashci, alhci, aconi , apmlt, atrpt
  !set_constants(scyr , rhoi , rhoo, grav , shci, lhci, coni , pmlt, trpt)
    scyr = ascyr
    rhoi = arhoi 
    rhoo = arhoo
    grav = agrav 
    shci = ashci
    lhci = alhci
    coni = aconi 
    pmlt = apmlt
    trpt = atrpt
end subroutine tarmac_set_constants

subroutine tarmac_update_temperature(temp, stemp, sflux, sdiric, btemp, mask, &
     bflux, rhs, thckold, thcknew, usig, fsig, dsig, time, dt, n)

  !update the mid-layer temperarures temp, surface temparature
  !stemp and the basal temperature btemp
  !
  !There are n layers, with sigma = 0 at the surface, 
  !sigma = sigma(i) at the middle of layer i and sigma = 1
  !at the base of the ice
  ! 
  use glunnamed

  implicit none
  integer :: n
  real(kind=8), dimension(1:n), intent(inout) :: temp,rhs,dsig
  real(kind=8), dimension(1:n+1), intent(inout) :: usig,fsig
  real(kind=8), intent(inout) :: stemp,btemp,bflux,sflux
  real(kind=8), intent(in) :: time,dt,thckold,thcknew
  integer, intent(in) :: mask
  logical, intent(in) :: sdiric
  !locals
  
  real(kind=8) :: chi, dtcfl,tthcknew,tthckold
  real(kind=8), dimension(1:n) :: rhsl
  
  integer l,nt,it,i,tt
 
  btemp=min(trpt - pmlt * rhoi * grav * thcknew - temp_eps,btemp)
  stemp=min(trpt - temp_eps,stemp)


  if (thcknew.gt.1.0d0) then
     chi = coni/(shci * rhoi) *  scyr 
     !work out a stable time step
     dtcfl = dt
     do i = 1,n
        dtcfl = min(dtcfl,(thckold+thcknew)*(fsig(i+1)-fsig(i))/(abs(usig(i+1)) + abs(usig(i))))
     end do
     nt = min(100,ceiling(dt/dtcfl))
     dtcfl = dt/dble(nt)
     tt = 0.0d0
     do it = 1,nt
        
        rhsl = rhs/dble(nt)
        
        !modify rhs to include interlayer advection across constant sigma faces
        call sigma_advect(rhsl, temp, stemp, btemp, usig , &
             fsig, dtcfl, mask,  n)
        !compute vertical diffusion & update the basal temperature
        tthckold = tt/dt * thcknew + (1.0d0-tt/dt) *thckold
        tt = tt + dtcfl
        tthcknew = tt/dt * thcknew + (1.0d0-tt/dt) *thckold
        call fo_diffusive_advance(temp, stemp, sflux, sdiric, btemp, bflux, rhsl, chi, &
             tthckold, tthcknew,fsig,dtcfl,mask,n)
       
     end do
  else
      !no ice
      if (sdiric) then
         sflux = bflux
     else
        stemp = trpt - temp_eps
     end if
     temp = stemp 
     where (temp > trpt - temp_eps)
        temp = trpt - temp_eps
     end where
     btemp = min(trpt - temp_eps,stemp)
  end if

  return

end subroutine tarmac_update_temperature

subroutine tarmac_compute_zvel(uz,uzs,husig,ux,uy, divhu, &
     fsig, csig, dsig, n, thck, dsx, dhx, dsy, dhy, dst, dht, & 
     smb, bmb)
  !compute the vertical velocity from the divergence
  !of horizontal velocity and the kinematic boundary
  !conditions. Also compute the contravariant \sigma 
  !component of the velocity u^{\sigma} = u.e^{\sigma}
  !which appears in cross-layer advection
  !
  !This means solving a first order ODE w' = f with
  !known w(s) and w(b), which is ill-posed unless 
  !w(s) - w(b) = integral(f dz,s,b) - which should 
  !be the case. To check this, we compute both 
  !w(1) - the surface velocity computed by integration, and
  !ws - the surface velocity boundary condition
  implicit none
  integer, intent(in) :: n
  real(kind=8), dimension(1:n+1), intent(in) :: ux,uy ! horizontal velocity at layer faces
  real(kind=8), dimension(1:n+1), intent(out) :: uz,husig ! vertical velocity & h*sigma velocity at   layer faces
  real(kind=8), intent(out) :: uzs ! uz at surface by kinematic bc
  real(kind=8), dimension(1:n+1), intent(in) :: fsig ! sigma at layer faces
  real(kind=8), dimension(1:n), intent(in) :: divhu ! horizontal part of div(Hu)
  real(kind=8), dimension(1:n), intent(in) :: csig, dsig ! sigma, dsigma at layer midoints
  real(kind=8), intent(in) :: thck,dsx,dsy,dhx,dhy,dst,dht !H, ds/dx, ds/dy, dH/dx, dH/dy
  real(kind=8), intent(in) :: smb, bmb ! surface & basal mass balance rates
  !locals
  integer :: i
  real(kind=8), dimension(1:n+1) :: gfx,gfy,gft
  real(kind=8) :: uzref
  !geoemtry factors
  gfx(1:n+1) = (dsx-fsig(1:n+1)*dhx)
  gfy(1:n+1) = (dsy-fsig(1:n+1)*dhy)
  gft(1:n+1) = (dst-fsig(1:n+1)*dht)

  !velocity at the base
  uz(n+1) =  gft(n+1) + ux(n+1)*gfx(n+1) + uy(n+1)*gfy(n+1) + bmb
  husig(n+1) =  - bmb + gft(n+1)

 
  !velocity at the surface
  uzs = dst + ux(1)*dsx + uy(1)*dsy - smb

  

  !integration
  do i = n , 1, -1
     uz(i) = uz(i+1) - divhu(i)*dsig(i) & 
          + (gfx(i)*ux(i) - gfx(i+1)*ux(i+1)) &
          + (gfy(i)*uy(i) - gfy(i+1)*uy(i+1))  
     
     !h* u^{\sigma}
     husig(i) =  husig(i+1) + (divhu(i) + dht)*dsig(i) 
     
  end do
  

  return

end subroutine tarmac_compute_zvel
