! column_themodynamics.f90 : f90 subroutines used to carry 
! out thermodynamics calculations that do not depend on the horizontal mesh
! Copyright (C) 2012-2015 BISICLES contributors
!
! Please refer to Copyright.txt, in Chombo's root directory.
!
module column_thermodynamics
  implicit none

  !physical constants, all SI units
  real(kind=8) scyr,  rhoi , rhoo,  grav , shci, lhci, coni , conm,  pmlt, trpt
  ! seconds in a year, density of ice, density of water, acceleration due to gravity 
  ! specific heat capacity of ice, latent heat of fusion of water, thermal conductivity of cold ice,  
  ! regularizing conductivity for temperate ice, 
  ! factor for dependence of melting point on pressure, triple point of water (K)

  integer, parameter :: groundedmaskval = 1, floatingmaskval = 2, openseamaskval = 4, openlandmaskval = 8
  real(kind=8), parameter :: temp_eps = 1.0e-3
  real(kind=8), parameter :: zero_debug = 0.0d0
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


  subroutine tdmasolve(n,x,a,b,c,d)
    ! solve the n x n tridiagonal system
    ! overwriting the coefficients a,b,c and d
    ! diagonal elements must not be zero
    ! |b1 c1 0 0 0  ...    0| |x1|   |d1|
    ! |a2 b2 c2 0 0 0  ... 0| |x2| = |d2|
    ! |0  a3 b3 c3 0 ......0| |x3|   |d3|
    ! etc      
    integer, intent(in) :: n
    real(kind=8), dimension(1:n), intent(inout) :: x,a,b,c,d
    !locals
    integer :: i
    real(kind=8) :: t

    c(1) = c(1)/b(1)
    d(1) = d(1)/b(1)

    do i = 2,n
       t = b(i) - a(i)*c(i-1)
       c(i) = c(i) / t
       d(i) = (d(i)-a(i)*d(i-1))/t
    end do

    x(n) = d(n)

    do i = n-1, 1,-1
       x(i) = (d(i)-c(i)*x(i+1))
    end do


  end subroutine tdmasolve

  subroutine sigma_advect(rhs, energy, senergy, benergy, husig ,  & 
       fsig, dt, mask, n)
    !increment rhs with interlayer advection
    !TODO : replace first order upwind fluxes with something nice
    integer :: n
    real(kind=8), dimension(1:n), intent(inout) :: rhs
    real(kind=8), dimension(1:n), intent(in) :: energy ! internal energy density
    real(kind=8), dimension(1:n+1), intent(in) :: husig ! thickness * u^{\sigma}
    real(kind=8), dimension(1:n+1), intent(in) :: fsig ! sigma at layer faces
    real(kind=8), intent(in) :: senergy, benergy, dt ! surface & basal temperature, time step
    integer, intent(in) :: mask
    !locals
    integer l
    real(kind=8) flux

    !flux across the surface : Dirichlett condition
    flux = dt * (fo_upwind(husig(1) ,senergy,energy(1)) )
    rhs(1) = rhs(1) + flux / (fsig(2)-fsig(1))

    do l = 2, n
       flux = dt * (fo_upwind(husig(l) ,energy(l-1),energy(l)) )
       rhs(l-1) = rhs(l-1) - flux / (fsig(l)-fsig(l-1))
       rhs(l) = rhs(l) + flux / (fsig(l+1)-fsig(l))
    end do

    !flux across the base
    flux =  dt * (fo_upwind(husig(n+1) ,energy(n),  benergy) )
    rhs(n) = rhs(n) - flux / (fsig(n+1)-fsig(n))

  end subroutine sigma_advect


  subroutine fo_diffusive_advance(energy,senergy,sflux,sdiric, & 
       benergy,bflux,rhs,thckold,thcknew,fsig,dt,mask,n)
    ! solve the equation 
    ! H*E - dt * 1/H * d/dsigma (q) = Ho*Eo + rhs
    !
    !Non-advective fluxes follow Aschwanden 2012, with
    !     { - (Ki + Kr) dE/dz               if E < Epmp
    ! q = {   
    !     { - Ki d/dz(Epmp) - Kr dE/dz  if E > Epmp
    !   
    ! Ki = ki/ci and Kr = km/ci
    ! ki is the ice conductivity, km is the (small) mositure conductivity ci is the specific heat capacity of ice, 
    ! Tpmp is the pressure melting point, Epmp = ci * Tpmp, kr is a regularization parameter 
    ! with either a Dirichlett boundary condition E = sE  
    ! or a flux boundary condition - chi(z) dE/dsigma = sflux at sigma = 0
    ! and a flux boundary condition - chi(sigma) dE/dsigma = bflux at sigma = 1
    !
    !
    ! (first order - backward Euler integration)
    implicit none
    integer :: n
    real(kind=8), intent(inout) :: senergy,benergy,sflux,bflux
    real(kind=8), intent(in) :: thckold,thcknew,dt
    real(kind=8), dimension(1:n), intent(in) :: rhs
    real(kind=8), dimension(1:n), intent(inout) :: energy
    real(kind=8), dimension(1:n+1), intent(in) :: fsig ! sigma, conductivity at layer faces
    integer, intent(in) :: mask
    logical, intent(in) :: sdiric
    !locals
    real(kind=8), dimension(1:n) :: a,c,b,r,a2,c2,b2,r2 ! TDMA coefficients
    real(kind=8), dimension(1:n) :: csig,cdsig,epmp,drain,energy2
    real(kind=8), dimension(1:n+1) :: fdsig,  kcc, ktt, kct, ktc
    real(kind=8) :: bmb, eps,k, kr,bepmp,sepmp, kb
    integer :: i
    logical :: btemperate
   
    eps = 1.0e-3 

    do i = 1,n
       cdsig(i) = fsig(i+1)-fsig(i)
       csig(i) = 0.5d0*(fsig(i+1)+fsig(i))
    end do
    do i = 2,n
       fdsig(i) = 0.5d0 * cdsig(i) + 0.5d0 * cdsig(i-1)
    end do
    fdsig(1) = cdsig(1)
    fdsig(n+1) = cdsig(n)

    !pressure melting point
    epmp = shci * (trpt - pmlt * 0.5d0 * (thckold+thcknew) * csig * rhoi * grav) ! at layer midpoints
    sepmp = shci * trpt !at surface
    bepmp = shci * (trpt - pmlt * 0.5d0 * (thckold+thcknew) * rhoi * grav) !at base

  

    ! conduction coeffs
    kr =  2.0 * dt / (thcknew+thckold) * conm/(shci * rhoi) *  scyr
    k = 2.0 * dt / (thcknew+thckold) * coni/(shci * rhoi) *  scyr - kr
    !face centered conduction coeffcients
    kcc = 0.0d0
    ktt = 0.0d0 
    kct = 0.0d0
    ktc = 0.0d0

    !regularization
    kcc(1:n+1) =  kcc(1:n+1) + kr /  fdsig(1:n+1)

    !cold-cold interfaces
    where ( (energy(1:n-1).lt.epmp(1:n-1)) .and. (energy(2:n).lt.epmp(2:n)) )
       kcc(2:n) = kcc(2:n) + k /  fdsig(2:n)
    end where

    !cold-temperate interfaces
    where ( (energy(1:n-1).lt.epmp(1:n-1)) .and. (energy(2:n).ge.epmp(2:n)) )
       kct(2:n) =  kct(2:n) +  k /  fdsig(2:n)
    end where
    
    !temperate-cold interfaces
    where ( (energy(1:n-1).ge.epmp(1:n-1)) .and. (energy(2:n).lt.epmp(2:n)) )
       ktc(2:n) =  ktc(2:n) +  k /  fdsig(2:n)
    end where

    !temperate-temperate interfaces
    where ( (energy(1:n-1).ge.epmp(1:n-1)) .and. (energy(2:n).ge.epmp(2:n)) )
       ktt(2:n) =  ktt(2:n) +  k /  fdsig(2:n)
    end where

    !TDMA coefficients
    !interior ice layers 2 to n-1
    do i = 2, n-1
       a(i) =  (- kcc(i)   - kct(i)  )/cdsig(i)
       c(i) =  (- kcc(i+1) - ktc(i+1))/cdsig(i)
       b(i) = thcknew + (kcc(i) + ktc(i) + kcc(i+1) + kct(i+1) )/cdsig(i)
       r(i) = thckold * energy(i) + rhs(i) &
            - epmp(i-1)/cdsig(i)* (-ktt(i) - ktc(i)) &
            - epmp(i+1)/cdsig(i)* (-ktt(i+1) - kct(i+1)) &
            - epmp(i)/cdsig(i) * (ktt(i) + kct(i) + ktt(i+1) + ktc(i+1))
    end do

    ! top layer
    i = 1
    if (sdiric) then
       !Dirichlett boundary condition, E(0) + E(1) = 2 Es. 
       if (senergy.lt.sepmp) then
          !cold surface
          if ( energy(i).lt.epmp(i) ) then
             !cold-cold
             kcc(i) = kcc(i) + k /  fdsig(i)
          else
             !cold-temperate
             kct(i) =  kct(i) +  k /  fdsig(i)
          end if
       else
          !temperate surface
           if ( energy(i).lt.epmp(i) ) then
             !temperate-cold
             ktc(i) = ktc(i) + k /  fdsig(i)
          else
             !temperate-temperate
             ktt(i) =  ktt(i) +  k /  fdsig(i)
          end if   
       end if
      
       a(i) =  0.0d0
       c(i) =  (- kcc(i+1) - ktc(i+1))/cdsig(i)
       b(i) = thcknew + (2.0d0 * kcc(i) + ktc(i) + kct(i) & 
            + kcc(i+1) + kct(i+1) )/cdsig(i)
       r(i) = thckold * energy(i) + rhs(i) &
            + 2.0d0*senergy/cdsig(i)  * (kcc(i) + kct(i)) & 
            - sepmp/cdsig(i)* (-ktt(i) - ktc(i)) &
            - epmp(i+1)/cdsig(i)* (-ktt(i+1) - kct(i+1)) ! &
            !- epmp(i)/cdsig(i) * (ktt(i) + kct(i) + ktt(i+1) + ktc(i+1)) &
                
    else
       !Flux boundary condition
       a(i) =  0.0d0
       c(i) =  (- kcc(i+1) - ktc(i+1))/cdsig(i)
       b(i) = thcknew + (kcc(i+1) + kct(i+1) )/cdsig(i)
       r(i) = thckold * energy(i) + rhs(i) &
            - epmp(i+1)/cdsig(i)* (-ktt(i+1) - kct(i+1)) &
            - epmp(i)/cdsig(i) * ( ktt(i+1) + ktc(i+1)) & 
            +  dt * sflux/cdsig(1)
    end if

    ! bottom layer, flux boundary condition
    i = n 
   
    btemperate = (mask.eq.floatingmaskval)
    if (btemperate) then
       !basal flux will be kb/dt * (bepmp - energy(i))
       kb = 2.0d0 * (k + kr) / ( fdsig(i+1))
       a(i) = ( -kcc(i) - kct(i) )/cdsig(i)
       c(i) = 0.0d0
       b(i) = thcknew + (kcc(i) + ktc(i) + kb)/ cdsig(i) 
       r(i) = thckold * energy(i) + rhs(i) &
            - epmp(i-1)/cdsig(i)* (-ktt(i) - ktc(i)) &
            - epmp(i)/cdsig(i) * (ktt(i) + kct(i)) &
            + bepmp/cdsig(i) * kb
    else
       a(i) = ( -kcc(i) - kct(i) )/cdsig(i)
       c(i) = 0.0d0
       b(i) = thcknew + (kcc(i) + ktc(i))/cdsig(i) 
       r(i) = thckold * energy(i) + rhs(i) &
            - epmp(i-1)/cdsig(i)* (-ktt(i) - ktc(i)) &
            - epmp(i)/cdsig(i) * (ktt(i) + kct(i)) &
            +  dt * bflux/cdsig(i)
       
    end if


    a2 = a
    b2 = b
    c2 = c
    r2 = r
    call tdmasolve(n,energy2,a2,b2,c2,r2)
    
    !drainage
    where (energy2 .gt. epmp + 0.01 * lhci)
       drain = 0.5d0 * dt
    elsewhere
       drain = 0.0d0
    end where

    b = b + thcknew*drain
    r = r + thcknew*drain*epmp
    call tdmasolve(n,energy,a,b,c,r)
    
    drain = (energy2 - energy)/lhci
 
    !compute surface energy density or heat flux
    if (sdiric) then
       sflux = 0.0 ! FIX ME
    else
       senergy = 1.5d0 * energy(1) - 0.5d0 * energy(2)
    end if

    !compute basal energy density and modify heat flux
    if (btemperate) then
       bflux = kb/dt * (bepmp - energy(i))
       benergy = bepmp
    else
       kb = coni/(shci * rhoi) *  scyr
       benergy = energy(n) + bflux / kb * thcknew * 0.5d0 * (fsig(n+1)-fsig(n))
    end if

    !basal melt rate
    bmb = sum(drain) / dt * cdsig(n) * thcknew

    return

  end subroutine fo_diffusive_advance

end module column_thermodynamics


subroutine column_thermodynamics_set_constants( ascyr , arhoi , arhoo, agrav , ashci, alhci, aconi , aconm, apmlt, atrpt)
  use column_thermodynamics
  implicit none
  real(kind=8) ascyr,  arhoi , arhoo,  agrav , ashci, alhci, aconi , aconm, apmlt, atrpt
  !set_constants(scyr , rhoi , rhoo, grav , shci, lhci, coni , pmlt, trpt)
  scyr = ascyr
  rhoi = arhoi 
  rhoo = arhoo
  grav = agrav 
  shci = ashci
  lhci = alhci
  coni = aconi 
  conm = aconm
  pmlt = apmlt
  trpt = atrpt
end subroutine column_thermodynamics_set_constants

subroutine column_thermodynamics_update_internal_energy(energy, senergy, sflux, sdiric, benergy, mask, &
     bflux, rhs, thckold, thcknew, usig, fsig, dsig, time, dt, n)

  !update the mid-layer internal energy density (energy), 
  !the surface  energy density (senergy), and the basal energy density (benergy)
  !
  !Advective fluxes are computed from the cross-layer contravariant component of the velocity, usig
 
  !
  !There are n layers, with sigma = 0 at the surface, 
  !sigma = sigma(i) at the middle of layer i and sigma = 1
  !at the base of the ice
  ! 
  use column_thermodynamics

  implicit none
  integer :: n
  real(kind=8), dimension(1:n), intent(inout) :: energy,rhs,dsig
  real(kind=8), dimension(1:n+1), intent(inout) :: usig,fsig
  real(kind=8), intent(inout) :: senergy,benergy,bflux,sflux
  real(kind=8), intent(in) :: time,dt,thckold,thcknew
  integer, intent(in) :: mask
  logical, intent(in) :: sdiric
  !locals

  real(kind=8) :: dtcfl,tthcknew,tthckold,tt,melt,bepmp
  real(kind=8), dimension(1:n) :: rhsl,drain,moist,csig
  ! face-centered conductivity, internal energy, pressure-melting-point, additional heat flux
  integer l,nt,it,i
  

  if (thcknew.gt.1.0d0) then
     
     do i = 1,n
        csig(i) = 0.5*(fsig(i+1)+fsig(i))
        if (energy(i).gt.6.0d+5) then
           write(*,*) 'odd, ', energy(i)
        end if
     end do
    
     !work out a stable time step
     dtcfl = dt
     do i = 1,n
        dtcfl = min(dtcfl,(thckold+thcknew)*(fsig(i+1)-fsig(i))/(abs(usig(i+1)) + abs(usig(i))))
     end do
     nt = ceiling(dt/dtcfl)
     dtcfl = dt/dble(nt)
     tt = 0.0d0
     do it = 1,nt

        rhsl = rhs/dble(nt)

        !modify rhs to include interlayer advection across constant sigma faces
        !
        
        bepmp = shci * (trpt - pmlt * 0.5 * (thckold+thcknew) * rhoi * grav)
        benergy = energy(n) + (energy(n)-energy(n-1))/(csig(n)-csig(n-1)) &
             * (fsig(n+1) - csig(n))
        call sigma_advect(rhsl, energy, senergy, benergy, usig , &
             fsig, dtcfl, mask,  n)

        !compute vertical diffusion & update the internal energy
        tthckold = tt/dt * thcknew + (1.0d0-tt/dt) *thckold
        tt = tt + dtcfl
        tthcknew = tt/dt * thcknew + (1.0d0-tt/dt) *thckold

        if (energy(n).gt.6.0d+5) then
           write(*,*) dt, dtcfl,'odd, ', energy(n)
        end if

        call fo_diffusive_advance(energy, senergy, sflux, sdiric, benergy, bflux, rhsl, &
             tthckold, tthcknew,fsig,dtcfl,mask,n)

     end do
  else
     !no ice
     senergy = max(senergy, shci*(trpt - temp_eps))
     if (sdiric) then
        sflux = bflux
     end if
     energy = senergy 
     benergy = senergy
  end if

  return

end subroutine column_thermodynamics_update_internal_energy

subroutine column_thermodynamics_compute_zvel(uz,uzs,husig,ux,uy, divhu, &
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
  husig(n+1) =  - bmb! + gft(n+1)


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

end subroutine column_thermodynamics_compute_zvel
