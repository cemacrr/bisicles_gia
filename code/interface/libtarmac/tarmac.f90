!A collection of subroutines which are entirely ignorant
!of the horizontal mesh

module glunnamed
  implicit none


  ! a bunch of constants copied from glimmer_physcon
  real(kind=8),parameter :: scyr = 31556926.0d0      !*FD Number of seconds in a year (s). Note
                                                 !*FD that this is for a 365.242 day year, and might
                                                 !*FD need changing.
  real(kind=8),parameter :: pi = 3.1415926535897d0   !*FD Value of $\pi$.
  real(kind=8),parameter :: rhoi = 910.0d0           !*FD The density of ice (kg m$^{-3}$)
  real(kind=8),parameter :: rhom = 3300.0d0          !*FD The density of magma(?) (kg m$^{-3}$)
  real(kind=8),parameter :: rhoo = 1028.0d0          !*FD The density of the ocean (kg m$^{-3}$)
  real(kind=8),parameter :: rhow = 1000.0d0          !*FD The density of fresh water (kg m$^{-3}$)
  real(kind=8),parameter :: rhos = 2600.0d0          !*FD The density of solid till (kg m$^{-3}$)
  real(kind=8),parameter :: f = - rhoo / rhoi
  real(kind=8),parameter :: grav = 9.81d0            !*FD The acceleration due to gravity (m s$^{-2}$)
  integer, parameter :: gn = 3                   !*FD The power dependency of Glenn's flow law.
  real(kind=8),parameter :: arrmlh = 1.733d3         !*FD Constant of proportionality in Arrhenius relation
                                                 !*FD in \texttt{patebudd}, for $T^{*}\geq263$K.
                                                 !*FD (Pa$^{-3}$ s$^{-1}$) 
  real(kind=8),parameter :: arrmll = 3.613d-13       !*FD Constant of proportionality in Arrhenius relation
                                                 !*FD in \texttt{patebudd}, for $T^{*}<263$K.
                                                 !*FD (Pa$^{-3}$ s$^{-1}$) 
  real(kind=8),parameter :: gascon = 8.314d0         !*FD The gas ideal constant $R$ (J mol$^{-1}$ K$^{-1}$)
  real(kind=8),parameter :: actenh = 139.0d3         !*FD Activation energy in Glenn's flow law for $T^{*}\geq263$K. (J mol$^{-1}$)
  real(kind=8),parameter :: actenl = 60.0d3          !*FD Activation energy in Glenn's flow law for $T^{*}<263$K. (J mol$^{-1}$)
  real(kind=8),parameter :: shci = 2009.0d0          !*FD Specific heat capacity of ice (J kg$^{-1}$ K$^{-1}$)
  real(kind=8),parameter :: lhci = 335.0d3           !*FD Latent heat of melting of ice (J kg$^{-1}$) 
  real(kind=8),parameter :: coni = 2.1d0             !*FD Thermal conductivity of ice (W m$^{-1}$ K$^{-1}$)
  real(kind=8),parameter :: pmlt = 9.7456d-8         !*FD Factor for dependence of melting point on pressure (K Pa$^{-1}$)
  real(kind=8),parameter :: trpt = 273.15d0          !*FD Triple point of water (K)

  integer, parameter :: groundedmaskval = 1, floatingmaskval = 2, openseamaskval = 4, openlandmaskval = 8

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

    subroutine sigma_advect(rhs, temp, stemp, btemp, usig , lthck, & 
         fsig, dt, mask, n)
      !increment rhs with interlayer advection
      !TODO : replace first order upwind fluxes with something nice
      integer :: n
      real(kind=8), dimension(1:n), intent(inout) :: rhs
      real(kind=8), dimension(1:n), intent(in) :: temp,lthck ! temperature, layer thickness
      real(kind=8), dimension(1:n+1), intent(in) :: usig ! u^{\sigma}
      real(kind=8), dimension(1:n+1), intent(in) :: fsig ! sigma at layer faces
      real(kind=8), intent(in) :: stemp, btemp, dt ! surface & basal temperature, time step
      integer, intent(in) :: mask
      !locals
      integer l
      real(kind=8) flux
      
      !flux across the surface : Dirichlett condition
      flux = dt * (fo_upwind(usig(1) ,stemp,temp(1)) )
      rhs(1) = rhs(1) + flux / (fsig(2)-fsig(1))

      do l = 2, n
         flux = dt * (fo_upwind(usig(l) ,temp(l-1),temp(l)) )
         rhs(l-1) = rhs(l-1) - flux / (fsig(l)-fsig(l-1))
         rhs(l) = rhs(l) + flux / (fsig(l+1)-fsig(l))
      end do

      !flux across the base
      if (mask.eq.floatingmaskval) then
         ! Dirichlett condition
         !flux = dt * (fo_upwind(usig(n+1) ,temp(n),btemp) )
         flux = dt * (fo_upwind(usig(n+1) ,temp(n),temp(n)) )
         
      else
         ! no advection flux
         flux =  0.0
      end if
      rhs(n) = rhs(n) - flux / (fsig(n+1)-fsig(n))

    end subroutine sigma_advect

    subroutine fo_diffusive_advance(temp,stemp,btemp,bflux,rhs,chi, &
         thckold,thcknew,fsig,dt,mask,n)
      ! solve the equation 
      ! H_nT_n - dt * 1/H * \chi * d^2 T_n / dsigma^2 = H_To + rhs
      ! with Robin boundary conditions 
      ! (first order - backward Euler integration)
      integer :: n
      real(kind=8), intent(inout) :: stemp,btemp
      real(kind=8), intent(in) :: bflux,chi,thckold,thcknew,dt
      real(kind=8), dimension(1:n), intent(in) :: rhs
      real(kind=8), dimension(1:n), intent(inout) :: temp
      real(kind=8), dimension(1:n+1), intent(in) :: fsig ! sigma at layer faces
      integer, intent(in) :: mask
      !locals
      real(kind=8), dimension(1:n) :: a,c,b,r ! TDMA coefficients
      real(kind=8), dimension(1:n) :: cdsig
      real(kind=8), dimension(1:n) :: fdsig
      real(kind=8) :: chistar, tpmp, eps, bmb, dtempdz
      logical :: diric
      integer :: i

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

      ! top layer : Dirichlett boundary condition temp = stemp
      a(1) = 0.
      b(1) =  thcknew + 3. * chistar / (cdsig(1)*fdsig(1))
      c(1) = - chistar / (cdsig(1)*fdsig(1))
      r(1) = thckold * temp(1) +  rhs(1) + &
           2. * stemp * chistar / (cdsig(1)*fdsig(1))

      ! bottom layer
      eps = 1.0e-12 ! hold temperature eps K below the pressure melting point
      tpmp = trpt - pmlt * rhoi * grav * thcknew - eps
      !diric = ((mask.eq.floatingmaskval) .or. (btemp .ge. tpmp))
      
      diric = .false.
      if ((btemp .ge. tpmp))  then 
         diric = .true.
         btemp = tpmp 
         !compute basal melt rate
         dtempdz = temp(n) + pmlt * rhoi * grav * thcknew * 0.5*(fsig(i) + fsig(i+1))
         dtempdz = 2.0 / thcknew * (dtempdz - trpt) / cdsig(n)
         bmb = (coni * scyr * dtempdz + bflux * (shci*rhoi))/(rhoi * lhci)
         !write(*,*) coni * scyr * dtempdz, bflux* (shci*rhoi),  bmb
      end if

      if (diric) then
         ! ice at pressure melting point 
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

      if (.not. diric) then 
         btemp = 1.5 * temp(n) - 0.5 * temp(n-1)
         !write(*,*) chistar * (btemp-temp(n)) / ( 0.5* fdsig(n+1)), dt * bflux
         if (btemp .ge. tpmp) btemp = tpmp   
      end if

       do i = 1, n
           tpmp = trpt - pmlt * rhoi * grav * thcknew * 0.5*(fsig(i) + fsig(i+1)) - eps
           if (temp(i) .ge. tpmp) temp(i) = tpmp 
       end do


    end subroutine fo_diffusive_advance

end module glunnamed




subroutine tarmac_update_temperature(temp, stemp, btemp, mask, &
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
  real(kind=8), intent(inout) :: stemp,btemp
  real(kind=8), intent(in) :: time,dt,thckold,thcknew,bflux
  integer, intent(in) :: mask
  !locals
  real(kind=8), dimension(1:n+1) :: lthck
  real(kind=8) :: chi
  
  integer l
 
  if ((mask.eq.groundedmaskval).or.(mask.eq.floatingmaskval)) then
     chi = coni/(shci * rhoi) *  scyr 
     lthck(1:n) = 0.5*(thcknew+thckold) * dsig(1:n)
     lthck(1:n) = dsig(1:n)
     !modify rhs to include interlayer advection across constant sigma faces
     call sigma_advect(rhs, temp, stemp, btemp, usig , &
          lthck, fsig, dt, mask,  n)
     !compute vertical diffusion & update the basal temperature
     !chi = 0.0d0;
     call fo_diffusive_advance(temp, stemp, btemp, bflux, rhs, chi, &
          thckold, thcknew,fsig,dt,mask,n)
  else
     !no ice
     temp = stemp
     btemp = stemp
  end if

end subroutine tarmac_update_temperature

subroutine tarmac_compute_zvel(uz,uzs,usig,ux,uy, divhu, &
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
  real(kind=8), dimension(1:n+1), intent(out) :: uz,usig ! vertical velocity at layer faces
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
  usig(n+1) =  - bmb + gft(n+1)

  !velocity at the surface
  uzs = dst + ux(1)*dsx + uy(1)*dsy - smb

  

  !integration
  do i = n , 1, -1
     uz(i) = uz(i+1) - divhu(i)*dsig(i) & 
          + (gfx(i)*ux(i) - gfx(i+1)*ux(i+1)) &
          + (gfy(i)*uy(i) - gfy(i+1)*uy(i+1))  
     
     !u^{\sigma}
     usig(i) =  usig(i+1) + (divhu(i) + dht)*dsig(i)
     
  end do
  
end subroutine tarmac_compute_zvel
