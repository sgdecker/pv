module diagnostics
  implicit none
  save
  private

  public :: vint, hydro_interp, mix_ratio_to_spec_hum, temp_from_theta_p,  &
            surface_pres, density, sea_pres, lifted_index, mean_rh,  &
            unstag_hght, cape, sbcape, mlcape, lcl_lfc, storm_motion,  &
            storm_rel_hel, calcdbz, theta_phi, pvor
  public :: RGas, Grav, Kappa, Eps

  ! Some meteorological constants (set to the values used by WRF)
  real, parameter :: RGas = 287., Grav = 9.81, KappaR = 3.5,  &
                     SpecHeatPresDry = KappaR * RGas,  &
                     Kappa = RGas / SpecHeatPresDry, Gamma = .0065,  &
                     LatHeatVap = 2.5e6, RVap = 461.6, Eps = RGas / RVap,  &
                     EpsR = RVap / RGas, Missing = -9999., EpsRM1 = EpsR - 1

  interface mix_ratio_to_spec_hum
     module procedure mrtsh3d
     module procedure mrtsh4d
  end interface

  interface temp_from_theta_p
     module procedure tftp0d
     module procedure tftp1d
     module procedure tftp2d
     module procedure tftp3d
     module procedure tftp4d
  end interface

  interface density
     module procedure dens3d
     module procedure dens4d
  end interface

contains  ! Public Functions ==================================================

  ! vint takes a vertical profile (prof) and interpolates it to regular levels
  ! starting at first with an interval of delta using a profile of values for
  ! the interpolating vertical coordinate variable (vCoordVals) defined at the
  ! same points as prof.  The interpolated output is returned as interpProf.
  ! interpType = 1 --> linear interpolation, interpType = 2 --> interpolate
  ! with respect to -ln(vertical coordinate value).
  subroutine vint(prof, vCoordVals, first, delta, interpType, interpProf)
    real, dimension(:), intent(in)  :: prof, vCoordVals
    integer,            intent(in)  :: first, delta, interpType
    real, dimension(:), intent(out) :: interpProf
    
    real, dimension(2) :: vCoordBounds, profBounds
    real :: level
    integer :: nz, np, k

    nz = size(prof)
    np = size(interpProf)

    ! Validate input
    if (nz /= size(vCoordVals)) then
       print *, "Error: Array extent mismatch in vint!"
       stop
    end if

    ! Initialize output
    interpProf = Missing

    ! Interpolate
    do k = 1, np
       level = first + (k-1) * delta
       ! Skip this level if we don't have data for it
       if (level < minval(vCoordVals) .or. level > maxval(vCoordVals)) cycle

       call find_bounds(vCoordVals, prof, level, vCoordBounds, profBounds)

       if (interpType == 2) then  ! Handle logarithmic interpolation case
          vCoordBounds = -log(vCoordBounds)
          level = -log(level)
       end if

       interpProf(k) = lin_int(vCoordBounds(1), vCoordBounds(2),  &
            profBounds(1), profBounds(2), level)
    end do
  end subroutine vint

  ! ===========================================================================

  ! hydro_interp hydrostatically interpolates geopotential heights on pressure
  ! surfaces underground.  Provided as input are the height (zs; m), pressure
  ! (ps; hPa), and temperature (ts; K) in the lowest WRF model layer, as well
  ! as the lowest (in altitude) pressure surface desired (pBot) and the
  ! pressure level spacing (dp).  The routine takes the currently computed
  ! geopotential heights on pressure surfaces (z), and interpolates to those
  ! grid points determined to be underground.
  subroutine hydro_interp(zs, ps, ts, pBot, dp, z)
    real, dimension(:,:,:),   intent(in)    :: zs, ps, ts
    integer,                  intent(in)    :: pBot, dp
    real, dimension(:,:,:,:), intent(inout) :: z

    real, parameter :: LapseRate = .0065, Expo = RGas * LapseRate / Grav
    
    integer :: nz, k, p

    nz = size(z,3)

    ! Validate
    if (any(shape(zs) /= shape(ps)) .or. any(shape(zs) /= shape(ts)) .or.  &
         size(z,1) /= size(zs,1) .or. size(z,2) /= size(zs,2) .or.  &
         size(z,4) /= size(zs,3)) stop "Array mismatch in hydro_interp!"

    do k = 1, nz
       if (all(z /= Missing)) exit  ! Are we above the highest terrain?
       
       p = pBot + (k-1) * dp

       where (z(:,:,k,:) == Missing)
          z(:,:,k,:) = zs - ts / LapseRate * ((p/ps)**Expo - 1)
       end where
    end do
  end subroutine hydro_interp

  ! ===========================================================================

  ! Both mrtsh3d and mrtsh4d convert mixing ratio (kg/kg) to specific
  ! humidity (kg/kg).
  function mrtsh3d(w) result (q)
    real, dimension(:,:,:), intent(in)             :: w
    real, dimension(size(w,1),size(w,2),size(w,3)) :: q

    q = w / (1 + w)
  end function mrtsh3d
    
  ! ...........................................................................

  function mrtsh4d(w) result (q)
    real, dimension(:,:,:,:), intent(in)                     :: w
    real, dimension(size(w,1),size(w,2),size(w,3),size(w,4)) :: q
    
    q = w / (1 + w)
  end function mrtsh4d

  ! ===========================================================================

  ! tftp0d, tftp1d, tftp2d, tftp3d, and tftp4d calculate the temperature (K) 
  ! given the potential temperature (th; K) and pressure (p; Pa or hPa
  ! depending on si flag).
  function tftp0d(th, p, si) result (t)
    real,    intent(in)           :: th, p
    logical, intent(in), optional :: si
    real                          :: t

    real :: p0r
    
    ! Set appropriate p0
    p0r = 1e-5  ! Default is SI units (p0 = 100000)
    if (present(si)) then
       if (.not. si) p0r = .001  ! (p0 = 1000)
    end if
    
    ! Calculate temperature (exp and log faster than ** on my machine)
    t = th * exp(Kappa * log(p * p0r))
  end function tftp0d

  ! ...........................................................................

  function tftp1d(th, p, si) result (t)
    real, dimension(:), intent(in)           :: th, p
    logical,            intent(in), optional :: si
    real, dimension(size(p,1))               :: t

    real :: p0r
    
    ! Make sure arrays match up
    if (any(shape(th) /= shape(p))) stop "Array extent mismatch in tftp2d!"
    
    ! Set appropriate p0
    p0r = 1e-5  ! Default is SI units (p0 = 100000)
    if (present(si)) then
       if (.not. si) p0r = .001  ! (p0 = 1000)
    end if
    
    ! Calculate temperature
    t = th * (p * p0r) ** Kappa
  end function tftp1d

  ! ...........................................................................

  function tftp2d(th, p, si) result (t)
    real, dimension(:,:), intent(in)           :: th, p
    logical,              intent(in), optional :: si
    real, dimension(size(p,1),size(p,2))       :: t

    real :: p0r

    ! Make sure arrays match up
    if (any(shape(th) /= shape(p))) stop "Array extent mismatch in tftp2d!"
    
    ! Set appropriate p0
    p0r = 1e-5  ! Default is SI units (p0 = 100000)
    if (present(si)) then
       if (.not. si) p0r = .001  ! (p0 = 1000)
    end if
    
    ! Calculate temperature
    t = th * (p * p0r) ** Kappa
  end function tftp2d
       
  ! ...........................................................................

  function tftp3d(th, p, si) result (t)
    real, dimension(:,:,:), intent(in)             :: th, p
    logical,                intent(in), optional   :: si
    real, dimension(size(p,1),size(p,2),size(p,3)) :: t
    
    real :: p0r
    
    ! Make sure arrays match up
    if (any(shape(th) /= shape(p))) stop "Array extent mismatch in tftp3d!"
    
    ! Set appropriate p0
    p0r = 1e-5  ! Default is SI units (p0 = 100000)
    if (present(si)) then
       if (.not. si) p0r = .001  ! (p0 = 1000)
    end if
    
    ! Calculate temperature
    t = th * (p * p0r) ** Kappa
  end function tftp3d

  ! ...........................................................................

  function tftp4d(th, p, si) result (t)
    real, dimension(:,:,:,:), intent(in)                     :: th, p
    logical,                  intent(in), optional           :: si
    real, dimension(size(p,1),size(p,2),size(p,3),size(p,4)) :: t
    
    real :: p0r
    
    ! Make sure arrays match up
    if (any(shape(th) /= shape(p))) stop "Array extent mismatch in tftp4d!"
    
    ! Set appropriate p0
    p0r = 1e-5  ! Default is SI units (p0 = 100000)
    if (present(si)) then
       if (.not. si) p0r = .001  ! (p0 = 1000)
    end if
    
    ! Calculate temperature
    t = th * (p * p0r) ** Kappa
  end function tftp4d
  
  ! ===========================================================================

  ! surface_pres calculates the surface pressure (Pa) given the following data:
  ! ht-   Geopotential height (m) at lowest two w (full) levels  (x,y,z,t)
  ! p1-   Pressure (Pa) at lowest mass (half) level              (x,y,t)
  ! q1-   Specific humidity (kg/kg) at lowest mass level         (x,y,t)
  ! t1-   Temperature (K) at lowest mass level                   (x,y,t)
  ! sfcq- Specific humidity (kq/kq) at surface/2-m above ground  (x,y,t)
  ! It is assumed that the array extents match appropriately.
  function surface_pres(ht, p1, q1, t1, sfcq)
    real, dimension(:,:,:,:), intent(in)              :: ht
    real, dimension(:,:,:), intent(in)                :: p1, q1, t1, sfcq
    real, dimension(size(p1,1),size(p1,2),size(p1,3)) :: surface_pres

    real, dimension(size(t1,1),size(t1,2),size(t1,3)) :: tvirt

    tvirt = (1 + .608 * .5 * (sfcq + q1)) * t1

    surface_pres = p1 * exp((Grav / (RGas * tvirt))  &
         * .5 * (ht(:,:,2,:) - ht(:,:,1,:)))
  end function surface_pres

  ! ===========================================================================

  ! dens3d and dens4d calculate the density (kg/m**3) using the ideal gas law
  ! given:
  ! p- Pressure (Pa)
  ! t- Temperature (K)
  ! q- Water vapor mixing ratio (kg/kg)
  ! It is assumed that the array extents match appropriately.
  function dens3d(p, t, q)
    real, dimension(:,:,:), intent(in)             :: p, t, q
    real, dimension(size(p,1),size(p,2),size(p,3)) :: dens3d

    dens3d = p / (t * (RGas + q * RVap))
    dens3d = dens3d * (1 + q)
  end function dens3d

  ! ...........................................................................

  function dens4d(p, t, q)
    real, dimension(:,:,:,:), intent(in)                     :: p, t, q
    real, dimension(size(p,1),size(p,2),size(p,3),size(p,4)) :: dens4d

    dens4d = p / (t * (RGas + q * RVap))
    dens4d = dens4d * (1 + q)
  end function dens4d

  ! ===========================================================================

  ! sea_pres is taken from the WRF NCL tools.  The inputs are:
  ! zs-    Geopotential height (m) on w (full) levels       (x,y,z,t)
  ! theta- Potential temperature (K) on mass (half) levels  (x,y,z,t)
  ! p-     Pressure (Pa) on mass levels                     (x,y,z,t)
  ! q-     Specific humidity (kg/kg) on mass levels         (x,y,z,t)
  ! It is assumed that the array extents match appropriately.
  function sea_pres(zs, theta, p, q)
    ! Estimate sea level pressure.
    real, dimension(:,:,:,:), intent(in)           :: zs, theta, p, q
    real, dimension(size(p,1),size(p,2),size(p,4)) :: sea_pres

    real,    dimension(size(p,1),size(p,2),size(p,3),size(p,4)) :: t, z
    real,    dimension(size(p,1),size(p,2)) :: tSurf, tSeaLevel
    integer, dimension(size(p,1),size(p,2)) :: level
    
    ! Specific constants for assumptions made in this routine:
    real, parameter :: Tc = 273.16+17.5, PConst = 10000
    
    !  Local variables:
    real    :: plo, phi, tlo, thi, zlo, zhi, pAtPConst, tAtPConst, zAtPConst
    integer :: nx, ny, nz, nt
    integer :: i, j, k, klo, khi, l
    logical :: l1, l2, l3, found
    
    nx = size(p,1)
    ny = size(p,2)
    nz = size(p,3)
    nt = size(p,4)
    
    ! Convert potential temperature to temperature
    t = temp_from_theta_p(theta, p)
    
    ! Unstagger z
    z(:,:,1:nz,:) = .5 * (zs(:,:,1:nz,:) + zs(:,:,2:nz+1,:))
    
    ! Grand loop over all times
    do l = 1, nt
       ! Find least zeta level that is PConst Pa above the surface.  We later
       ! use this level to extrapolate a surface pressure and temperature,
       ! which is supposed to reduce the effect of the diurnal heating cycle in
       ! the pressure field.
       
       do j = 1 , ny
          do i = 1 , nx
             level(i,j) = -1
             
             k = 1
             found = .false.
             do while( (.not. found) .and. (k <= nz))
                if ( p(i,j,k,l) < p(i,j,1,l)-PConst ) then
                   level(i,j) = k
                   found = .true.
                end if
                k = k+1
             end do
             
             if ( level(i,j) == -1 ) then
                print '(A,I4,A)','Troubles finding level ',  &
                     nint(PConst)/100,' above ground.'
                print '(A,I4,A,I4,A)',  &
                     'Problems first occur at (',i,',',j,')'
                print '(A,F6.1,A)',  &
                     'Surface pressure = ',p(i,j,1,l)/100,' hPa.'
                stop 'Error_in_finding_100_hPa_up'
             end if
                          
          end do
       end do
       
       ! Get temperature PConst Pa above surface.  Use this to extrapolate 
       ! the temperature at the surface and down to sea level.
       
       do j = 1 , ny
          do i = 1 , nx
             
             klo = max ( level(i,j) - 1 , 1      )
             khi = min ( klo + 1        , nz - 1 )
             
             if ( klo == khi ) then
                print '(A)','Trapping levels are weird.'
                print '(A,I3,A,I3,A)','klo = ',klo,', khi = ',khi,  &
                     ': and they should not be equal.'
                stop 'Error_trapping_levels'
             end if
             
             plo = p(i,j,klo,l)
             phi = p(i,j,khi,l)
             tlo = t(i,j,klo,l) * (1. + 0.608 * q(i,j,klo,l) )
             thi = t(i,j,khi,l) * (1. + 0.608 * q(i,j,khi,l) )
             zlo = z(i,j,klo,l)         
             zhi = z(i,j,khi,l)
             
             pAtPConst = p(i,j,1,l) - PConst
             TAtPConst = thi - (thi-tlo) * log(pAtPConst/phi) * log(plo/phi)
             zAtPConst = zhi - (zhi-zlo) * log(pAtPConst/phi) * log(plo/phi)
             
             tSurf(i,j) = tAtPConst * (p(i,j,1,l)/pAtPConst)**(Gamma*RGas/Grav)
             tSeaLevel(i,j) = tAtPConst + Gamma * zAtPConst
             
          end do
       end do

       ! If we follow a traditional computation, there is a correction to the
       ! sea level temperature if both the surface and sea level temnperatures
       ! are *too* hot.
       
       do j = 1 , ny
          do i = 1 , nx
             l1 = tSeaLevel(i,j) <  Tc 
             l2 = tSurf    (i,j) <= Tc
             l3 = .not. l1
             if ( l2 .and. l3 ) then
                tSeaLevel(i,j) = Tc
             else
                tSeaLevel(i,j) = Tc - 0.005 * (tSurf(i,j) - Tc)**2
             end if
          end do
       end do
       
       ! The grand finale: ta da!
       sea_pres(:,:,l) = p(:,:,1,l) *  &
            exp((2.*Grav*z(:,:,1,l)) / (RGas*(tSeaLevel + tSurf)))

    end do
  end function sea_pres

  ! ===========================================================================

  ! lifted_index calculates the lifted index (K) based on the supplied profiles
  ! of potential temperature (theta; K), pressure (p; hPa), and mixing ratio
  ! (w; kg/kg).  Array indices are ordered from the ground up.  The lowest 4
  ! model levels are mixed to produce the hypothetical parcel, which is then
  ! lifted adiabatically to saturation, and pseudoadiabatically thereafter.
  function lifted_index(theta, p, w) result (li)
    real, dimension(:), intent(in) :: theta, p, w
    real                           :: li

    real, dimension(2) :: pBounds, thBounds
    real               :: dens, thMix, wMix, sm, tk, pLCL, thes, tParcel, tEnv
    integer            :: nz, k
    
    nz = size(theta)
    
    ! Check that array extents match
    if (nz /= size(p) .or. nz /= size(w))  &
         stop "Array mismatch in lifted_index."
    
    ! Set lifted index to 0 if we don't have enough data
    if (nz < 4 .or. p(1) < 500.) then
       li = 0
    else
       ! Determine temperature of pseudoadiabatically lifted parcel at 500 mb
       ! Step 1: Mix the lowest 4 model levels
       thMix = 0; wMix = 0; sm = 0
       do k = 1, 4
          dens = density_thp(theta(k), p(k))  ! Effect of water vapor ignored
          thMix = thMix + dens * theta(k)
          wMix = wMix + dens * w(k)
          sm = sm + dens
       end do
       thMix = thMix / sm
       wMix = wMix / sm

       ! Step 2: Lift parcel to LCL
       tk = temp_from_theta_p(thMix, p(1), si=.false.)
       ! This equation is derived from Curry & Webster (1999, p. 174)
       pLCL = 1000 * (temp_lcl(tk, rh(wMix, p(1), tk)) / thMix)**KappaR
       
       ! Step 3: Lift parcel to 500 mb
       thes = theta_es(thMix, pLCL)
       tParcel = temp_from_theta_p(theta_es_to_theta_fast(thes, 500.), 500.,  &
            si=.false.)

       ! Get environmental temperature at 500 mb
       call find_bounds(p, theta, 500., pBounds, thBounds)
       tEnv = temp_from_theta_p(lin_int(pBounds(1), pBounds(2), thBounds(1),  &
            thBounds(2), 500.), 500., si=.false.)
       
       li = tEnv - tParcel
    end if
  end function lifted_index

  ! ===========================================================================

  ! mean_rh computes the mean relative humidity (%) between pBot hPa and pTop
  ! hPa based on the supplied profiles of potential temperature (theta; K),
  ! pressure (p; hPa), and mixing ratio (w; kg/kg).  Also provided is the
  ! surface pressure (pSfc; hPa).  Array indices are ordered from the ground
  ! up.  The mean is weighted by the thickness of each layer.
  function mean_rh(theta, p, w, pSfc, pBot, pTop) result (mrh)
    real, dimension(:), intent(in) :: theta, p, w
    real,               intent(in) :: pSfc, pBot, pTop
    real                           :: mrh

    real, dimension(0:size(theta)-1) :: pMid
    real                             :: pThick, relHum, totPresThick, totRH
    integer                          :: nz, k

    nz = size(theta)

    ! Check that array extents match
    if (nz /= size(p) .or. nz /= size(w)) stop "Array mismatch in mean_rh."
    
    ! Calculate pressure at midpoints
    pMid(0) = pSfc
    pMid(1:nz-1) = .5 * (p(1:nz-1) + p(2:nz))

    ! Calculate mean RH
    totPresThick = 0; totRH = 0
    do k = 1, nz-1
       if (pMid(k) > pBot) cycle
       if (pMid(k-1) < pTop) exit

       ! Calculate layer thickness
       if (pMid(k-1) > pBot) then
          pThick = pBot - pMid(k)
       else if (pMid(k) < pTop) then
          pThick = pMid(k-1) - pTop
       else
          pThick = pMid(k-1) - pMid(k)
       end if

       relHum = rh(w(k), p(k), temp_from_theta_p(theta(k), p(k), si=.false.))
       totRH = totRH + pThick * relHum
       totPresThick = totPresThick + pThick
    end do
    
    if (totPresThick > .01) then
       mrh = 100 * totRH / totPresThick
    else
       mrh = 0
    end if
  end function mean_rh

  ! ===========================================================================

  ! unstag_hght unstaggers the geopotential height from w levels (zs) given
  ! by znw to u levels (z) given by znu.  mu and pt provide additional
  ! information needed to do an interpolation based on ln p.
  function unstag_hght(zs, znw, znu, mu, pt) result (z)
    real, dimension(:), intent(in) :: zs, znw, znu
    real,               intent(in) :: mu, pt
    real, dimension(size(znu))     :: z

    real, dimension(size(znw)) :: pw
    real, dimension(size(znu)) :: pu
    integer :: k

    ! Validate input
    if (size(zs) /= size(znw) .or. size(znw) /= size(znu) + 1)  &
         print *, "Warning: Array extent mismatch in unstag_hght!"

    ! Initialize pressures
    do k = 1, size(znw)
       if (k < size(znw)) pu(k) = pd(znu(k))
       pw(k) = pd(znw(k))
    end do
    
    ! Interpolate to half levels
    do k = 1, size(znu)
!orig   z(k) = lin_int(-log(pw(k)), -log(pw(k+1)), zs(k), zs(k+1), -log(pu(k)))
       z(k) = lin_int(pw(k), pw(k+1), zs(k), zs(k+1), pu(k))
    end do

  contains
    
    real function pd(eta)
      real, intent(in) :: eta
      
      pd = pt + mu * eta
    end function pd
  end function unstag_hght

  ! ===========================================================================

  ! cape calculates the convective available potential energy (J/kg) given
  ! profiles of potential temperature (theta; K), pressure (p; hPa), 
  ! mixing ratio (w; kg/kg), and height (z; m).  By default, it returns the
  ! most unstable CAPE, but if lplOutI is set to .true., the height above
  ! ground level of the parcel whose CAPE is most unstable is returned instead.
  ! In that case, the terrain height at the profile location (ter) is also
  ! required.
  real function cape(theta, p, w, z, ter, lplOutI)
    real, dimension(:), intent(in) :: theta, p, w, z
    real,    optional,  intent(in) :: ter
    logical, optional,  intent(in) :: lplOutI

    real, dimension(size(z)) :: tProf
    real, dimension(2)       :: pBounds, zBounds
    real    :: lpl, tk, pLCL, thes, zLCL, parcelCape
    integer :: nz, i
    logical :: lplOut
    
    ! Initialize
    cape = 0; lpl = 0
    nz = size(z)
    lplOut = .false.
    if (present(lplOutI)) then
       lplOut = lplOutI
       if (.not. present(ter)) stop "For LPL, need terrain height!"
    end if

    ! Validate
    if (nz /= size(theta) .or. nz /= size(p) .or. nz /= size(w))  &
         stop "Error: Array extent mismatch in cape!"
    
    ! Compute environmental profile with virtual temperature correction.
    ! Note that applying the correction here doesn't affect LCL values below.
    tProf = temp_from_theta_p(theta, p, si=.false.) * (1 + EpsRM1 * w)

    ! We loop over all parcels originating in the lowest 300 hPa
    do i = 1, nz-1
       if (p(i) < p(1) - 300) exit
       
       ! Find LCL for this parcel
       tk = temp_from_theta_p(theta(i), p(i), si=.false.)
       pLCL = 1000 * (temp_lcl(tk, rh(w(i), p(i), tk)) / theta(i))**KappaR
       if (pLCL > p(i)) pLCL = p(i)
       thes = theta_es(theta(i), pLCL)
       call find_bounds(p, z, pLCL, pBounds, zBounds)
       zLCL = lin_int(pBounds(1), pBounds(2), zBounds(1), zBounds(2), pLCL)

       parcelCape = calc_cape(tProf(i:nz), p(i:nz), z(i:nz), theta(i), w(i),  &
            thes, pLCL, zLCL)

       ! Check if largest CAPE so far
       if (parcelCape > cape) then
          cape = parcelCape
          lpl = z(i) - ter
       end if
    end do
    if (lplOut) cape = lpl
  end function cape
  
  ! ===========================================================================

  ! sbcape calculates the convective available potential energy (J/kg) of the
  ! parcel lifted from the lowest model layer given profiles of potential
  ! temperature (theta; K), pressure (p; hPa), mixing ratio (w; kg/kg), and
  ! height (z; m).  By default, it returns the CAPE, but if cinOutI is set to
  ! .true., the convective inhibition (J/kg) is returned instead.
  real function sbcape(theta, p, w, z, cinOutI)
    real, dimension(:), intent(in) :: theta, p, w, z
    logical, optional,  intent(in) :: cinOutI
    
    real, dimension(size(z)) :: tProf
    real, dimension(2)       :: pBounds, zBounds
    real    :: tk, pLCL, zLCL, thes
    integer :: nz
    logical :: cinOut
    
    ! Initialize
    sbcape = 0
    nz = size(z)
    cinOut = .false.
    if (present(cinOutI)) cinOut = cinOutI

    ! Validate
    if (nz /= size(theta) .or. nz /= size(p) .or. nz /= size(w))  &
         stop "Error: Array extent mismatch in sbcape!"
    
    ! Compute environmental profile with virtual temperature correction.
    ! Note that applying the correction here doesn't affect LCL values below.
    tProf = temp_from_theta_p(theta, p, si=.false.) * (1 + EpsRM1 * w)

    ! Find LCL for lowest parcel
    tk = temp_from_theta_p(theta(1), p(1), si=.false.)
    pLCL = 1000 * (temp_lcl(tk, rh(w(1), p(1), tk)) / theta(1))**KappaR
    if (pLCL > p(1)) pLCL = p(1)
    thes = theta_es(theta(1), pLCL)
       
    ! Calculate what we want
    if (cinOut) then
       sbcape = calc_cin(tProf, p, z, theta(1), w(1), thes, pLCL)
    else
       call find_bounds(p, z, pLCL, pBounds, zBounds)
       zLCL = lin_int(pBounds(1), pBounds(2), zBounds(1), zBounds(2), pLCL)
       sbcape = calc_cape(tProf, p, z, theta(1), w(1), thes, pLCL, zLCL)
    end if
  end function sbcape

  ! ===========================================================================

  ! mlcape calculates the convective available potential energy (J/kg) of the
  ! parcel mixed from the lowest 4 model layers given profiles of potential
  ! temperature (theta; K), pressure (p; hPa), mixing ratio (w; kg/kg), and
  ! height (z; m).  By default, it returns the CAPE, but if cinOutI is set to
  ! .true., the convective inhibition (J/kg) is returned instead.
  real function mlcape(theta, p, w, z, cinOutI)
    real, dimension(:), intent(in) :: theta, p, w, z
    logical, optional,  intent(in) :: cinOutI
    
    real, dimension(size(z)) :: tProf
    real, dimension(2)       :: pBounds, zBounds
    real    :: thMix, wMix, total, dens, tk, pLCL, zLCL, thes
    integer :: nz, k
    logical :: cinOut
    
    ! Initialize
    mlcape = 0
    nz = size(z)
    cinOut = .false.
    if (present(cinOutI)) cinOut = cinOutI

    ! Validate
    if (nz /= size(theta) .or. nz /= size(p) .or. nz /= size(w))  &
         stop "Error: Array extent mismatch in mlcape!"

    ! Set result to 0 if we don't have enough data
    if (nz < 4) then
       mlcape = 0
       return
    end if

    ! Mix the lowest 4 model levels
    thMix = 0; wMix = 0; total = 0
    do k = 1, 4
       dens = density_thp(theta(k), p(k))  ! Effect of water vapor ignored
       thMix = thMix + dens * theta(k)
       wMix = wMix + dens * w(k)
       total = total + dens
    end do
    thMix = thMix / total
    wMix = wMix / total
    
    ! Compute environmental profile with virtual temperature correction.
    ! Note that applying the correction here doesn't affect LCL values below.
    tProf = temp_from_theta_p(theta, p, si=.false.) * (1 + EpsRM1 * w)

    ! Find LCL for mixed parcel
    tk = temp_from_theta_p(thMix, p(1), si=.false.)
    pLCL = 1000 * (temp_lcl(tk, rh(wMix, p(1), tk)) / thMix)**KappaR
    if (pLCL > p(1)) pLCL = p(1)
    thes = theta_es(thMix, pLCL)

    ! Calculate what we want
    if (cinOut) then
       mlcape = calc_cin(tProf, p, z, thMix, wMix, thes, pLCL)
    else
       call find_bounds(p, z, pLCL, pBounds, zBounds)
       zLCL = lin_int(pBounds(1), pBounds(2), zBounds(1), zBounds(2), pLCL) 
       mlcape = calc_cape(tProf, p, z, thMix, wMix, thes, pLCL, zLCL)
    end if
  end function mlcape

  ! ===========================================================================

  ! lcl_lfc returns the lifting condensation level or level of free convection.
  ! If lclFlag is true (false), the former (latter) is returned.  Provided as
  ! input are profiles of potential temperature (theta; K), pressure (p; hPa),
  ! mixing ratio (w; kg/kg), and geopotential height (z; m).  The terrain
  ! height (ter; m) is also provided as input.
  real function lcl_lfc(theta, p, w, z, ter, lclFlag)
    real, dimension(:), intent(in) :: theta, p, w, z
    real,               intent(in) :: ter
    logical,            intent(in) :: lclFlag
    
    real, dimension(size(z)) :: tProf, tTraj
    real, dimension(2)       :: pBounds, zBounds
    real    :: thMix, wMix, total, dens, tk, pLCL, zLCL, thes
    integer :: nz, k
    
    ! Initialize
    lcl_lfc = -9999  ! The GEMPAK missing value
    nz = size(z)
    
    ! Validate
    if (nz /= size(theta) .or. nz /= size(p) .or. nz /= size(w))  &
         stop "Error: Array extent mismatch in lcl_lfc!"
    
    ! Return immediately if we don't have enough data
    if (nz < 4) return
    
    ! Mix the lowest 4 model levels
    thMix = 0; wMix = 0; total = 0
    do k = 1, 4
       dens = density_thp(theta(k), p(k))  ! Effect of water vapor ignored
       thMix = thMix + dens * theta(k)
       wMix = wMix + dens * w(k)
       total = total + dens
    end do
    thMix = thMix / total
    wMix = wMix / total
    
    ! Find LCL for mixed parcel
    tk = temp_from_theta_p(thMix, p(1), si=.false.)
    pLCL = 1000 * (temp_lcl(tk, rh(wMix, p(1), tk)) / thMix)**KappaR
    if (pLCL > p(1)) pLCL = p(1)
    call find_bounds(p, z, pLCL, pBounds, zBounds)
    zLCL = lin_int(pBounds(1), pBounds(2), zBounds(1), zBounds(2), pLCL) 
    
    ! Calculate what we want
    if (lclFlag) then
       lcl_lfc = zLCL - ter
    else
       thes = theta_es(thMix, pLCL)    

       ! Compute environmental profile and ascent trajectory using the
       ! virtual temperature correction.
       tProf = temp_from_theta_p(theta, p, si=.false.) * (1 + EpsRM1 * w)
       tTraj = calc_traj(p, tProf(1), theta(1), w(1), thes, pLCL)
       
       ! Find LFC height (from bottom up)
       do k = 2, nz-1
          if (p(k) >= pLCL .or. tTraj(k) <= tProf(k)) cycle
          
          ! We just passed through LFC, but we need to allow for case where the
          ! LCL is in a region of positive area already due to a superadiabatic
          ! environmental lapse rate.
          if (tTraj(k-1) - tProf(k-1) >= 0) then
             lcl_lfc = zLCL
          else
             lcl_lfc = lin_int(tTraj(k-1)-tProf(k-1), tTraj(k)-tProf(k),  &
                  z(k-1), z(k), 0.)
          end if
          exit
       end do
       if (lcl_lfc /= -9999) lcl_lfc = lcl_lfc - ter
    end if
  end function lcl_lfc

  ! ===========================================================================

  ! storm_motion computes either the u or v component of the supercell storm
  ! motion vector derived by Bunkers et al. (2000).  uComp = true (false)
  ! means the u (v) component will be output.  Provided for input are the u and
  ! v components of the wind (m), the 10-m winds (u and v; m), and the
  ! terrain height (ter; m).
  function storm_motion(u, v, z, u10m, v10m, ter, uComp)
    real, dimension(:,:,:,:), intent(in) :: u, v, z
    real, dimension(:,:,:),   intent(in) :: u10m, v10m
    real, dimension(:,:),     intent(in) :: ter
    logical,                  intent(in) :: uComp
    real, dimension(size(u,1),size(u,2),size(u,4)) :: storm_motion

    real, parameter :: D = 7.5
    integer, parameter :: NumLevs = 13, Dz = 500  ! 0 to 6 km
    real, dimension(NumLevs), parameter :: Heights = (/ (Dz * k,  &
                                                        k = 0, NumLevs-1) /)
    real, dimension(3), parameter :: HeadHeights = (/ 5500, 5750, 6000 /),  &
                                     TailHeights = (/ 0, 250, 500 /)

    real, dimension(size(u,1),size(u,2),size(u,4)) :: uShear, vShear, magShear
    integer :: k, nx, ny, nt

    nx = size(u,1); ny = size(u,2); nt = size(u,4)

    ! Validate
    if (any(shape(u) /= shape(z)) .or. any(shape(v) /= shape(z)) .or.  &
         nx /= size(u10m,1) .or. ny /= size(u10m,2) .or.  &
         nt /= size(u10m,3) .or. nx /= size(v10m,1) .or.  &
         ny /= size(v10m,2) .or. nt /= size(v10m,3) .or.  &
         nx /= size(ter,1) .or. ny /= size(ter,2))  &
         stop "Array mismatch in storm_motion!"

    ! Calculate Bunkers wind shear and its magnitude
    uShear = mean_wind(u, z, u10m, ter, HeadHeights) -  &
         mean_wind(u, z, u10m, ter, TailHeights)
    vShear = mean_wind(v, z, v10m, ter, HeadHeights) -  &
         mean_wind(v, z, v10m, ter, TailHeights)
    magShear = sqrt(uShear**2 + vShear**2)

    ! Calculate appropriate storm motion vector.
    if (uComp) then
       storm_motion = mean_wind(u, z, u10m, ter, Heights) +  &
            D * VShear / MagShear
    else
       storm_motion = mean_wind(v, z, v10m, ter, Heights) -  &
            D * UShear / MagShear
    end if
  end function storm_motion

  ! ===========================================================================

  ! storm_rel_hel computes the storm relative helicity (m2/s2) for the layer
  ! defined by heights (m), given as input the u and v components of the wind
  ! (m), the geopotential height (z; m), the u and v components of the wind at
  ! 10m (u10m, v10m; m), and the terrain height (ter; m).
  function storm_rel_hel(u, v, z, u10m, v10m, ter, heights) result (srh)
    real, dimension(:,:,:,:), intent(in) :: u, v, z
    real, dimension(:,:,:),   intent(in) :: u10m, v10m
    real, dimension(:,:),     intent(in) :: ter
    real, dimension(:),       intent(in) :: heights
    real, dimension(size(u,1),size(u,2),size(u,4)) :: srh

    real, dimension(size(z,1),size(z,2),size(z,3),size(z,4))       :: hagl
    real, dimension(size(u,1),size(u,2),size(heights),size(u,4))   :: uInt,  &
                                                                      vInt,  &
                                                                      xVort,  &
                                                                      yVort,  &
                                                                      srhs
    real, dimension(size(u,1),size(u,2),size(u,4)) :: uStorm, vStorm
    real, dimension(2) :: haglBounds, windBounds
    integer :: nx, ny, nz, nt, numLevs, i, j, k, t
    
    nx = size(u,1); ny = size(u,2); nz = size(u,3); nt = size(u,4)
    numLevs = size(heights)

    ! Validate
    if (any(shape(u) /= shape(z)) .or. any(shape(v) /= shape(z)) .or.  &
         nx /= size(u10m,1) .or. ny /= size(u10m,2) .or.  &
         nt /= size(u10m,3) .or. nx /= size(v10m,1) .or.  &
         ny /= size(v10m,2) .or. nt /= size(v10m,3) .or.  &
         nx /= size(ter,1) .or. ny /= size(ter,2))  &
         stop "Array mismatch in storm_motion!"
    
    hagl = z -  &
        spread(source=spread(source=ter, dim=3, ncopies=nz), dim=4, ncopies=nt)

    ! Interpolate winds to constant height grid
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             do k = 1, numLevs
                if (heights(k) < hagl(i,j,1,t)) then
                   uInt(i,j,k,t) = u10m(i,j,t)
                   vInt(i,j,k,t) = v10m(i,j,t)
                else
                   call find_bounds(hagl(i,j,:,t), u(i,j,:,t), heights(k),  &
                        haglBounds, windBounds)
                   uInt(i,j,k,t) = lin_int(haglBounds(1), haglBounds(2),  &
                        windBounds(1), windBounds(2), heights(k))
                   call find_bounds(hagl(i,j,:,t), v(i,j,:,t), heights(k),  &
                        haglBounds, windBounds)
                   vInt(i,j,k,t) = lin_int(haglBounds(1), haglBounds(2),  &
                        windBounds(1), windBounds(2), heights(k))
                end if
             end do
          end do
       end do
    end do
    
    ! Neglect vertical motion component of horizontal vorticity, a valid
    ! approximation on the synoptic and mesoscales.
    xVort = -calc_ddz(vInt, heights)
    yVort = calc_ddz(uInt, heights)

    uStorm = storm_motion(u, v, z, u10m, v10m, ter, .true.)
    vStorm = storm_motion(u, v, z, u10m, v10m, ter, .false.)
    
    srhs = (u - spread(uStorm, dim=3, ncopies=numLevs)) * xVort +  &
         (v - spread(vStorm, dim=3, ncopies=numLevs)) * yVort

    srh = 0
    do k = 1, numLevs-1
       srh = srh + .5*(heights(k+1)-heights(k))*(srhs(:,:,k,:)+srhs(:,:,k+1,:))
    end do
  end function storm_rel_hel

  ! ===========================================================================

  ! calcdbz diagnoses reflectivity (dBZ) provided inputs of air density 
  ! (rhoair; kg/m^3), temperature (tmk; K), rain water mixing ratio
  ! (qra; kg/kg), snow mixing ratio (qsn; kg/kg), and graupel mixing ratio
  ! (qgr; kg/kg).  This function is largely taken from that in RIP4.  Note that
  ! to be more accurate the parameters in this function should depend on the
  ! microphysical scheme being used.
  function calcdbz(rhoair, tmk, qra, qsn, qgr) result (dbz)
    real, dimension(:,:,:,:), intent(in) :: rhoair, tmk, qra, qsn, qgr
    real, dimension(size(tmk,1),size(tmk,2),size(tmk,3),size(tmk,4)) :: dbz
    
    real, parameter :: r1 = 1.e-15, ron2 = 1.e10, gon = 5.e7, ron_min = 8.e6, &
                       ron_qr0 = 0.00010, ron_delqr0 = 0.25*ron_qr0,  &
                       ron_const1r = (ron2 - ron_min)*0.5,  &
                       ron_const2r = (ron2 + ron_min)*0.5,  &
                       gamma_seven = 720., rhowat = 1000.,  &
                       rho_r = rhowat,  &  ! 1000. kg m^-3
                       rho_s = 100.,  &    ! kg m^-3
                       rho_g = 400.,  &    ! kg m^-3
                       alpha = 0.224, celkel = 273.15

    real, dimension(size(tmk,1),size(tmk,2),size(tmk,3),size(tmk,4)) :: ronv, &
                                                                        sonv, &
                                                                        gonv
    real :: pi, factor_r, factor_s, factor_g

    ! Constants
    pi = 4. * atan(1.)
    factor_r = gamma_seven * 1.e18 * (1./(pi*rho_r))**1.75
    factor_s = gamma_seven * 1.e18 * (1./(pi*rho_s))**1.75  &
         * (rho_s/rhowat)**2 * alpha
    factor_g = gamma_seven * 1.e18 * (1./(pi*rho_g))**1.75  &
         * (rho_g/rhowat)**2 * alpha
    
    ! Validate
    if (any(shape(tmk) /= shape(rhoair)) .or. any(shape(tmk) /= shape(qra))  &
         .or. any(shape(tmk) /= shape(qsn))  &
         .or. any(shape(tmk) /= shape(qgr))) stop "Array mismatch in calcdbz!"

    ! Calculate variable intercept parameters
    sonv = min(2.0e8, 2.0e6 * exp(-0.12 * min(-0.001, tmk-celkel)))

    where (qgr <= r1)
       gonv = gon
    elsewhere
       gonv = max(1.e4, min(2.38 * (pi*rho_g / (rhoair*qgr))**0.92, gon))
    end where

    where (qra <= r1)
       ronv = ron2
    elsewhere
       ronv = ron_const1r * tanh((ron_qr0-qra) / ron_delqr0) + ron_const2r
    end where

    ! Total equivalent reflectivity factor (in mm^6 mm^-3) is the sum of the
    ! reflectivity for each hydrometeor species (rain, snow, and graupel):
    where (tmk > celkel)
       dbz = factor_r * (rhoair*qra)**1.75 / ronv**.75  &
            + (factor_s * (rhoair*qsn)**1.75 / sonv**.75  &
            + factor_g * (rhoair*qgr)**1.75 / gonv**.75) / alpha
    elsewhere
       dbz = factor_r * (rhoair*qra)**1.75 / ronv**.75  &
            + factor_s * (rhoair*qsn)**1.75 / sonv**.75  &
            + factor_g * (rhoair*qgr)**1.75 / gonv**.75
    end where

    ! Convert reflectivity factor to dBZ while ensuring the result is > 0.
    dbz = 10. * log10(max(dbz, 1.))
  end function calcdbz

  ! ===========================================================================
  
  subroutine theta_phi(phi, mu, etaz, etam, pt, th)
    real, dimension(:,:,0:,:), intent(in) :: phi
    real, dimension(:,:,:),    intent(in) :: mu
    real, dimension(0:),       intent(in) :: etaz
    real, dimension(:),        intent(in) :: etam
    real,                      intent(in) :: pt
    real, dimension(size(phi,1),size(phi,2),size(etam),size(phi,4)) :: th

    real, parameter :: p0 = 100000.

    integer :: nx, ny, nz, nt, i, j, k, t
    
    nx = size(phi,1)
    ny = size(phi,2)
    nz = size(etam)
    nt = size(phi,4)

    forall (i = 1:nx, j = 1:ny, k = 1:nz, t = 1:nt)
       th(i,j,k,t) = (phi(i,j,k,t) - phi(i,j,k-1,t)) * (p0 / (pt + etam(k) * mu(i,j,t)))**Kappa  &
            / (RGas * log((pt + etaz(k-1) * mu(i,j,t)) / (pt + etaz(k) * mu(i,j,t))))
    end forall
  end subroutine theta_phi

  ! ===========================================================================
  
  subroutine pvor(u, v, hght, theta, rho, f, mf, etau, etaw, dx, q)
    real, dimension(:,:,:,:), intent(in) :: u, v, hght, theta, rho
    real, dimension(:,:),     intent(in) :: f, mf
    real, dimension(:),       intent(in) :: etau, etaw
    real,                     intent(in) :: dx
    real, dimension(size(v,1),size(u,2),size(u,3),size(u,4)), intent(out) :: q

    real, dimension(:,:,:,:), allocatable :: dthetadeta, dzdx, dzdy, dvdeta, dudeta, vort, detadz
    real, dimension(3) :: c
    real :: dy, detal, detah, mfavg
    real :: dvdeavg, dudeavg, dzdxavg, dzdyavg, detadzavg
    integer :: nx, ny, nz, nt, i, j, k, t
    
    nx = size(q,1)
    ny = size(q,2)
    nz = size(q,3)
    nt = size(q,4)
    dy = dx

    print *, "Calculating intermediate quantities globally..."

    ! Calculate the intermediate quantities globally
    allocate(dthetadeta(2:nx-1,2:ny-1,2:nz-1,nt), detadz(nx,ny,2:nz-1,nt),  &
         dvdeta(nx,2:ny,2:nz-1,nt), dudeta(2:nx,ny,2:nz-1,nt))
    do k = 2, nz-1
       detal = etau(k) - etau(k-1)
       detah = etau(k+1) - etau(k)
       c(1) = -detah / (detal * (detal+detah))
       c(2) = (detah**2 - detal**2) / (detal * detah * (detal+detah))
       c(3) = detal / (detah * (detal+detah))
       dthetadeta(2:nx-1,2:ny-1,k,:) = c(1)*theta(2:nx-1,2:ny-1,k-1,:)  &
            + c(2)*theta(2:nx-1,2:ny-1,k,:) + c(3)*theta(2:nx-1,2:ny-1,k+1,:)
       dvdeta(:,2:ny,k,:) = c(1)*v(:,2:ny,k-1,:) + c(2)*v(:,2:ny,k,:)  &
            + c(3)*v(:,2:ny,k+1,:)
       dudeta(2:nx,:,k,:) = c(1)*u(2:nx,:,k-1,:) + c(2)*u(2:nx,:,k,:)  &
            + c(3)*u(2:nx,:,k+1,:)      
       detadz(:,:,k,:) = (etaw(k+1)-etaw(k)) / (hght(:,:,k+1,:)-hght(:,:,k,:))
    end do

    allocate(dzdx(2:nx,ny,2:nz,nt))
    forall (i = 2:nx, k = 2:nz, t = 1:nt)
       dzdx(i,:,k,t) = .5 * (mf(i-1,:) + mf(i,:)) * (hght(i,:,k,t) - hght(i-1,:,k,t)) / dx
    end forall

    allocate(dzdy(nx,2:ny,2:nz,nt))
    forall (j = 2:ny, k = 2:nz, t = 1:nt)
       dzdy(:,j,k,t) = .5 * (mf(:,j-1) + mf(:,j)) * (hght(:,j,k,t) - hght(:,j-1,k,t)) / dy
    end forall
    
    allocate(vort(2:nx,2:ny,2:nz-1,nt))
    do t = 1, nt
       do k = 2, nz-1
          do j = 2, ny
             do i = 2, nx
                mfavg = .25 * (mf(i-1,j-1) + mf(i,j-1) + mf(i-1,j) + mf(i,j))
                detadzavg = .25 * (detadz(i-1,j-1,k,t) +  &
                     detadz(i,j-1,k,t) + detadz(i-1,j,k,t) + detadz(i,j,k,t))
                dzdxavg = .25 * (dzdx(i,j-1,k,t) + dzdx(i,j,k,t)  &
                     + dzdx(i,j-1,k+1,t) + dzdx(i,j,k+1,t))
                dzdyavg = .25 * (dzdy(i-1,j,k,t) + dzdy(i,j,k,t)  &
                     + dzdy(i-1,j,k+1,t) + dzdy(i,j,k+1,t))
                dvdeavg = .5 * (dvdeta(i-1,j,k,t) + dvdeta(i,j,k,t))
                dudeavg = .5 * (dudeta(i,j-1,k,t) + dudeta(i,j,k,t))
                vort(i,j,k,t) = mfavg * ((v(i,j,k,t) - v(i-1,j,k,t))/dx  &
                     + detadzavg * (dzdyavg*dudeavg - dzdxavg*dvdeavg)  &
                     - (u(i,j,k,t) - u(i,j-1,k,t))/dy)
             end do
          end do
       end do
    end do

    print *, "Calculating PV..."

    ! We've got all the pieces... now combine them!
    forall (i = 2:nx-1, j = 2:ny-1, k = 2:nz-1, t = 1:nt)
       q(i,j,k,t) = (detadz(i,j,k,t) * dthetadeta(i,j,k,t)  &
            * (f(i,j) + .25 * (vort(i,j,k,t) + vort(i+1,j,k,t)  &
            + vort(i,j+1,k,t) + vort(i+1,j+1,k,t)))) / rho(i,j,k,t)
    end forall
    
    deallocate(dthetadeta, dzdx, dzdy, dvdeta, dudeta, vort, detadz)

    ! Set boundary values to interior neighbors
    print *, "Filling boundaries..."
    q(2:nx-1,(/1,ny/),2:nz-1,:) = q(2:nx-1,(/2,ny-1/),2:nz-1,:)
    q((/1,nx/),:,2:nz-1,:) = q((/2,nx-1/),:,2:nz-1,:)
    q(:,:,(/1,nz/),:) = q(:,:,(/2,nz-1/),:)
  end subroutine pvor

  ! ===========================================================================
  ! Private functions (called only from this module) go below this line.
  ! ===========================================================================

  ! density_thp uses the ideal gas law to compute the density of dry air given
  ! its potential temperature (theta; K) and pressure (p; hPa)
  real function density_thp(theta, p)
    real, intent(in) :: theta, p
    
    real :: t
    
    t = temp_from_theta_p(theta, p, si=.false.)
    density_thp = 100 * p / (RGas * t)  ! * 100 since p is in hPa
  end function density_thp

  ! ===========================================================================

  ! theta_es computes the saturated equivalent potential temperature given
  ! potential temperature (theta; K) and pressure (p; hPa)
  real function theta_es(theta, p)
    real, intent(in) :: theta, p

    real :: tmpk, es, ws

    tmpk = temp_from_theta_p(theta, p, si=.false.)
    es = sat_vap_pres(tmpk)
    ws = mixr(es, p)

    theta_es = theta * exp(LatHeatVap * ws / (SpecHeatPresDry * tmpk))
  end function theta_es

  ! ===========================================================================

  ! theta_es_to_theta_fast computes the potential temperature given the
  ! saturated equivalent potential temperature (thes; K) and the
  ! pressure (p; hPa), using a lookup table.
  function theta_es_to_theta_fast(thes, p) result(theta)
    real, intent(in) :: thes, p
    real             :: theta
    
    integer, dimension(2), parameter :: PLim = (/ 50, 1020 /),  &
                                        TLim = (/ 2500, 3500 /)
    real, dimension(PLim(1):PLim(2),TLim(1):TLim(2)), save :: lookup
    real    :: t
    integer :: i, j
    logical :: firstTime = .true.
    
    ! Establish lookup table.
    if (firstTime) then
       firstTime = .false.
       do j = TLim(1), TLim(2)
          t = .1 * j
          do i = PLim(1), PLim(2)
             lookup(i,j) = theta_es_to_theta(t, real(i))
          end do
       end do
    end if
    
    ! Are we within lookup table?
    t = 10 * thes
    if (t > TLim(1) .and. t < TLim(2) .and. p > PLim(1) .and. p < PLim(2)) then
       i = int(p)
       j = int(t)
       theta = bilin_int(lookup(i,j), lookup(i,j+1), lookup(i+1,j),  &
            lookup(i+1,j+1), p-i, t-j)
    else
       theta = theta_es_to_theta(thes, p)
    end if
  end function theta_es_to_theta_fast
  
  ! theta_es_to_theta computes the potential temperature given the saturated
  ! equivalent potential temperature (thes; K) and the pressure (p; hPa),
  ! using the secant method.
  function theta_es_to_theta(thes, p) result(theta)
    real, intent(in) :: thes, p
    real             :: theta

    real,    parameter :: Delta = 1e-4, Epsil = 1e-4
    integer, parameter :: MaxIter = 200

    real    :: a, b, fa, fb, s
    integer :: k

    ! See if we're so dry that theta_es approximately equals theta
    a = thes
    if (abs(theta_es(a, p) - a) < .05) then
       theta = a
       return
    end if

    ! Make initial guess
    a = temp_from_theta_p(thes, p, si=.false.)
    b = a - 5
    fa = f(a)
    fb = f(b)
    
    ! Here's the secant method, at the end of which a is the temperature.
    do k = 2, MaxIter
       if (abs(fa) > abs(fb)) then
          call swap(a, b)
          call swap(fa, fb)
       end if
       s = (b - a) / (fb - fa)
       b = a
       fb = fa
       a = a - fa * s
       fa = f(a)
       if (abs(fa) < Epsil .or. abs(b-a) < Delta) exit
    end do
    
    theta = a * (1000. / p)**Kappa
    
  contains
    
    ! This function has a zero when the temperature satisfies the equation
    ! for saturated equivalent potential temperature.
    real function f(t)
      real, intent(in) :: t
      
      real :: es
      real :: ws
      
      es = sat_vap_pres(t)
      ws = mixr(es, p)
      
      ! Compare f to function theta_es
      f = thes - t * (1000./p)**kappa *  &
           exp(LatHeatVap * ws / (SpecHeatPresDry * t))
    end function f
  end function theta_es_to_theta
  
  ! ===========================================================================

  ! mixr calculates the mixing ratio (kg/kg) given the vapor pressure (e) and
  ! pressure (p). e and p must have the same units.
  real function mixr(e, p)
    real, intent(in) :: e, p

    mixr = Eps * e / (p - e)
  end function mixr

  ! ===========================================================================

  ! vap_pres calculates the vapor pressure given the mixing ratio (w; kg/kg)
  ! and pressure (p).  The vapor pressure will be in the same units as p.
  real function vap_pres(w, p)
    real, intent(in) :: w, p
    
    vap_pres = w * p / (Eps + w)
  end function vap_pres

  ! ===========================================================================

  ! sat_vap_pres uses the data from Curry & Webster (1999, p.113) to compute
  ! the saturation vapor pressure (hPa) given temperature (tmpk; K).
  real function sat_vap_pres(tmpk)
    real, intent(in) :: tmpk

    real, dimension(7), parameter :: a = (/ 6.1117675,  &
         .443986062, .143053301e-1, .265027242e-3, .302246994e-5,  &
         .203886313e-7, .638780966e-10 /)
    
    real    :: psum, tdiff
    integer :: i
    
    tdiff = tmpk - 273.15
    
    if (tdiff < -50) then
       sat_vap_pres = 6.11 * exp(LatHeatVap / RVap * (1. / 273.16 - 1. / tmpk))
    else
       psum = 0
       do i = 7, 2, -1
          psum = psum + a(i) * tdiff**(i-1)
       end do
       sat_vap_pres = psum + a(1)
    end if
  end function sat_vap_pres

  ! ===========================================================================

  ! temp_lcl calculates the temperature at the LCL given a parcel's temperature
  ! (tmpk; K) and relative humidity (rh; 0<rh<1) using Curry & Webster (1999,
  ! Eq. 6.30).
  real function temp_lcl(tmpk, rh)
    real, intent(in) :: tmpk, rh

    real :: a, b
    
    a = 1. / (tmpk - 55)
    b = log(rh) / 2840.
    temp_lcl = 1. / (a - b) + 55
  end function temp_lcl

  ! ===========================================================================

  ! rh calculates the relative humidity (0<rh<1) given the mixing ratio
  ! (w; kg/kg), pressure (p; hPa), and temperature (tmpk; K).
  real function rh(w, p, tmpk)
    real, intent(in) :: w, p, tmpk

    real :: es, e

    es = sat_vap_pres(tmpk)
    e = vap_pres(w, p)
    
    rh = e/es
  end function rh

  ! ===========================================================================

  ! swap interchanges the values of a and b.
  subroutine swap(a, b)
    real, intent(inout) :: a, b
    
    real :: temp
    
    temp = a
    a = b
    b = temp
  end subroutine swap

  ! ===========================================================================

  ! Given profiles of variables a (aProf) and b (bProf) and a target value of
  ! variable a (a), find_bounds outputs two adjacent values in aProf that
  ! bound a (aBound) and the associated values of b (bBound).
  subroutine find_bounds(aProf, bProf, a, aBound, bBound)
    real, dimension(:), intent(in)  :: aProf, bProf
    real,               intent(in)  :: a
    real, dimension(2), intent(out) :: aBound, bBound

    integer :: i, aSize
    logical :: found
    
    ! Make sure input arrays match
    aSize = size(aProf)
    if (aSize /= size(bProf)) then
       print *, "aProf and bProf have different sizes!", aSize, size(bProf)
       stop
    end if

    ! Find bounds
    found = .false.
    do i = 1, aSize-1
       if (aProf(i) >= a .and. aProf(i+1) < a .or.  &
            aProf(i) <= a .and. aProf(i+1) > a) then
          found = .true.
          aBound(1) = aProf(i)
          aBound(2) = aProf(i+1)
          bBound(1) = bProf(i)
          bBound(2) = bProf(i+1)
          exit
       end if
    end do

    if (.not. found) then
       print *, "Could not find bounds!"
       print *, aProf, a
       stop
    end if
  end subroutine find_bounds

  ! ===========================================================================

  ! lin_int performs linear interpolation given two data points (x1,y1) and 
  ! (x2,y2) and a known value for x.
  real function lin_int(x1, x2, y1, y2, x)
    real, intent(in) :: x1, x2, y1, y2, x
    
    if (x1 /= x2) then
       lin_int = y1 + (y2-y1)*(x-x1)/(x2-x1)
    else
       lin_int = y1
    end if
  end function lin_int

  ! ===========================================================================

  ! The first four values are the known quantities in the four corners.
  ! a and b are the fraction of the way we are between the left and right
  ! and lower and upper sides, respectively.
  real function bilin_int(ll, lr, ul, ur, a, b)
    real, intent(in) :: ll, lr, ul, ur, a, b
    
    real, dimension(4) :: w
    
    w(1) = (1-a) * (1-b)
    w(2) = a * (1-b)
    w(3) = (1-a) * b
    w(4) = a * b
    
    bilin_int = w(1) * ll + w(2) * lr + w(3) * ul + w(4) * ur
  end function bilin_int

  ! ===========================================================================

  ! calc_cape returns the convective available potential energy (J/kg) given
  ! an environmental temperature profile of tProf (K), profiles of pressure 
  ! (p; hPa) and height (z; m), and a parcel that, at its LCL, has a potential
  ! temperature of theta (K), mixing ratio of w (kg/kg), (saturated) equivalent
  ! potential temperature of thes (K), pressure of pLCL (hPa), and height of
  ! zLCL (m).
  real function calc_cape(tProf, p, z, theta, w, thes, pLCL, zLCL)
    real, dimension(:), intent(in) :: tProf, p, z
    real,               intent(in) :: theta, w, thes, pLCL, zLCL

    real, dimension(size(tProf)) :: tTraj, x, fx
    integer :: nz, k, kk, lev

    ! Initialize
    nz = size(z)
    calc_cape = 0
    
    ! Compute ascent trajectory and apply virtual temperature correction    
    tTraj = calc_traj(p, tProf(1), theta, w, thes, pLCL)
    
    ! Find LFC height (from bottom up)
    do k = 2, nz-1
       if (p(k) >= pLCL .or. tTraj(k) <= tProf(k)) cycle
       
       ! We just passed through LFC, but we need to allow for case where the
       ! LCL is in a region of positive area already due to a superadiabatic
       ! environmental lapse rate.
       if (tTraj(k-1) - tProf(k-1) >= 0) then
          x(k-1) = zLCL
       else
          x(k-1) = lin_int(tTraj(k-1)-tProf(k-1), tTraj(k)-tProf(k),  &
               z(k-1), z(k), 0.)
       end if
          
       ! Now find EL (from top down)
       ! Allow for the case where there is still positive area at the top of
       ! the sounding.
       if (tTraj(nz) >= tProf(nz)) then
          x(nz) = z(nz)
          kk = nz-1
       else
          do kk = nz-1, k, -1
             if (tTraj(kk) <= tProf(kk)) cycle
             exit
          end do
          ! We just passed through EL
          x(kk+1) = lin_int(tTraj(kk)-tProf(kk), tTraj(kk+1)-tProf(kk+1),  &
               z(kk), z(kk+1), 0.)
       end if
          
       ! Compute CAPE using composite trapezoid rule
       x(k:kk) = z(k:kk)
       where (tTraj(k-1:kk+1) > tProf(k-1:kk+1))
          fx(k-1:kk+1) = (tTraj(k-1:kk+1) - tProf(k-1:kk+1)) / tProf(k-1:kk+1)
       elsewhere
          fx(k-1:kk+1) = 0
       end where
       calc_cape = .5 * Grav *  &
            sum( (/ ( (x(lev)-x(lev-1)) * (fx(lev-1)+fx(lev)),  &
                      lev = k, kk+1 ) /) )
       exit
    end do
  end function calc_cape

  ! ===========================================================================

  ! calc_traj returns a temperature profile (traj; K) for a parcel undergoing
  ! pseudoadiabatic ascent.  The parcel is defined by an initial temperature
  ! (start; K), potential temperature (theta; K), mixing ratio (w; kg/kg),
  ! equivalent potential temperature (thes; K), and LCL pressure (pLCL; hPa),
  ! while the profile's pressure structure is provided as p (hPa).
  function calc_traj(p, start, theta, w, thes, pLCL) result(traj)
    real, dimension(:), intent(in) :: p
    real,               intent(in) :: start, theta, w, thes, pLCL
    real, dimension(size(p))       :: traj
    
    integer :: nz, k
    
    ! Initialize
    nz = size(p)
    
    ! Compute ascent trajectory and apply virtual temperature correction
    traj(1) = start
    do k = 2, nz
       if (p(k) >= pLCL) then
          traj(k) = temp_from_theta_p(theta, p(k), si=.false.) *  &
               ( 1 + EpsRM1 * w )
       else
          traj(k) = temp_from_theta_p(theta_es_to_theta_fast(thes, p(k)),  &
               p(k), si=.false.)
          traj(k) = traj(k) *  &
               ( 1 + EpsRM1 * mixr( sat_vap_pres(traj(k)), p(k) ) )
       end if
    end do
  end function calc_traj

  ! ===========================================================================
  
  ! calc_cin returns the convective inhibition (J/kg) given an environmental
  ! temperature profile of tProf (K), profiles of pressure (p; hPa) and height
  ! (z; m), and a parcel that, at its LCL, has a potential temperature of theta
  ! (K), mixing ratio of w (kg/kg), (saturated) equivalent potential
  ! temperature of thes (K), and pressure of pLCL (hPa).  The result is a
  ! nonnegative number.
  real function calc_cin(tProf, p, z, theta, w, thes, pLCL)
    real, dimension(:), intent(in) :: tProf, p, z
    real,               intent(in) :: theta, w, thes, pLCL
    
    real, dimension(size(tProf)) :: tTraj, x, fx
    integer :: nz, k, lev

    ! Initialize
    nz = size(z)
    calc_cin = 0
    
    ! Compute ascent trajectory and apply virtual temperature correction    
    tTraj = calc_traj(p, tProf(1), theta, w, thes, pLCL)
    
    ! Find LFC height (from bottom up)
    do k = 2, nz-1
       if (p(k) >= pLCL .or. tTraj(k) <= tProf(k)) cycle
       
       ! Calculate CIN
       ! LCL must be in a region of negative area for there to be CIN.
       if (tTraj(k-1) - tProf(k-1) < 0) then
          x(1:k-1) = z(1:k-1)
          x(k) = lin_int(tTraj(k-1)-tProf(k-1), tTraj(k)-tProf(k),  &
               z(k-1), z(k), 0.)
          where (tTraj(1:k) > tProf(1:k))
             fx(1:k) = 0
          elsewhere
             fx(1:k) = (tProf(1:k) - tTraj(1:k)) / tProf(1:k)
          end where
          calc_cin = .5 * Grav *  &
               sum( (/ ( (x(lev)-x(lev-1)) * (fx(lev-1)+fx(lev)),  &
                         lev = 2, k ) /) )
       end if
       exit
    end do
  end function calc_cin  

  ! ===========================================================================

  ! mean_wind computes a mean wind component in the layer given by heights
  ! (m AGL).  Provided as input are a wind component (wind; m/s), geopotential
  ! heights (z; m), the wind component at 10 m (wind10m; m/s), and the terrain
  ! height (ter; m).
  function mean_wind(wind, z, wind10m, ter, heights)
    real, dimension(:,:,:,:), intent(in) :: wind, z
    real, dimension(:,:,:),   intent(in) :: wind10m
    real, dimension(:,:),     intent(in) :: ter
    real, dimension(:),       intent(in) :: heights
    real, dimension(size(wind,1),size(wind,2),size(wind,4)) :: mean_wind

    real, dimension(size(z,1),size(z,2),size(z,3),size(z,4))     :: hagl
    real, dimension(size(z,1),size(z,2),size(heights),size(z,4)) :: windInterp
    real, dimension(2) :: haglBounds, windBounds
    integer :: nx, ny, nz, nt, numLevs, i, j, k, t

    nx = size(wind,1); ny = size(wind,2); nz = size(wind,3); nt = size(wind,4)
    numLevs = size(heights)

    ! Validate
    if (any(shape(wind) /= shape(z)) .or. nx /= size(wind10m,1) .or.  &
         ny /= size(wind10m,2) .or. nt /= size(wind10m,3) .or.  &
         nx /= size(ter,1) .or. ny /= size(ter,2))  &
         stop "Array mismatch in mean_wind6!"

    hagl = z -  &
        spread(source=spread(source=ter, dim=3, ncopies=nz), dim=4, ncopies=nt)
    
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             do k = 1, numLevs
                if (Heights(k) < hagl(i,j,k,t)) then
                   windInterp(i,j,k,t) = wind10m(i,j,t)
                else
                   call find_bounds(hagl(i,j,:,t), wind(i,j,:,t), Heights(k), &
                        haglBounds, windBounds)
                   windInterp(i,j,k,t) = lin_int(haglBounds(1), haglBounds(2),&
                        windBounds(1), windBounds(2), Heights(k))
                end if
             end do
          end do
       end do
    end do
    mean_wind = sum(windInterp, dim=3) / numLevs
  end function mean_wind

  ! ===========================================================================

  ! calc_ddz takes the vertical derivative of the input wind field given a set
  ! of vertical coordinate values (heights).
  function calc_ddz(wind, heights)
    real, dimension(:,:,:,:), intent(in) :: wind
    real, dimension(:),       intent(in) :: heights
    real, dimension(size(wind,1),size(wind,2),size(wind,3),size(wind,4)) ::  &
         calc_ddz

    integer :: nx, ny, nz, nt, i, j, t
    
    nx = size(wind,1); ny = size(wind,2); nt = size(wind,4)
    nz = size(heights)

    ! Validate
    if (nz /= size(wind,3)) stop "Array mismatch in calc_ddz!"
    
    do t = 1, nt
       do j = 1, ny
          do i = 1, nx
             calc_ddz(i,j,(/1,nz/),t) = (wind(i,j,(/2,nz/),t) -  &
                  wind(i,j,(/1,nz-1/),t)) / (heights((/2,nz/)) -  &
                  heights((/1,nz-1/)))
             calc_ddz(i,j,2:nz-1,t) = (wind(i,j,3:nz,t)-wind(i,j,1:nz-2,t)) / &
                  (heights(3:nz)-heights(1:nz-2))
          end do
       end do
    end do
  end function calc_ddz
end module diagnostics
