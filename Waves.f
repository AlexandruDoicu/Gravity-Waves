      program waves
      implicit none
      integer mmass, hwm_mod
      real    Tinf_scl
      real    aap(7), ap
      character*72 fileIRI
c
      integer year, month, day, hour, minute, second
      integer iyd
      real    LatRef, LonRef, altR
      real    f10p7, fbar
      real    sec, hrl
      real    d(9), temp(2)
      real    app(2), whm07(2)
      real*8  DipDeg, Dip, SinDip, CosDip
c
      integer pth, pthe, ptn, pto, ptn2, ptno, pto2, ptar
      real*8  re, pie, gzero, bolt, amu, Avogadro,
     &        Rgas, po180
      real*8  Prandtl, CpAverage, ThermCond
      real*8  M_H, M_He, M_N, M_O, M_N2, M_O2, M_Ar
      real*8  THETAv_N2, THETAv_O2
      real*8  gmid, ntot, mu_loc
      real*8  yH, yHe, yN, yO, yN2, yO2, yAr, Mbar
      real*8  Cv_mono, Cv_diat_base, Cv_vib_N2, Cv_vib_O2
      real*8  Cv_vib_molar, Cv_molar, Rspec
      real*8  v_from_R, v_from_H, H_m, rel_diff
      real*8  dz, alt
      logical ConstantGravAccel, TestGasConstant, TestSoundSpeed
c
      integer  IdzT, i, NModes, ishift
      real*8   Omega, Omega0, OmegaIm, Omega2
      real*8   kT
      real*8   TempMax, uT0Max
      real*8   z, zmax, dzT
      logical  write_data, more, fail
c
      integer  NOmegaIm, iOmegaIm, iSol
      real*8   AbsOmegaImSF, AbsOmegaIm, RMSOmegaIm
      real*8   AbsOmegaImMax, AbsOmegaImMin, DOmegaIm, DAbsOmegaIm
      real*8   OmegaImMin, OmegaImMax, dOm
      real*8   UMax, timeUMax, zUMax
      real*8   WMax, timeWMax, zWMax
      real*8   TMax, timeTMax, zTMax
c
      integer  nzWT, iT
      complex*16, allocatable:: URhat(:), UThat(:)
      complex*16, allocatable:: TWhat(:)
      complex*16, allocatable:: pWhat(:)
      complex*16, allocatable:: roWhat(:)
c
      complex*16, allocatable:: UThat0(:)
      complex*16, allocatable:: URhat0(:)
      complex*16, allocatable:: pWhat0(:)
      complex*16, allocatable:: TWhat0(:)
      complex*16, allocatable:: roWhat0(:)
c
      complex*16, allocatable:: UThatPlus(:)
      complex*16, allocatable:: URhatPlus(:)
      complex*16, allocatable:: pWhatPlus(:)
      complex*16, allocatable:: TWhatPlus(:)
      complex*16, allocatable:: roWhatPlus(:)
c
      complex*16, allocatable:: UThatMinus(:)
      complex*16, allocatable:: URhatMinus(:)
      complex*16, allocatable:: pWhatMinus(:)
      complex*16, allocatable:: TWhatMinus(:)
      complex*16, allocatable:: roWhatMinus(:)
c
      real*8, allocatable:: dennW0(:,:)
      real*8, allocatable:: tnW0(:)
      real*8, allocatable:: roW0(:)
      real*8, allocatable:: dAirW0(:)
      real*8, allocatable:: pW0(:)
      real*8, allocatable:: uT0(:), uP0(:)
      real*8, allocatable:: HScaleW0(:)
      real*8, allocatable:: HDensW0(:)
      real*8, allocatable:: MUW0(:)
      real*8, allocatable:: RgasW0(:)
      real*8, allocatable:: gW0(:)
      real*8, allocatable:: CvW0(:)
      real*8, allocatable:: CpW0(:)
      real*8, allocatable:: GammaW0(:)
      real*8, allocatable:: CsoundW0(:)
      real*8, allocatable:: alphaT(:)
c
      real*8, allocatable:: deniW0(:), deniW0S(:)
      real*8, allocatable:: denOW0(:), denOW0S(:)
      real*8, allocatable:: nu0niW(:), nu0inW(:)
      real*8, allocatable:: DiffCoef0W(:)
      real*8, allocatable:: DiffVelocityW0(:), DiffVelocityW0S(:)
      real*8, allocatable:: DDiffVelocityW0(:)
      real*8, allocatable:: DdeniW0(:), DdeniW0S(:)
c
      integer, allocatable:: izT(:)
      real*8, allocatable:: zT(:)
c
      integer    nzW
      integer    TypeSolMet, TypeLinModel
      integer    TypeBC, qBC, qNEV
      integer    TypeSclLayEq
      logical    DoIonDrag
      logical    DoNormEigVctZGEEV
      logical    DimLessStVct
      logical    DoScalingFT, DoScalingTime
      logical    DoSpatialFourierTransform
      real*8     zmin, HwT, zminDiffVelocity
      real*8     LambdaK, LambdaT, epsT, epsU
      real*8     Lambdazmin, LambdaZMax
      real*8     t0BC, u0BC, w0BC
      real*8     nu0niMax, FctIonDrag
      complex*16, allocatable:: KRadCT(:,:), EValT(:,:), EVctT(:,:,:)
      complex*16, allocatable:: AmplitudeT(:,:)
c
      complex*16, allocatable:: eT(:,:)
      complex*16, allocatable:: eT0(:,:)
      complex*16, allocatable:: eTPlus(:,:)
      complex*16, allocatable:: eTMinus(:,:)
c
      real*8, allocatable:: DroW0(:), DpW0(:), DtnW0(:),
     &                      DuT0(:), DHDensW0(:), DHDensW0S(:), DgW0(:)
c
      integer    NFFT, NPeriod
      real*8     RatioSigmaOmega, SigmaOmega, DOmega
      real*8     TimeMin, TimeShift
      real*8     MaxEVal1, MinEval4
      real*8     LregDen, LregDen1
c
      real*8     ErrU12, ErrW12, ErrT12
      real*8     ErrU13, ErrW13, ErrT13
      real*8     ErrU23, ErrW23, ErrT23
      logical    DoRMSErrors, ContinueWP
c
      real*8, allocatable:: OmegaK( : ), TimeK( : )
      complex*16, allocatable:: UTFFT(:,:), URFFT(:,:)
      complex*16, allocatable:: TFFT(:,:)
c
      complex*16, allocatable:: UTTime(:,:), URTime(:,:)
      complex*16, allocatable:: TTime(:,:)
c
      integer    NSol
      real*8, allocatable:: UMaxSol(:), WMaxSol(:)
      real*8, allocatable:: TMaxSol(:), OmegaImSol(:)
      complex*16, allocatable:: UTTimeSol(:,:,:), URTimeSol(:,:,:)
      complex*16, allocatable:: TTimeSol(:,:,:)
c
      integer  NZShift, NOShift, j, iO
      real*8   dZShift, dmin, dmax, dminG, eps
      integer, allocatable:: iZMap( : ), iOMap(:)
      real*8, allocatable:: dist(:,:)
      logical  DoGlobalCausalityCondition, crossG
c
      integer*4   count0, count1
      integer*4   count_rate, count_max
      real        time0, time1
c
      namelist /Parameters/
     &          LambdaK, zmin, zmax, nzW,

     &          TypeSolMet, TypeLinModel,
     &          DoIonDrag, zminDiffVelocity,
     &          Prandtl, qNEV, epsT, epsU, w0BC,

     &          NFFT, RatioSigmaOmega, NPeriod, TimeMin,

     &          NZShift, dZShift, NOmegaIm, AbsOmegaImSF,
     &          AbsOmegaImMax, AbsOmegaImMin, DOmegaIm,

     &          fileIRI, mmass, Tinf_scl, ap, hwm_mod
c
      call SYSTEM_CLOCK(count0, count_rate, count_max)
      time0 = count0 * 1.0 / count_rate
c
      write_data = .false.
c
      pie   = dacos(-1.d0)
      po180 = pie / 180.d0
c
      call SYSTEM_CLOCK(count0, count_rate, count_max)
      time0 = count0 * 1.0 / count_rate
c     -------------------------------------------------------------------------
c                             Read namelist
c     -------------------------------------------------------------------------
      open ( unit=10, file='WaveData.namelist'  )
      read(10,Parameters)
      close(unit=10)
c     -------------------------------------------------------------------------
c                    Read IRI-data file (day = day of the month)
c     -------------------------------------------------------------------------
      allocate( deniW0(nzW), deniW0S(nzW) )
      call READ_IRI_IONFILE
     i   ( fileIRI,
     i     nzW, zmin, zmax,
     o     year, month, day, hour, minute, second,
     o     LatRef, LonRef,
     o     DipDeg, f10p7, fbar,
     O     deniW0 )
c
c     --- sine and cosine of magnetic dipole angle ---
      Dip = DipDeg * po180
      SinDip = sin(Dip)
      CosDip = cos(Dip)
c     -------------------------------------------------------------------------
c                           Neutral-winds setup
c     -------------------------------------------------------------------------
      call IRI2HWM07_SETUP
     i   ( year, month, day, hour, minute, second, LonRef,
     o     iyd, sec, hrl )
c     -------------------------------------------------------------------------
c                        Control flags and parameters
c     -------------------------------------------------------------------------
      DimLessStVct = .false. !if false, the derivatives in the state vector have
                             !the dimension 1/m
      DoNormEigVctZGEEV = .false. !if false, the ZGEEV-eigenvector normalization
                                  !is not used
      TypeSclLayEq = 1 !if 1, the scaling of the layer equation is performed with
                       !respect to the positive and negtive values of the real
                       !part of the eigenvalues
      TypeBC = 2     ! 1: lower BC is imposed on modes;
                     ! 2: lower BC is imposed on state vector components
      qBC = qNEV !in the case TypeBC = 2, the boundary condition is imposed on
                 !the qBC component of the state vector, where
                 !qBC = 1 for u, qBC = 2 for w, and qBC = 3 for T
      NModes = 3 !number of wave modes:gravity, thermal and viscosity waves
c
c     --- regularization parameters for smoothing the density scale height ---
      LregDen  = 1.d0 !in km
      LregDen1 = 0.1d0
c
c     --- global causality condition at sample frequencies ---
      DoGlobalCausalityCondition  = .false. !if true, the global causality
                                            !condition is used, but the computation
                                            !time increases
      NOShift = 1   !number of sample frequencies at which the
                    !the global causality condition is checked
c
c     --- RMS errors for single-frequency waves ---
      DoRMSErrors = .false.  !it can be true for single-frequency waves,
                             !but after all solution methods GMMA, GMMN, and SMMA
                             !have been applied
c
c     --- min and max assumed vertical wavelengths for estimating wave period ---
      Lambdazmin = 125.d0
      LambdaZMax = 250.d0
c
c     --- spatial Fourier transform over the altitude ---
      DoSpatialFourierTransform = .false. !if true, a directory named
                                          !'FourierTransformSF' must be created
c     -------------------------------------------------------------------------
c                              Basic Parameters
c     -------------------------------------------------------------------------
      re    = 6371.d0       !Earth radius in km
      gzero = 980.665d0     !in cm/s**2
      bolt  = 1.38044d-16   !in cm**2 * g/(s**2 * K)
      amu   = 1.66054d-24   !atomic mass unit in g
      Avogadro = 6.02214086d+23 !in 1/mol
      Rgas  = 8.31446261815324d0 !gas constant in J/(mol*K)
      ConstantGravAccel = .true.
c
c     --- molar mass in kg/mol ---
      M_H  = 1.00794d-3
      M_He = 4.002602d-3
      M_N  = 14.0067d-3
      M_O  = 15.999d-3
      M_N2 = 28.0134d-3
      M_O2 = 31.9988d-3
      M_Ar = 39.948D-3
c
C     --- vibrational characteristic temperatures [K] for diatomics
      THETAv_N2 = 3395.d0
      THETAv_O2 = 2273.d0
c
c     --- pointers-neutrals ---
      pth   = 1    ! h
      pthe  = 5    ! he
      ptn   = 7    ! n
      pto   = 2    ! o
      ptn2  = 6    ! n2
      ptno  = 3    ! no
      pto2  = 4    ! o2
      ptar  = 8
c
      app(1) = ap
      app(2) = ap
c
      do i = 1,7
        aap(i) = ap
      enddo
c     ---------------------------------------------------------------------------
c                             Wave grid parameters
c     ---------------------------------------------------------------------------
      HwT = zmax - zmin     !altitude range
      dz = HwT / (nzW - 1)  !in km
      IdzT = 2              !2, 4
c
      zmax = zmin + HwT      !in km
c
      dzT  = IdzT * dz       !in km
      nzWT = nzW / IdzT      !number of isothermal regions
c     --------------------------------------------------------------------------
c                         Allocate wave grid quantities
c     ---------------------------------------------------------------------
      allocate( dennW0(nzW,8) )
      allocate( tnW0(nzW) )
      allocate( roW0(nzW) )
      allocate( dAirW0(nzW) )
      allocate( pW0(nzW) )
      allocate( uT0(nzW), uP0(nzW) )
      allocate( HScaleW0(nzW) ) !scale height in km
      allocate( HDensW0(nzW) )  !density scale height in m
      allocate( RgasW0(nzW) )   !specific gas constant in J/(kg*K)
      allocate( gW0(nzW) )      !gravitational acceleration in m/s**2
      allocate( MUW0(nzW) ) !mean molecular mass (dim.less)
      allocate( CvW0(nzW) ) ! c_v [J/(kg*K)]
      allocate( CpW0(nzW) ) ! c_p [J/(kg*K)]
      allocate( GammaW0(nzW) )   ! gamma = Cp/Cv [-]
      allocate( CsoundW0(nzW) )  ! speed of sound [m/s]
      allocate( alphaT(nzWT) )
c
      allocate( izT(nzWT+1), zT(nzWT+1) )
c
      allocate( denOW0(nzW), denOW0S(nzW) )
      allocate( nu0niW(nzW), nu0inW(nzW) )
      allocate( DiffCoef0W(nzW) )
      allocate( DiffVelocityW0(nzW), DiffVelocityW0S(nzW) )
      allocate( DDiffVelocityW0(nzW) )
      allocate( DdeniW0(nzW), DdeniW0S(nzW) )
c
      allocate( KRadCT(2*NModes, nzWT), EValT(2*NModes, nzWT) )
      allocate( EVctT(2*NModes,2*NModes, nzWT) )
      allocate( AmplitudeT(2*NModes,nzWT) )

      allocate( eT(2*NModes,nzWT) )
      allocate( eT0(2*NModes,nzWT) )
      allocate( eTPlus(2*NModes,nzWT) )
      allocate( eTMinus(2*NModes,nzWT) )
c
      allocate( UThat(nzWT) )
      allocate( URhat(nzWT) )
      allocate( pWhat(nzWT) )
      allocate( TWhat(nzWT) )
      allocate( roWhat(nzWT) )
c
      allocate( UThat0(nzWT) )
      allocate( URhat0(nzWT) )
      allocate( pWhat0(nzWT) )
      allocate( TWhat0(nzWT) )
      allocate( roWhat0(nzWT) )
c
      allocate( UThatPlus(nzWT) )
      allocate( URhatPlus(nzWT) )
      allocate( pWhatPlus(nzWT) )
      allocate( TWhatPlus(nzWT) )
      allocate( roWhatPlus(nzWT) )
c
      allocate( UThatMinus(nzWT) )
      allocate( URhatMinus(nzWT) )
      allocate( pWhatMinus(nzWT) )
      allocate( TWhatMinus(nzWT) )
      allocate( roWhatMinus(nzWT) )

      allocate( DroW0(nzW) )
      allocate( DpW0(nzW) )
      allocate( DtnW0(nzW) )
      allocate( DuT0(nzW) )
      allocate( DgW0(nzw) )
      allocate( DHDensW0(nzW) )
      allocate( DHDensW0S(nzW) )
c
      write(6,*)
      write(6,'(a)') 'IRI data:'
      write(6,'(a,i4,1x,i2,1x,i2)') '- Date (Y M D) = ',
     &      year, month, day
      write(6,'(a,i2,1x,i2,1x,i2)') '- Time (H M S) = ',
     &      hour, minute, second
      write(6,'(a,3x,i5)')'- IYD (YYDDD, day of year)   = ', iyd
      write(6,'(a,f12.3)') '- SEC (UT seconds)           = ', sec
      write(6,'(a,f12.6)') '- HRL (solar local time hrs) = ', hrl
      write(6,'(a,f10.3)') '- LatRef = ', LatRef
      write(6,'(a,f10.3)') '- LonRef = ', LonRef
      write(6,'(a,f10.3)') '- DipDeg = ', DipDeg
      write(6,'(a,f10.3)') '- F10P7  = ', f10p7
      write(6,'(a,f10.3)') '- FBAR   = ', fbar
c
      write(6,*)
      write(6,'(A)')    'Discretization grid:'
      write(6,'(2(a,f7.2))') '- Altitude bounds [km]: zmin = ', zmin,
     &      ',  zmax = ', zmax
      write(6,'(A,I5)') '- Number of altitude grid points = ', nzW
      write(6,'(A,I5)') '- Number of isothermal regions   = ', nzWT
      write(6,'(A,F7.2)')
     & '- Altitude discretization step [km]       = ', dz
      write(6,'(A,F7.2)')
     & '- Isothermal-region step in altitude [km] = ', dzT
      write(6,*)
c     --------------------------------------------------------------------------
c                           Atmospheric parameters
c     --------------------------------------------------------------------------
      TempMax =-1.d+10
      uT0Max  =-1.d+10
      TestGasConstant  = .true.
      TestSoundSpeed = .false.
      CpAverage = 0.d0
      do i = 1, nzW
          alt = zmin + (i-1)*dz !in km
C         ----------------------------------------------------------------------
C                     Neutral Atmosphere Empirical Model GTD7 in
C                              nrlsise00_modified.f
C         Input
C            iyd - YEAR AND DAY AS YYDFctD (day of year from 1 to 365 (or 366))
C                  (Year ignored in current model)
C            sec - UT(SEC)
C            alts(i,j) - ALTITUDE(KM)
C            LatRef(i,j) - GEODETIC LATITUDE(DEG)
C            LonRef(i,j) - GEODETIC LONGITUDE(DEG)
C            hrl - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
C            fbar - 81 day AVERAGE OF F10.7 FLUX (centered on day DFctD)
C            f10p7 - DAILY F10.7 FLUX FOR PREVIOUS DAY
C            aap - MAGNETIC INDEX(DAILY)
C            mmass - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
C                    CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
C                    MASS 17 IS Anomalous O ONLY.)
C
C        Output
C            d(1) - he number density(cm-3)
C            d(2) - o number density(cm-3)
C            d(3) - n2 number density(cm-3)
C            d(4) - o2 number density(cm-3)
C            d(5) - ar number density(cm-3)
C            d(6) - total mass density(g/cm3)
C            d(7) - h number density(cm-3)
C            d(8) - n number density(cm-3)
C            d(9) - anomalous O (see msis)
C            temp(1) - exospheric temperature
C            temp(2) - temperature at alt alts(i,j)
C         ----------------------------------------------------------------------
          altR = sngl(alt)
          call gtd7 ( iyd, sec, altR, LatRef, LonRef,
     i                hrl, fbar, f10p7, aap, mmass, Tinf_scl,
     o                d, temp )
c         ----------------------------------------------------------------------
c                          Number densities in 1/cm**3
c         ----------------------------------------------------------------------
          dennW0(i,pth ) = dble(d(7))        ! H
          dennW0(i,pthe) = dble(d(1))        ! He
          dennW0(i,ptn ) = dble(d(8))        ! N
          dennW0(i,pto ) = dble(d(2))        ! O
          dennW0(i,ptn2) = dble(d(3)) + 1.d-30! N2
          dennW0(i,pto2) = dble(d(4)) + 1.d-30! O2
          dennW0(i,ptar) = dble(d(5))        ! Ar

          tnW0(i) = dble(temp(2))             ! temperature [K]
c
C         --- NO parameterization (remains separate; NOT included in totals)
          dennW0(i,ptno) = 0.4d0 * exp(-3700.d0/tnW0(i))
     &                  * dennW0(i,pto2)
     &                  + 5.0d-7 * dennW0(i,pto)
c         -----------------------------------------------------------------------
c                       Atomic oxigen number density in 1/m**3
c         ----------------------------------------------------------------------
          denOW0(i) = 1.d+6 * dennW0(i,pto )
c         -----------------------------------------------------------------------
c                  Total mass density in g/cm**3 (= 10**3 kg/m**3 )
c         ----------------------------------------------------------------------
          roW0(i) = dble(d(6))
c         ----------------------------------------------------------------------
c                            Total number density in 1/cm**3
c         ----------------------------------------------------------------------
          ntot = 0.d0
          ntot = ntot + dennW0(i,pth )  ! H
          ntot = ntot + dennW0(i,pthe)  ! He
          ntot = ntot + dennW0(i,ptn )  ! N
          ntot = ntot + dennW0(i,pto )  ! O
          ntot = ntot + dennW0(i,ptn2)  ! N2
          ntot = ntot + dennW0(i,pto2)  ! O2
          ntot = ntot + dennW0(i,ptar)  ! Ar
          ntot = max(ntot, 1.d-300)     ! safety

          dAirW0(i) = ntot              ! total number density [cm^-3]
c         ----------------------------------------------------------------------
C                 Mean molecular mass mu (dimensionless; g/mol numerically)
C                      mu_loc  = Avogadro * rho / ntot [g/mol]
C                      MUW0(i) = mu_loc * 1e-3         [kg/mol]
c                where Avogadro = 6.02214086d+23 1/mol
c         ----------------------------------------------------------------------
          mu_loc  = Avogadro * roW0(i) / ntot !in g/mol
          MUW0(i) = 1.d-3 * mu_loc !in kg/mol
c         ----------------------------------------------------------------------
C           Mass-specific heat capacity at constant volume, CvW0(i) [J/(kg*K)]
c         ----------------------------------------------------------------------
C
C         --- molar fractions from available species (no Ar, no O* here) ---
          yH  = dennW0(i,pth ) / ntot
          yHe = dennW0(i,pthe) / ntot
          yN  = dennW0(i,ptn ) / ntot
          yO  = dennW0(i,pto ) / ntot
          yN2 = dennW0(i,ptn2) / ntot
          yO2 = dennW0(i,pto2) / ntot
          yAr = dennW0(i,ptar) / ntot
c
C         --- mixture mean molar mass [kg/mol] ---
          Mbar = yH*M_H + yHe*M_He + yN*M_N + yO*M_O
     &         + yN2*M_N2 + yO2*M_O2 + yAr*M_Ar
c
C         --- base contributions: monatomic = 3/2 Rgas, diatomic = 5/2 Rgas ---
c         --- where Rgas  = 8.31446261815324d0 J/(mol*K)   ---
          Cv_mono      = 1.5d0 * Rgas
          Cv_diat_base = 2.5d0 * Rgas
c
C         --- vibrational corrections for N2 and O2
          Cv_vib_N2 = Cv_vib_molar(tnW0(i), THETAv_N2, Rgas)
          Cv_vib_O2 = Cv_vib_molar(tnW0(i), THETAv_O2, Rgas)
c
C         --- Cv_molar: molar Cv of the mixture [J/(mol*K)] ---
          Cv_molar = (yH + yHe + yN + yO + yAr) * Cv_mono
     &              + yN2 * (Cv_diat_base + Cv_vib_N2)
     &              + yO2 * (Cv_diat_base + Cv_vib_O2)
c
C         --- convert to mass-specific: Cv_mass = Cv_molar / Mbar [J/(kg*K)] ---
          CvW0(i) = Cv_molar / Mbar  !J/(kg*K)
c         ----------------------------------------------------------------------
C                            Specific gas constant
c         ----------------------------------------------------------------------
          RgasW0(i) = Rgas / Mbar  !J/(kg*K)
          Rspec = Rgas / MUW0(i) ! J/(kg*K)
c         ----------------------------------------------------------------------
C                 Check specific gas constant: Rspec = RgasW0(i)
c         ----------------------------------------------------------------------
          if ( TestGasConstant ) then
            rel_diff = ABS(Rspec - RgasW0(i)) / MAX(RgasW0(i), 1.d-12)
            if (rel_diff .gt. 1.d-2) then
              WRITE(6,'(A,I6,3(1X,1PE11.4))')
     &       'Gas-constant error:', i, RgasW0(i), Rspec, rel_diff
            end if
c
            rel_diff = ABS(Mbar - MUW0(i)) / MAX(MUW0(i), 1.d-12)
            if (rel_diff .gt. 1.d-2) then
              WRITE(6,'(A,I6,3(1X,1PE11.4))')
     &       'Molar-mass error:', i, MUW0(i), Mbar, rel_diff
            end if
          end if
c         ----------------------------------------------------------------------
c                      Convert total mass density in kg/m**3
c         ----------------------------------------------------------------------
          roW0(i) = 1.d+3 * roW0(i)  !in kq/m**3
c         ----------------------------------------------------------------------
C                         Presure in N/m**2 = kg / (m*s**2)
c                              p = rho * Rspec * T
c         ----------------------------------------------------------------------
          pW0(i) = roW0(i) * RgasW0(i) * tnW0(i)  !N/m**2
c         ----------------------------------------------------------------------
C         Gravitational acceleration in m/s**2 from gzero = 980.665d0 cm/s**2
c         ----------------------------------------------------------------------
          gmid   = gzero * ( re / (re + alt) )**2 !in cm/s**2
          gW0(i) = 1.d-2 * gmid !in m/s**2
          if ( ConstantGravAccel ) then
            gW0(i) = 1.d-2 * gzero !in m/s**2
          end if
c        ------------------------------------------------------------------------
c                            Derivative (1/g)dg/dz in 1/m
c        ------------------------------------------------------------------------
          DgW0(i) = -2.0d-3 / (re + alt)   ! in 1/m
          if ( ConstantGravAccel )  DgW0(i) = 0.d0
c         ----------------------------------------------------------------------
C                              Scale height in km
C                       H = bolt * T / (amu * mu_loc * g)
c         or
c                       H = Rspec * T / g(z),
c         since Rspec = bolt / (amu * mu_loc) = Rgas / mu_loc
c         ----------------------------------------------------------------------
c         HScaleW0(i) = 1.e-5 * bolt * tnW0(i)
c     &               / ( mu_loc * amu * gmid )
c
          HScaleW0(i) = 1.d-3 * Rspec * tnW0(i) / gW0(i) !in km
c         ----------------------------------------------------------------------
c                           gamma = 1 + Rspec / Cv, where
c                     Rspec = Rgas / Mbar, with Mbar = MUW0(i)
c         ----------------------------------------------------------------------
          GammaW0(i) = 1.d0 + Rspec / CvW0(i)
c         ----------------------------------------------------------------------
C         Mass-specific heat capacity at constant pressure, CpW0(i) [J/(kg*K)]
c         ----------------------------------------------------------------------
          CpW0(i)   = GammaW0(i) * CvW0(i)
          CpAverage = CpAverage + CpW0(i)
c         ----------------------------------------------------------------------
C                              Sound speed in m/s
c         ----------------------------------------------------------------------
          v_from_R = SQRT( GammaW0(i) * Rspec * tnW0(i) )
          CsoundW0(i) = v_from_R
c         ----------------------------------------------------------------------
c                              Check sound speed
c         ----------------------------------------------------------------------
          if (TestSoundSpeed ) then
            H_m  = 1.d+3 * HScaleW0(i) ! in m
            v_from_H = SQRT( GammaW0(i) * gW0(i) * H_m )
            rel_diff = ABS(v_from_R - v_from_H) / MAX(v_from_R, 1.d-12)
            if (rel_diff .gt. 1.d-2) then
              WRITE(6,'(A,I6,3(1X,1PE11.4))')
     &       'Sound-speed error:', i, v_from_R, v_from_H, rel_diff
            end if
          end if
c         ----------------------------------------------------------------------
c                                Maximum temperature
c         ----------------------------------------------------------------------
          TempMax = max( TempMax, tnW0(i) )
c         ----------------------------------------------------------------------
C             Horizontal wind model HWM93 covering all altitude regions
c         INPUT:
C           iyd - YEAR AND DAY AS YYDDD
C           sec - UT(SEC)  (Not important in lower atmosphere)
C           ALT - ALTITUDE(KM)
C           GLAT - GEODETIC LATITUDE(DEG)
C           GLONG - GEODETIC LONGITUDE(DEG)
C           hrl - LOCAL APPARENT SOLAR TIME(HRS)
C           fbar - 3 MONTH AVERAGE OF F10.7 FLUX (Use 150 in lower atmos.)
C           f10p7 - DAILY F10.7 FLUX FOR PREVIOUS DAY ( " )
C           app - Two element array with
C                app(1) = MAGNETIC INDEX(DAILY) (use 4 in lower atmos.)
C                app(2)=CURRENT 3HR ap INDEX (used only when SW(9)=-1.)
C           Note:  Ut, Local Time, and Longitude are used independently in the
C                  model and are not of equal importance for every situation.
C                  For the most physically realistic calculation these three
C                  variables should be consistent.
C         OUTPUT
C           whm07(1) = MERIDIONAL (m/sec + Northward) = uTheta
C           whm07(2) = ZONAL (m/sec + Eastward) = uPhi
C         ----------------------------------------------------------------------
          if (hwm_mod .eq. 7) then
            call hwm07 ( iyd, sec, altR, LatRef, LonRef,
     .                   hrl, fbar, f10p7, app, whm07 )
          else if (hwm_mod .eq. 14) then
            call hwm14 ( iyd,sec, altR, LatRef, LonRef,
     .                   hrl, fbar, f10p7, app, whm07 )
          else
            call gws5 ( iyd, sec, altR, LatRef, LonRef,
     .                  hrl, fbar, f10p7, app, whm07 )
          endif
c         ----------------------------------------------------------------------
c                                Velocities
c                ux0 = - uT0: x-direction toward geographic south;
c                uy0 =   uP0: y-direction toward east
c         ----------------------------------------------------------------------
          uT0(i) = - dble(whm07(1)) ! in m/sec
          uP0(i) =   dble(whm07(2)) ! in m/sec
c         ----------------------------------------------------------------------
c                        Maximum horizontal velocity
c         ----------------------------------------------------------------------
          uT0Max = max( uT0Max, abs(uT0(i)) )
      enddo
      CpAverage = CpAverage / dble(nzW)
      ThermCond = 3.34d-7 * CpAverage / Prandtl
c     -------------------------------------------------------------------------
c              Boundary values for temperature, and velocities
c     -------------------------------------------------------------------------
      t0BC = epsT * TempMax !in K
      u0BC = epsU * uT0Max  !in m/s
c     -------------------------------------------------------------------------
c                    Altitude indices of each isothermal region
c                         izT(nzWT)´ and zT(nzWT)
c     -------------------------------------------------------------------------
      izT(1) = 1
      zT(1)  = zmin
      iT = 1
      do i = 2, nzW
        alt = zmin + (i-1)*dz
        if ( mod(i,IdzT) .eq. 1 ) then
          iT = it + 1
          zT(it)  = alt
          izT(it) = i
        end if
      end do
      if ( write_data ) then
        do iT = 1, nzWT-1
          i = izT(iT) + IdzT / 2
          write(6,'(100(2x,1pe11.4))') zT(it), zT(it+1),
     &    zmin + (izT(it)-1)*dz, zmin + (izT(it+1)-1)*dz,
     &    zmin + (i-1)*dz,
     &    0.5 * ( zT(it) + zT(it+1) )
        end do
        read(5,*)
      end if
c     -------------------------------------------------------------------------
c                          Wavenumbers kT in 1/km
c     -------------------------------------------------------------------------
      kT = 2.d0 * pie / LambdaK  !in 1/km
c     -------------------------------------------------------------------------
c                        Wave period in minutes
c     -------------------------------------------------------------------------
      call WavePeriod
     i   ( nzW, LambdaK, uT0, gW0,
     i     GammaW0, CsoundW0, HScaleW0,
     i     Lambdazmin, LambdaZMax,
     o     LambdaT )
c     -------------------------------------------------------------------------
c                           Frequency in 1/s
c     ------------------------------------------------------------------------
      Omega  = 2.d0 * pie / (60.d0 * LambdaT)
      Omega2 = Omega * Omega
c     -------------------------------------------------------------------------
c                          Reference frequency Omega0 = Omega
c     -------------------------------------------------------------------------
      Omega0 = Omega
c     -------------------------------------------------------------------------
c                         Parameter alphaT(nzWT)
c     -------------------------------------------------------------------------
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        alphaT(iT) = 1.d0 / ( kT * HScaleW0(i) )
      end do
c     ---------------------------------------------------------------------------
c                           Derivatives in SI
c     ---------------------------------------------------------------------------
c
c     --- smooth deniW0 (density of O+)---
      call csplineD
     i  ( .false., nzW, deniW0, LregDen,
     o     deniW0S )
      do i = 1, nzW
        deniW0(i) = deniW0S(i)
      end do
c
c     --- smooth denOW0 (density of O ) ---
      call csplineD
     i  ( .false., nzW, denOW0, LregDen,
     o     denOW0S )
      do i = 1, nzW
        denOW0(i) = denOW0S(i)
      end do
c
c     --- ion-drag parameters ---
      nu0niMax = - 1.d+20
      do i = 1, nzW
c
c       --- neutral-ion frequency in 1/s ---
        nu0niW(i) = 7.22d-17 * tnW0(i)**0.37d0 * deniW0(i)
        nu0niMax = max( nu0niMax, nu0niW(i) )
c
c       --- ion-neutral frequency in 1/s ---
        nu0inW(i) = nu0niW(i) * denOW0(i) / deniW0(i)
C
c       --- diffusion coefficient in m**2/s
        DiffCoef0W(i) = 1.03d+3 * tnW0(i) / nu0inW(i)
      end do
      FctIonDrag = nu0niMax / Omega
c
c     --- compute derivatives per meter ---
      do i = 2, nzW - 1
        DtnW0(i)   = 1.d-3 * 0.5* ( tnW0(i+1)   - tnW0(i-1) ) / dz
        DuT0(i)    = 1.d-3 * 0.5* ( uT0(i+1)    - uT0(i-1) )  / dz
        DdeniW0(i) = 1.d-3 * 0.5* ( deniW0(i+1) - deniW0(i-1) ) / dz
      end do
      DtnW0(1)   = DtnW0(2)
      DtnW0(nzW) = DtnW0(nzW-1)
c
      DuT0(1)   = DuT0(2)
      DuT0(nzW) = DuT0(nzW-1)
c
      DdeniW0(1)   = DdeniW0(2)
      DdeniW0(nzW) = DdeniW0(nzW-1)
c
c     --- smooth DdeniW0 ---
      call csplineD
     i  ( .false., nzW, DdeniW0, LregDen1,
     o     DdeniW0S )
      do i = 1, nzW
        DdeniW0(i) = DdeniW0S(i)
      end do
c
c     --- pressure and mass density derivatives per meter ---
      DroW0 = 0.d0
      DpW0  = 0.d0
      do i = 1, nzW
        DroW0(i) = - roW0(i) * ( 1.e-3 / HScaleW0(i)
     &             + DtnW0(i) / tnW0(i) )
        DpW0(i)  = - gW0(i) * roW0(i)
      end do
c
c     --- density scale height in m ---
      if ( TypeLinModel .eq. 2 ) then
        do i = 1, nzW
          HDensW0(i)  = pW0(i) / ( roW0(i) * gW0(i) ) !in m
          DHDensW0(i) = 0.d0
        end do
      else
        do i = 1, nzW
          HDensW0(i) = - roW0(i) / DroW0(i) !in m
c
c         DHDensW0(i) = -gw0(i) * HDensW0(i)**2 * roW0(i)
c     &  * ( DgW0(i) - DtnW0(i) / tnW0(i) ) / pw0(i)
        end do
c
c       --- derivative of HDensW0 (dimensionless ) with smoothing ---
        do i = 2, nzW - 1
          DHDensW0(i) = 1.d-3 * 0.5*( HDensW0(i+1) - HDensW0(i-1) ) / dz
        end do
        DHDensW0(1)   = DHDensW0(2)
        DHDensW0(nzW) = DHDensW0(nzW-1)
        call csplineD
     i    ( .false., nzW, DHDensW0, LregDen,
     o       DHDensW0S )
        do i = 1, nzW
          DHDensW0(i) = DHDensW0S(i)
        end do
      end if
c
c     --- diffusion velocity along magnetic field ---
      do i = 1, nzW
        DiffVelocityW0(i) = DiffCoef0W(i)
     &         * ( DdeniW0(i) / deniW0(i) + DtnW0(i) / tnW0(i) )
     &         + gW0(i) / nu0inW(i)
        DiffVelocityW0(i) =  SinDip * DiffVelocityW0(i)
      end do
      call csplineD
     i  ( .false., nzW, DiffVelocityW0, LregDen,
     o     DiffVelocityW0S )
      do i = 1, nzW
        DiffVelocityW0(i) = DiffVelocityW0S(i)
      end do
c
c     --- derivative of diffusion velocity ---
      do i = 2, nzW - 1
        DDiffVelocityW0(i) = 1.d-3 * 0.5* ( DiffVelocityW0(i+1)
     &                     - DiffVelocityW0(i-1) ) / dz
      end do
      DDiffVelocityW0(1)   = DDiffVelocityW0(2)
      DDiffVelocityW0(nzW) = DDiffVelocityW0(nzW-1)
c
c     --- switching function for diffusion velocity ---
      if ( DoIonDrag ) then
        do i = 1, nzW
          alt = zmin + (i-1)*dz !in km
          if ( alt .le. zminDiffVelocity ) then
            DiffCoef0W(i) = 0.d0
            gW0(i) = 0.d0
            DiffVelocityW0(i)  = 0.d0
            DDiffVelocityW0(i) = 0.d0
          end if
        end do
      end if
c     -------------------------------------------------------------------------
c                        Write data to screen
c     -------------------------------------------------------------------------
      write(6,'(a)') 'Frequency and time period:'
      write(6,'(2(a,1pe11.4))')
     &   '- Time frequency [1/s]                = ', Omega
      write(6,'(2(a,1pe11.4))')
     &   '- Maximum ion-drag frequency [1/s]    = ', nu0niMax
       write(6,'(2(a,1pe11.4))')
     &   '- Ion-drag frequency / Time frequency = ', FctIonDrag
      write(6,'(2(a,1pe11.4))')
     &   '- Time period [min]                   = ', LambdaT
      write(6,*)
c
      write(6,'(a)') 'Lower boundary values of perturbed quantities:'
      write(6,'(2(a,1pe11.4))')
     &   '- Temperature [K]           = ', t0BC
      write(6,'(2(a,1pe11.4))')
     &   '- Horizontal velocity [m/s] = ', u0BC
      write(6,'(2(a,1pe11.4))')
     &   '- Vertical velocity [m/s]   = ', w0BC
      if ( qNEV .eq. 1 ) then
        write(6,'(a)') '- Horizontal-velocity criterion is chosen'
      else if ( qNEV .eq. 2 ) then
        write(6,'(a)') '- Vertical-velocity criterion is chosen'
      else
        write(6,'(a)') '- Temperature criterion is chosen'
      end if
      write(6,*)
c     ---------------------------------------------------------------------------
c                   Write atmospheric parameters to files
c     ---------------------------------------------------------------------------
      call  PrintAtmosphericParameters
     i    ( nzW, nzWT, IdzT, izT, zmin, zminDiffVelocity, dz,
     i      tnW0, roW0, pW0, uT0, HScaleW0,
     i      HDensW0, CvW0, GammaW0, CsoundW0,
     i      deniW0S, nu0niW, nu0inW, DiffVelocityW0,
     i      DroW0, DpW0, DtnW0, DuT0, DHDensW0,
     i      DdeniW0, DDiffVelocityW0 )
c     -------------------------------------------------------------------------
c
c
c                            Single-frequency wave
c
c
c     -------------------------------------------------------------------------
      OmegaIm = - AbsOmegaImSF
c     .......................................................................
c                           Eigenvalues and Eigenvectors
c     .......................................................................
      call EigenValuesEigenVectors
     i   ( NModes, nzW, nzWT, IdzT, izT,
     i     Omega, OmegaIm, Omega0, kT, alphaT, Prandtl,
     i     TypeLinModel, DoIonDrag,
     i     DimLessStVct, DoNormEigVctZGEEV, qNEV,
C
     i     tnW0, roW0, pW0, uT0,
     i     CvW0, gW0, HDensW0, GammaW0, CsoundW0,
     i     DroW0, DtnW0, DuT0, DHDensW0,
     i     SinDip, CosDip, nu0niW, deniW0,
     i     nu0inW, DiffCoef0W, DiffVelocityW0,
     i     DdeniW0, DDiffVelocityW0,
C
     o     KRadCT, EValT, EVctT )
c
      call PrintWAvenumbersAndEigenvaluesSF
     i   ( TypeLinModel, NModes, nzWT, IdzT, izT, zmin, dz,
     i     KRadCT, EvalT )
C     .......................................................................
C                                 Method 1
c                    Global-Matrix Method for Amplitudes
C     .......................................................................
      if ( TypeSolMet .eq. 1 ) then
        call WaveAmplitudesModelGMMA
     i     ( TypeSclLayEq, TypeBC, qBC, NModes, nzWT, kT, zT,
     i       EValT, EVctT,
     o       AmplitudeT )
c
        call WaveParametersModelGMMAComplete
     i     ( NModes, nzWT, nzW, IdzT, izT,
     i       Omega, OmegaIm, Omega0, kT,
     i       EValT, EVctT, AmplitudeT,
     i       TypeBC, DimLessStVct, .true.,
     i       t0BC, u0BC, w0BC, qBC,
     i       tnW0, roW0, pW0, uT0, HDensW0,
     o       UThat, UThat0, UThatPlus, UThatMinus,
     o       URhat, URhat0, URhatPlus, URhatMinus,
     o       TWhat, TWhat0, TWhatPlus, TWhatMinus,
     o       pWhat, roWhat )
c
        call PrintWaveParametersModelGMMAComplete
     i     ( TypeLinModel, nzWT, izT, zmin, dz, OmegaIm,
     i       UThat, UThat0, UThatPlus, UThatMinus,
     i       URhat, URhat0, URhatPlus, URhatMinus,
     i       TWhat, TWhat0, TWhatPlus, TWhatMinus,
     i       pWhat, roWhat )
C     .......................................................................
C                              ·Method 2
c         Global-Matrix Method for Nodal (grid-point) Values
C     .......................................................................
      else if ( TypeSolMet .eq. 2 ) then
        call WaveAmplitudesModelGMMN
     i     ( TypeSclLayEq, TypeBC, qBC, NModes, nzWT, kT, zT,
     i       EValT, EVctT,
     o       eT, eT0, eTPlus, eTMinus )
c
        call WaveParametersModelGMMN
     i     ( NModes, nzWT, nzW, IdzT, izT,
     i       Omega,  OmegaIm, Omega0, kT,
     i       eT, DimLessStVct, .true.,
     i       t0BC, u0BC, w0BC, qBC,
     i       tnW0, roW0, pW0, uT0, HDensW0,
     o       UThat, URhat, TWhat, pWhat, roWhat )
c
        call PrintWaveParametersModelGMMNandSMMA
     i     ( TypeLinModel, TypeSolMet, nzWT, izT, zmin, dz, OmegaIm,
     i       UThat, URhat, TWhat, pWhat, roWhat )
C     .......................................................................
C                               Method 3
c                  Scattering-Matrix Method for Amplitudes
C     .......................................................................
      else  if ( TypeSolMet .eq. 3 ) then
        call WaveAmplitudesModelSMMA
     i     ( TypeSclLayEq, TypeBC, qBC, NModes, nzWT, kT, zT,
     i       EValT, EVctT,
     o       AmplitudeT )

        call WaveParametersModelSMMA
     i     ( NModes, nzWT, nzW, IdzT, izT,
     i       Omega, OmegaIm, Omega0, kT,
     i       EVctT, AmplitudeT,
     i       DimLessStVct, .true.,
     i       t0BC, u0BC, w0BC, qBC,
     i       tnW0, roW0, pW0,  uT0, HDensW0,
     o       UThat, URhat, TWhat, pWhat, roWhat )
c
        call PrintWaveParametersModelGMMNandSMMA
     i     ( TypeLinModel, TypeSolMet, nzWT, izT, zmin, dz, OmegaIm,
     i       UThat, URhat, TWhat, pWhat, roWhat )
      end if
C     .......................................................................
C               Spatial Fourier transform over the altitude
c     If true, a directory named 'FourierTransformSF' must be created.
C     .......................................................................
      if ( DoSpatialFourierTransform ) then
        call SpatialFourierTransform
     i     ( TypeSolMet, TypeLinModel,
     i       nzWT, IdzT, izT, dz,
     i       UThat, URhat, TWhat )
      end if
C     .......................................................................
C                        RMS errors in solution methods
C     .......................................................................
      if ( DoRMSErrors ) then
        call RMSErrors
     I     ( TypeLinModel, nzWT,
     o       ErrU12, ErrW12, ErrT12,
     o       ErrU13, ErrW13, ErrT13,
     o       ErrU23, ErrW23, ErrT23)
        write(6,*)
        write(6,'(a)') 'RMS errors in solution methods'
        write(6,'(3(a,1pe11.4))')
     &  '- GMMA-GMMN: error in U = ', ErrU12,
     &  '; error in W = ', ErrW12,'; error in T = ', ErrT12
c
        write(6,'(3(a,1pe11.4))')
     &  '- GMMA-SMMA: error in U = ', ErrU13,
     &  '; error in W = ', ErrW13,'; error in T = ', ErrT13
C
        write(6,'(3(a,1pe11.4))')
     &  '- GMMN-SMMA: error in U = ', ErrU23,
     &  '; error in W = ', ErrW23,'; error in T = ', ErrT23
      end if
c     -------------------------------------------------------------------------
c
c
c                       Time-dependent wave packet
c
c
c     -------------------------------------------------------------------------
      write(6,*)
      write(6,'(a)')
     &'Enter TRUE to continue with the time-dependent wave packet '//
     &'calculation, or FALSE to stop.'
      read(5,*) ContinueWP
      if ( .not. ContinueWP ) goto 1
c
      allocate( UTFFT(nzWT, NFFT) )
      allocate( URFFT(nzWT, NFFT) )
      allocate( TFFT(nzWT, NFFT) )
c
      allocate( UTTime(nzWT, NFFT) )
      allocate( URTime(nzWT, NFFT) )
      allocate( TTime(nzWT, NFFT) )
c
      allocate( OmegaK( NFFT ), TimeK( NFFT ) )
c
      allocate( iZMap( NZShift ) )
      allocate( iOMap( NOShift ) )
      allocate( dist(NFFT,nzWT) )
c
      allocate( UMaxSol(NOmegaIm), WMaxSol(NOmegaIm) )
      allocate( TMaxSol(NOmegaIm), OmegaImSol(NOmegaIm) )
      allocate( UTTimeSol(nzWT, NFFT, NOmegaIm) )
      allocate( URTimeSol(nzWT, NFFT, NOmegaIm) )
      allocate( TTimeSol(nzWT, NFFT, NOmegaIm) )
c
c     --- Set type of scaliing ---
      DoScalingFT   = .false.
      DoScalingTime = .not. DoScalingFT
c     -----------------------------------------------------------------------
c                        Fourier-transform parameters
c     -----------------------------------------------------------------------
      call TimeFTParameters
     i   ( NFFT, Omega0, RatioSigmaOmega, NPeriod, TimeMin,
     o     TimeShift, OmegaK, TimeK, SigmaOmega, DOmega )
c     -----------------------------------------------------------------------
c                     Set of altitudes for frequency shift
c     -----------------------------------------------------------------------
      call AltitudeMap
     i   ( NZShift, dZShift, nzWT, izT, zmin, dz,
     o     iZMap )
c     -----------------------------------------------------------------------
c                   Set of frequencies for frequency shift
c     -----------------------------------------------------------------------
      call OmegaMap
     i   ( NFFT, NOShift,
     o     iOMap)
c     -----------------------------------------------------------------------
c                                  Step 1
c    Check the maximum imaginary frequency shift (e.g., AbsOmegaImMin = 1.e-6)
c          by estimating the error in the source function reconstruction
c     -----------------------------------------------------------------------
      write(6,*)
      write(6,'(a)') 'Step 1: Maximum imaginary frequency shift '//
     &'with RMS-error in source-function reconstruction:'
      OmegaIm = - AbsOmegaImMin
      call SourceFunctionReconstruction
     i   ( NFFT, TimeShift,
     i     Omega0, OmegaIm, RatioSigmaOmega, DOmega,
     i     OmegaK, TimeK, RMSOmegaIm)
      write(6,'(2(a,1pe11.4))')
     &'- RMS-error = ', RMSOmegaIm, ' for OmegaIm = ', OmegaIm
      if ( RMSOmegaIm .gt. 5.d-3 ) then
        write(6,'(a,1pe11.4)')
     & 'Source function reconstruction test not passed '//
     & 'for AbsOmegaImMin = ', AbsOmegaImMin
        stop
      else
        write(6,'(a,1pe11.4)')'- Maximum imaginary frequency shift = ',
     &        -AbsOmegaImMin
      end if
c     -----------------------------------------------------------------------
c                        If required, determine
c      the minimum imaginary frequency shift (e.g., AbsOmegaImMax) with
c              global causality condition at specified frequencies
c     -----------------------------------------------------------------------
      if ( DoGlobalCausalityCondition ) then
        write(6,*)
        write(6,'(a)') 'Minimum imaginary frequency shift '//
     &                 'with global causality condition:'
        AbsOmegaIm = -1.e+20
        fail = .false.
        do ishift = 1, NOShift
          iO = iOMap(ishift)
          Omega = OmegaK(iO)
          call GlobalCausalityCondition
     i       ( NModes, nzW, nzWT, IdzT, izT,
     i         Omega, AbsOmegaImMin, AbsOmegaImMax, DOmegaIm,
     i         Omega0, kT, alphaT, Prandtl,
     i         TypeLinModel, DoIonDrag,
     i         DimLessStVct, DoNormEigVctZGEEV, qNEV,
C
     i         tnW0, roW0, uT0,
     i         CvW0, gW0, HDensW0, GammaW0, CsoundW0,
     i         DroW0, DtnW0, DuT0, DHDensW0,
     i         SinDip, CosDip, nu0niW, deniW0,
     i         nu0inW, DiffCoef0W, DiffVelocityW0,
     i         DdeniW0, DDiffVelocityW0,
C
     o         fail, OmegaIm )
          if ( fail ) then
            write(6,'(a,1pe9.2)')
     &     'Global causality condition fails '//
     &     'at sample frequency = ', Omega
            write(6,'(a)')
     &     'Increase AbsOmegaImMax in the input namelist'
            stop
          else
            write(6,'(2(a,1pe11.4))')
     &      '- Imaginary frequency shift = ', OmegaIm,
     &      ' for sample frequency = ', Omega
            if ( abs(OmegaIm) .gt. AbsOmegaIm ) then
              AbsOmegaIm = abs(OmegaIm)
            end if
          end if
        end do
        AbsOmegaImMax = AbsOmegaIm
        write(6,'(a,1pe11.4)')'- Minimum imaginary frequency shift = ',
     &        -AbsOmegaImMax
      end if
c     -----------------------------------------------------------------------
c                                Step 2
c     Check the minimum imaginary frequency shift (e.g., AbsOmegaImMax = 1.e-4)
c            by estimating the error in the source function reconstruction
c     -----------------------------------------------------------------------
      write(6,*)
      write(6,'(a)') 'Step 2: Minimum imaginary frequency shift '//
     &'with RMS-error in source-function reconstruction:'
      more = .true.
      fail = .false.
      AbsOmegaIm = AbsOmegaImMax
      do while ( more )
        OmegaIm = - AbsOmegaIm
        call SourceFunctionReconstruction
     i     ( NFFT, TimeShift,
     i       Omega0, OmegaIm, RatioSigmaOmega, DOmega,
     i       OmegaK, TimeK, RMSOmegaIm)
c       write(6,'(2(a,1pe11.4))')
c       '- RMS-error = ', RMSOmegaIm, ' for OmegaIm = ', OmegaIm
        if ( RMSOmegaIm .lt. 5.d-3 ) then
          more = .false.
        else
          if ( AbsOmegaIm .ge. AbsOmegaImMin + DOmegaIm ) then
            AbsOmegaIm = AbsOmegaIm - DOmega
          else
            more = .false.
            fail = .true.
          end if
        end if
      end do
      if ( fail ) then
        write(6,'(a,1pe11.4)')
     & 'Source function reconstruction test not passed '//
     & 'for AbsOmegaImMax = ', AbsOmegaImMax
        stop
      else
        AbsOmegaImMax = AbsOmegaIm
        write(6,'(2(a,1pe11.4))')
     &  '- RMS-error = ', RMSOmegaIm, ' for OmegaIm = ', OmegaIm
        write(6,'(a,1pe11.4)')'- Minimum imaginary frequency shift = ',
     &        -AbsOmegaImMax
      end if
c     -----------------------------------------------------------------------
c                            Step 3
c                  Adjustment of internal OmegaIm-grid
c     -----------------------------------------------------------------------
      DAbsOmegaIm = AbsOmegaImMax - AbsOmegaImMin
      if ( DAbsOmegaIm .lt. 0.5d0*DOmegaIm ) then
        AbsOmegaImMin = AbsOmegaImMax - 0.5d0*DOmegaIm
        DOmegaIm = DOmegaIm / 8.d0

      else if ( (0.5d0*DOmegaIm.le. DAbsOmegaIm) .and.
     &          (DAbsOmegaIm .lt. 1.5d0*DOmegaIm) ) then
        DOmegaIm = DOmegaIm / 4.d0
      else if ( (1.5d0*DOmegaIm.le. DAbsOmegaIm) .and.
     &          (DAbsOmegaIm .lt. 2.5d0*DOmegaIm) ) then
        DOmegaIm = 2.d0 * DOmegaIm / 3.d0
      end if
c
      write(6,*)
      write(6,'(a)')
     &'Step 3: Adjustment of internal imaginary frequency interval:'
      write(6,'(3(a,1pe11.4))')
     &  '- OmegaImMin = ', -AbsOmegaImMax,
     &  ', OmegaMax = ',-AbsOmegaImMin,
     &  ', DOmegaIm = ', DOmegaIm
c     -----------------------------------------------------------------------
c                              Step 4
c             Maximum frequency shift (AbsOmegaImMin) with
c       layerwise causality condition at specified altitude levels
c     -----------------------------------------------------------------------
      write(6,*)
      write(6,'(a)') 'Step 4: Maximum imaginary frequency shift '//
     &               'with layerwise causality condition:'
      AbsOmegaIm = -1.e+20
      fail = .false.
      do ishift = 1, NZShift
        iT = iZMap(ishift)
        call LayerwiseCausalityCondition
     i     ( NFFT, iT, OmegaK,
     i       NModes, nzW, nzWT, IdzT, izT,
     i       AbsOmegaImMin, AbsOmegaImMax, DOmegaIm,
     i       Omega0, kT, alphaT, Prandtl,
     i       TypeLinModel, DoIonDrag,
     i       DimLessStVct, DoNormEigVctZGEEV, qNEV,
C
     i       tnW0, roW0, uT0,
     i       CvW0, gW0, HDensW0, GammaW0, CsoundW0,
     i       DroW0, DtnW0, DuT0, DHDensW0,
     i       SinDip, CosDip, nu0niW, deniW0,
     i       nu0inW, DiffCoef0W, DiffVelocityW0,
     i       DdeniW0, DDiffVelocityW0,
C
     o       fail, OmegaIm, dmin, dmax )
        z = zmin + (izT(iT) - 1)*dz
        if ( fail ) then
          write(6,'(a,1pe9.2)')
     &   'Layerwise causality condition fails '//
     &   'at altitude = ', z
          stop
        else
          if ( z .le. 200.d0 ) then
            write(6,'(2(a,1pe9.2),a)')
     &     '- Imaginary frequency shift = ', OmegaIm,
     &     ' at altitude = ', z,' yields:'
            write(6,'(4x,2(a,1pe9.2))')
     &     'MinDist(curves) = ', dmin,
     &     ',  MaxDist(curves) = ', dmax
          end if
          if ( abs(OmegaIm) .gt. AbsOmegaIm ) then
            AbsOmegaIm = abs(OmegaIm)
          end if
        end if
      end do
      AbsOmegaImMin = AbsOmegaIm
c
      OmegaImMin = -AbsOmegaImMax
      OmegaImMax = -AbsOmegaImMin
      write(6,'(2(a,1pe13.4))')
     &   'Imaginary frequency interval: OmegaImMin = ',
     &    OmegaImMin, ', OmegaImMax = ', OmegaImMax
c     -----------------------------------------------------------------------
c                               Step 5
c                        Loop over frequency shift
c     -----------------------------------------------------------------------
      dOm = 0.d0
      if ( NOmegaIm .gt. 1 ) then
        dOm = (OmegaImMax - OmegaImMin) / DBLE(NOmegaIm - 1)
      end if
      UTTimeSol = (0.d0,0.d0)
      URTimeSol = (0.d0,0.d0)
      TTimeSol  = (0.d0,0.d0)
      iSol = 0
      do iOmegaIm = 1, NOmegaIm
        if ( NOmegaIm .eq. 1 ) then
          OmegaIm = OmegaImMax
        else
          OmegaIm = OmegaImMin + dOm * DBLE(iOmegaIm - 1)
        end if
c       -----------------------------------------------------------------------
c                        Computation for OmegaIm
c       -----------------------------------------------------------------------
        write(6,*)
        write(6,'(a,1pe13.4,a)')
     & 'Step 5: Wave parameters for OmegaIm = ', OmegaIm,':'
        write(6,'(a)') '- Compute wave parameters in frequency domain'
        MaxEVal1 = -1.d+20
        MinEVal4 =  1.d+20
        do i = 1, NFFT
          Omega = OmegaK(i)
c
          call EigenValuesEigenVectors
     i       ( NModes, nzW, nzWT, IdzT, izT,
     i         Omega, OmegaIm, Omega0, kT, alphaT, Prandtl,
     i         TypeLinModel, DoIonDrag,
     i         DimLessStVct, DoNormEigVctZGEEV, qNEV,
C
     i         tnW0, roW0, pW0, uT0,
     i         CvW0, gW0, HDensW0, GammaW0, CsoundW0,
     i         DroW0, DtnW0, DuT0, DHDensW0,
     i         SinDip, CosDip, nu0niW, deniW0,
     i         nu0inW, DiffCoef0W, DiffVelocityW0,
     i         DdeniW0, DDiffVelocityW0,
C
     o         KRadCT, EValT, EVctT )
c
c         --- distance between eigenvalue curves ---
          do iT = 1, nzWT
            dist(i,iT) = dble(EValT(4,iT)) - dble(EValT(1,iT))
          end do
c
c         --- Global-Matrix Method for Amplitudes ---
          if ( TypeSolMet .eq. 1 ) then
            call WaveAmplitudesModelGMMA
     i         ( TypeSclLayEq, TypeBC, qBC, NModes, nzWT, kT, zT,
     i           EValT, EVctT,
     o           AmplitudeT )
c
            call WaveParametersModelGMMA
     i         ( NModes, nzWT, nzW, IdzT, izT,
     i           Omega, OmegaIm, Omega0, kT,
     i           EVctT, AmplitudeT,
     i           DimLessStVct, DoScalingFT,
     i           t0BC, u0BC, w0BC, qBC,
     i           tnW0, roW0, pW0, uT0, HDensW0,
     o           UThat, URhat, TWhat, pWhat, roWhat )
c
c         --- Global-Matrix Method for Nodal Values ---
          else if ( TypeSolMet .eq. 2 ) then
            call WaveAmplitudesModelGMMN
     i         ( TypeSclLayEq, TypeBC, qBC, NModes, nzWT, kT, zT,
     i           EValT, EVctT,
     o           eT, eT0, eTPlus, eTMinus )
c
            call WaveParametersModelGMMN
     i         ( NModes, nzWT, nzW, IdzT, izT,
     i           Omega,  OmegaIm, Omega0, kT,
     i           eT, DimLessStVct, DoScalingFT,
     i           t0BC, u0BC, w0BC, qBC,
     i           tnW0, roW0, pW0, uT0, HDensW0,
     o           UThat, URhat, TWhat, pWhat, roWhat )
c
c         --- Scattering-Matrix Method for Amplitudes ---
          else  if ( TypeSolMet .eq. 3 ) then
            call WaveAmplitudesModelSMMA
     i         ( TypeSclLayEq, TypeBC, qBC, NModes, nzWT, kT, zT,
     i           EValT, EVctT,
     o           AmplitudeT )

            call WaveParametersModelSMMA
     i         ( NModes, nzWT, nzW, IdzT, izT,
     i           Omega, OmegaIm, Omega0, kT,
     i           EVctT, AmplitudeT,
     i           DimLessStVct, DoScalingFT,
     i           t0BC, u0BC, w0BC, qBC,
     i           tnW0, roW0, pW0,  uT0, HDensW0,
     o           UThat, URhat, TWhat, pWhat, roWhat )
          end if
c
          do iT = 1, nzWT
            UTFFT(iT,i) = UThat(iT)
            URFFT(iT,i) = URhat(iT)
            TFFT(iT,i)  = TWhat(iT)
          end do
c
          do iT = 1, nzWT
            MaxEVal1 = max( MaxEVal1, dble(EValT(1,iT)) )
            MinEVal4 = min( MinEVal4, dble(EValT(4,iT)) )
          end do
        end do !end do NFFT (end do Omega)
c
c       --- check eigenvalue curve crossing for this OmegaIm ---
        crossG = .false.
        dminG  = 1.d+20
        eps    = 1.0d-10
        do iT = 1, nzWT
          dmin =  1.0d300
          dmax = -1.0d300
          do i = 1, NFFT
            dmin = min(dmin, dist(i,iT)) !minimum over Omega for altitude iT
            dmax = max(dmax, dist(i,iT))
          end do
          dminG = min(dminG,dmin)  !minimum over all altitudes
          if (dmin.lt.-eps .and. dmax.gt.eps) crossG = .true.
          if (abs(dmin).le.eps .or. abs(dmax).le.eps) crossG = .true.
          if ( crossG ) then
            write(6,'(a,1pe11.4)')
     &      '- Layerwise causality condition is '//
     &      'not satisfied (eigenvalue curves cross) at altitude:',
     &       zmin + (izT(iT)-1)*dz
             exit
          end if
        end do
        if ( .not. crossG ) then
          write(6,'(a,1pe11.4)')'- Layerwise causality condition is '//
     &   'satisfied (eigenvalue curves do not cross) with'
          write(6,'(a,1pe11.4)')'  a minimum '//
     &   'distance between curves over all altitudes of ', dminG
        end if
c
c       --- check global causality condition ---
        if ( MaxEVal1 .lt. MinEVal4 ) then
            write(6,'(a)')
     &      '- Global causality condition is satisfied:'
            write(6,'(4x,2(a,1pe11.4))')
     & 'Max(EVal(asc)) = ', MaxEVal1, ' < Min(EVal(desc)) = ', MinEVal4
        else
            write(6,'(a)')
     &      '- Global causality condition is not satisfied:'
            write(6,'(4x,2(a,1pe11.4))')
     & 'Max(EVal(asc)) = ', MaxEVal1, ' > Min(EVal(desc)) = ', MinEVal4
        end if
c       .....................................................................
c             Solutions that satisfy the LC condition at all altitudes
c       .....................................................................
        if ( .not. crossG ) then
          iSol = iSol + 1
c
          write(6,'(a)') '- Compute wave parameters in time domain'
          call WaveParametersInTime
     i       ( NFFT, TimeShift, nzWT,
     i         Omega0, OmegaIm, RatioSigmaOmega, DOmega,
     i         OmegaK, TimeK, UTFFT, URFFT, TFFT,
     i         DoScalingTime, t0BC, u0BC, w0BC, qBC,
     o         UTTime, URTime, TTime )
c
          write(6,'(a)') '- Maximum wave parameters in time domain:'
          call MaxWaveParametersInTime
     i       ( NFFT, nzWT, izT, zmin, dz,
     i         TimeK, UTTime, URTime, TTime,
     o         UMax, timeUMax, zUMax,
     o         WMax, timeWMax, zWMax,
     o         TMax, timeTMax, zTMax )
c
          write(6,'(4x,3(a,1pe11.4))')
     &   'Maximum horizontal velocity [m/s] = ', UMax,
     &   ' at time [hr] = ', timeUMax/60.d0, ' and level [km] = ', zUMax
c
          write(6,'(4x,3(a,1pe11.4))')
     &   'Maximum vertical velocity [m/s]   = ', WMax,
     &   ' at time [hr] = ', timeWMax/60.d0, ' and level [km] = ', zWMax
C
          write(6,'(4x,3(a,1pe11.4))')
     &   'Maximum temperature [K]           = ', TMax,
     &   ' at time [hr] = ', timeTMax/60.d0, ' and level [km] = ', zTMax
c
C         --- store solutions ---
          do j = 1, NFFT
            do iT = 2, nzWT
                UTTimeSol(iT,j,isol) = UTTime(iT,j)
                URTimeSol(iT,j,isol) = URTime(iT,j)
                TTimeSol(iT,j,isol)  = TTime(iT,j)
            end do
          end do
          UMaxSol(isol)    = UMax
          WMaxSol(isol)    = WMax
          TMaxSol(isol)    = TMax
          OmegaImSol(isol) = OmegaIm
        end if !end if admissible solutions (.not. crossG)
c
      end do !end do iOmegaIm
      NSol = isol
c
      write(6,*)
      if ( NSol .eq. 0 ) then
         write(6,'(A)')
     &  'No admissible solution satisfying the LC ' //
     &  'condition at all altitudes was found'
         stop
      else
c
c       --- select admissible solution ---
        UTTime = (0.d0,0.d0)
        URTime = (0.d0,0.d0)
        TTime  = (0.d0,0.d0)
        call SelectSolution
     i     ( nzWT, NFFT, NOmegaIm, NSol,
     i       UMaxSol, WMaxSol, TMaxSol, OmegaImSol,
     i       UTTimeSol, URTimeSol, TTimeSol,
     o       UMax, WMax, TMax, OmegaIm,
     o       UTTime, URTime, TTime )
        if ( NSol .eq. 1 ) then
           write(6,'(A)')
     &    'One admissible solution satisfying the LC ' //
     &    'condition at all altitudes was found;'
           write(6,'(A,1PE13.4)')
     &    'Frequency shift of this solution: OmegaIm = ',
     &     OmegaIm
        else
           write(6,'(A)')
     &    'Several admissible solutions satisfying the LC ' //
     &    'condition at all altitudes were found;'
           write(6,'(A,1PE13.4)')
     &    'Frequency shift of the selected solution: OmegaIm = ',
     &     OmegaIm
        end if
      end if
c
c     --- Gnuplot figures for velocities and temperature ---
      write(6,'(a)')
     &'- Maximum wave parameters in time domain '//
     &'for the final solution:'
      call PrintWaveParsModelWP
     i   ( NFFT, nzWT, izT, zmin, dz,
     i     TimeShift, Omega0, RatioSigmaOmega, TimeK,
     i     UTTime, URTime, TTime )
c
c     --- deallocate ---
      deallocate( UTFFT )
      deallocate( URFFT )
      deallocate( TFFT )
c
      deallocate( UTTime )
      deallocate( URTime )
      deallocate( TTime )
c
      deallocate( OmegaK, TimeK )
c
      deallocate( iZMap, iOMap, dist )
c
      deallocate( UMaxSol, WMaxSol )
      deallocate( TMaxSol, OmegaImSol )
      deallocate( UTTimeSol )
      deallocate( URTimeSol )
      deallocate( TTimeSol )
c
1     continue
c
      write(6,*)
      call SYSTEM_CLOCK(count1, count_rate, count_max)
      time1 = count1 * 1.0 / count_rate
      if ( TypeSolMet .eq. 1 ) then
        write(6,'(a,f7.2)')
     & 'Elapsed time for GMMA solution method [sec] = ',
     &  time1 - time0
      else if ( TypeSolMet .eq. 2 ) then
        write(6,'(a,f7.2)')
     & 'Elapsed time for GMMN solution method [sec] = ',
     &  time1 - time0
      else
        write(6,'(a,f7.2)')
     & 'Elapsed time for SMMA solution method [sec] = ',
     &  time1 - time0
      end if
c
c     --- deallocate ---
      deallocate( dennW0 )
      deallocate( tnW0 )
      deallocate( roW0 )
      deallocate( dAirW0 )
      deallocate( pW0 )
      deallocate( uT0, uP0 )
      deallocate( HScaleW0 )
      deallocate( HDensW0 )
      deallocate( RgasW0 )
      deallocate( gW0 )
      deallocate( MUW0 )
      deallocate( CvW0 )
      deallocate( CpW0 )
      deallocate( GammaW0 )
      deallocate( CsoundW0 )
      deallocate( alphaT )
!
      deallocate( izT, zT )
c
      deallocate( KRadCT, EValT )
      deallocate( EVctT )
      deallocate( AmplitudeT )

      deallocate( eT )
      deallocate( eT0 )
      deallocate( eTPlus )
      deallocate( eTMinus )
!
      deallocate( UThat )
      deallocate( URhat )
      deallocate( pWhat )
      deallocate( TWhat )
      deallocate( roWhat )
c
      deallocate( UThat0 )
      deallocate( URhat0 )
      deallocate( pWhat0 )
      deallocate( TWhat0 )
      deallocate( roWhat0 )
!
      deallocate( UThatPlus )
      deallocate( URhatPlus )
      deallocate( pWhatPlus )
      deallocate( TWhatPlus )
      deallocate( roWhatPlus )
!
      deallocate( UThatMinus )
      deallocate( URhatMinus )
      deallocate( pWhatMinus )
      deallocate( TWhatMinus )
      deallocate( roWhatMinus )

      deallocate( DroW0 )
      deallocate( DpW0 )
      deallocate( DtnW0 )
      deallocate( DuT0 )
      deallocate( DgW0 )
      deallocate( DHDensW0 )
      deallocate( DHDensW0S )

      deallocate( deniW0, deniW0S )
      deallocate( denOW0, denOW0S )
      deallocate( nu0niW, nu0inW )
      deallocate( DiffCoef0W )
      deallocate( DiffVelocityW0, DiffVelocityW0S )
      deallocate( DDiffVelocityW0 )
      deallocate( DdeniW0, DdeniW0S )

      end
c     *************************************************************************
      subroutine WavePeriod
     i         ( nzW, LambdaK, uT0, gW0,
     i           GammaW0, CsoundW0, HScaleW0,
     i           Lambdazmin, LambdaZMax,
     o           LambdaT )
      implicit none
      integer  nzW
      real*8   LambdaK, uT0(nzW), gW0(nzW)
      real*8   GammaW0(nzW), CsoundW0(nzW), HScaleW0(nzW)
      real*8   Lambdazmin, LambdaZMax
c
      real*8   LambdaT
c
c     --- local ---
      integer  i, ILambda1, ILambda2, ILambda3
      integer  ILambdaMin, ILambdaMax
      real*8   Pi, OmegaG, OmegaA, Csound
      real*8   kzmin, kZMax, kh, kX, kz, kz2, kh2, b, c, dd, alpha
      real*8   Omega, Omega1, Omega2, Omega3
      real*8   Omega1Min, Omega2Min, Omega3Min
      real*8   Lambda1, Lambda2, Lambda3
c
      Pi = acos(-1.d0)
      kh = 2.d0*Pi / LambdaK !in 1/km
      kX = 1.d-3 * kh !in 1/m
c
      kzmin  = 2.d0 * Pi / LambdaZMax !in 1/km
      kZMax  = 2.d0 * Pi / Lambdazmin !in 1/km
c
      Omega1Min = 1.d+20
      Omega2Min = 1.d+20
      Omega3Min = 1.d+20
      do i = 1, nzW
c
c       --- frequency of gravity waves ---
        OmegaG = sqrt(GammaW0(i) - 1.d0) * gW0(i) / CsoundW0(i)
c
c       --- frequency of acoustic waves ---
        OmegaA = GammaW0(i) * gW0(i) / (2.d0 * CsoundW0(i))
c
c       --- speed of sound ---
        Csound = 1.d-3 * CsoundW0(i) !in km/s
c       -------------------------------------------------------------------------
c                     Omega for kzmin (LambdaZMax = 250 km)
c       -------------------------------------------------------------------------
        kz  = kzmin
        kz2 = kz * kz
        kh2 = kh * kh
        b   = 0.5d0*( OmegaA**2 + Csound**2 * (kh2 + kz2) )
        c   = Csound**2 * kh2 * OmegaG**2
        dd  = b**2 - c
        if ( dd .lt. 0 ) then
          write(6,*) 'Time frequency cannot be computed'
          stop
        else
          Omega = b - sqrt(dd)
          if ( Omega .le. 0.d0 ) then
            write(6,*) 'Time frequency cannot be computed'
            stop
          else
            Omega = sqrt(Omega)
            if ( Omega .ge. OmegaG ) then
              write(6,*) 'The wave is not a gravity wave'
              stop
            end if
          end if
        end if
        Omega1 = Omega + kX * uT0(i)
        Omega1Min = min( Omega1Min, Omega1)
c       -------------------------------------------------------------------------
c                     Omega for kZMax (Lambdazmin = 50 km)
c       -------------------------------------------------------------------------
        kz  = kZMax
        kz2 = kz * kz
        kh2 = kh * kh
        b   = 0.5d0*( OmegaA**2 + Csound**2 * (kh2 + kz2) )
        c   = Csound**2 * kh2 * OmegaG**2
        dd  = b**2 - c
        if ( dd .lt. 0 ) then
          write(6,*) 'Time frequency cannot be computed'
          stop
        else
          Omega = b - sqrt(dd)
          if ( Omega .le. 0.d0 ) then
            write(6,*) 'Time frequency cannot be computed'
            stop
          else
            Omega = sqrt(Omega)
            if ( Omega .ge. OmegaG ) then
              write(6,*) 'The wave is not a gravity wave'
              stop
            end if
          end if
        end if
        Omega2 = Omega + kX * uT0(i)
        Omega2Min = min( Omega2Min, Omega2)
c       -------------------------------------------------------------------------
c                     Omega for kZ = 0 (LambdaZMax = infinity)
c       -------------------------------------------------------------------------
        alpha  = 1.d0 / ( kh * HScaleW0(i) )
        Omega  = OmegaG / sqrt( 1.d0 + alpha**2 / 4.d0 )
        Omega3 = Omega + kX * uT0(i)
        Omega3Min = min( Omega3Min, Omega3)
      end do
      Lambda1 = 2.d0*Pi / (60.d0*Omega1Min)
      Lambda2 = 2.d0*Pi / (60.d0*Omega2Min)
      Lambda3 = 2.d0*Pi / (60.d0*Omega3Min)
c
      ILambda1 = 10 * ceiling(Lambda1 / 10.0)
      ILambda2 = 10 * ceiling(Lambda2 / 10.0)
      ILambda3 = 10 * ceiling(Lambda3 / 10.0)
c
      write(6,'(a)') 'Time period of the wave:'
      write(6,'(a,2x,f7.2)')
     & '- Period in minutes for LambdaZ =   infinity: ', Lambda3
      write(6,'(2(a,f7.2))')'- Period in minutes for LambdaZ = ',
     &      LambdaZmax, ' km:   ', Lambda1
      write(6,'(2(a,f7.2))')'- Period in minutes for LambdaZ = ',
     &      Lambdazmin, ' km:   ', Lambda2
c
c     --- safety rule ---
      if ( LambdaK .lt. 390.d0 ) then
        ILambdaMin = 30
      else if ( (390.d0 .le. LambdaK) .and. (LambdaK .lt. 540.d0) ) then
        ILambdaMin =  40
      else if ( (540.d0 .le. LambdaK) .and. (LambdaK .lt. 640.d0) ) then
        ILambdaMin =  50
      else
        ILambdaMin =  60
      end if
c
      write(6,'(2(a,i4),a)')
     &'Enter the TIME PERIOD of the wave, preferably in the range ',
     & ILambdaMin, ' min  - ', ILambda2,' min'
      read(5,*) LambdaT
      write(6,*)
      end
c     *************************************************************************
      subroutine msistim ( iyr,iday,hrl,glong,iyd,secut )

!     msistim calculates time parameters for the
!     nrlmsise00 neutral atmosphere model.

!     the arguments are defined as follows:

!       iyr    the julian year
!       iday   the day of the year
!       hr     the local time in hours
!       glong  the geocentric longitude in degrees east
!       iyd    the year and day in the form yydd
!       secut  the universal time in seconds

      iyd    = 1000 * mod(iyr,100) + iday
      hrut   = hrl - glong /15.

      do while ( hrut .lt. 0.  )
        hrut = hrut + 24.
      enddo

      do while ( hrut. ge. 24. )
        hrut = hrut - 24.
      enddo

      secut  = hrut * 3600.

      return
      end
c     *************************************************************************
      subroutine AltitudeMap
     i         ( NZShift, dZShift, nzWT, izT, zmin, dz,
     o           iZMap )
      implicit none
      integer NZShift, nzWT
      integer izT(nzWT+1)
      real*8  dZShift, zmin, dz
c
      integer iZMap( NZShift )
c
      integer k, iT, i, iT0
      real*8  h, eps, z
c
      do k = 1, NZShift
        h = zmin + dble(k-1) * dZShift
        eps = 1.d+20
        do iT = 1, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz
          if ( abs(z-h) .lt. eps ) then
            eps = abs(z-h)
            iT0 = iT
          end if
        end do
        iZMap(k) = iT0
      end do
      end
c     ************************************************************************
      subroutine OmegaMap
     i         ( NFFT, NOShift,
     o           iOMap)
      implicit none
      integer NFFT, NOShift
c
      integer iOMap(NOShift)
c
      integer k
      real*8  step
c
      if (NOShift .eq. 1) then
        iOMap(1) = 1

      else if (NOShift .eq. 2) then
        iOMap(1) = 1
        iOMap(2) = NFFT

      else
c
C       --- NOShift > 2: equidistant indices including 1 and NFFT ---
        step = dble(NFFT - 1) / dble(NOShift - 1)
        do k = 1, NOShift
          iOMap(k) = 1 + int( dble(k-1)*step + 0.5d0 )
        end do
      end if
      end
C     *************************************************************************
C
C
C
c
C                        Eigenvalues and Eigenvectors
C
C
C
C
C     ****************************************·*********************************
      subroutine EigenValuesEigenVectors
     i         ( NModes, nzW, nzWT, IdzT, izT,
     i           Omega, OmegaIm, Omega0, kT, alphaT, Prandtl,
     i           TypeLinModel, DoIonDrag,
     i           DimLessStVct, DoNormEigVctZGEEV, qNEV,
C
     i           tnW0, roW0, pW0, uT0,
     i           CvW0, gW0, HDensW0, GammaW0, CsoundW0,
     i           DroW0, DtnW0, DuT0, DHDensW0,
     i           SinDip, CosDip, nu0niW, deniW0,
     i           nu0inW, DiffCoef0W, DiffVelocityW0,
     I           DdeniW0, DDiffVelocityW0,
C
     o           KRCT, EValT, EVctT )
C      ------------------------------------------------------------------------
C      The routine computes the wavenumbers, eigenvalues, and eigenvectors,
C      for each isothermal layer.
C      ------------------------------------------------------------------------
C
C     --- input variables ---
      implicit none
      integer NModes, nzW, nzWT, IdzT, TypeLinModel, qNEV
      integer izT(nzWT+1)
      real*8  Omega, OmegaIm, Omega0, kT, Prandtl
      real*8  alphaT(nzWT)
C
      logical DoIonDrag
      logical DimLessStVct, DoNormEigVctZGEEV
C
      real*8 tnW0(nzW)
      real*8 roW0(nzW)
      real*8 pW0(nzW)
      real*8 uT0(nzW)
      real*8 CvW0(nzW)
      real*8 gW0(nzW)
      real*8 HDensW0(nzW)
      real*8 GammaW0(nzW)
      real*8 CsoundW0(nzW)
c
      real*8 DroW0(nzW)
      real*8 DtnW0(nzW)
      real*8 DuT0(nzW)
      real*8 DHDensW0(nzW)
c
      real*8 SinDip, CosDip
      real*8 nu0niW(nzW)
      real*8 deniW0(nzW)
      real*8 nu0inW(nzW)
      real*8 DiffCoef0W(nzW)
      real*8 DiffVelocityW0(nzW)
      real*8 DdeniW0(nzW)
      real*8 DDiffVelocityW0(nzW)
C
C     --- output variables ---
      complex*16 KRCT(2*NModes, nzWT)
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
C
C     --- local variables ---
      integer    SimpleLinModel, i, j, iT, k, M, IDouble1, IDouble2
      real*8     alpha, normJ
      complex*16 AVRi
      logical    CheckDispPolynom, DoAnalyticEVal, PrintInfo
C
      real*8     g0, Cv, kX, kX2
      real*8     mu0, Dmu0, muK0, ro0, DRo0, u0, Du0, T0, DT0, p0
      real*8     HDens, Gamma, Csound2, HDH
      real*8     nu0ni, nu0in, DiffCoef, uD0, deni, Ddeni, DuD0
      real*8     Dmu0ro0, DT0T
      real*8     DEValR
      real*8     scale
      complex*16 im, one, zero
      complex*16 OmegaC, Omega0C, OmegaDop
      complex*16 FX(6), FZ(6), P(6)
      logical    DoubleEVal
C
      integer, allocatable:: MapJ(:)
      complex*16, allocatable:: KRCN(:), EValN(:), EVctN(:,:)
C
      complex*16, allocatable:: Amat(:,:)
      real*8, allocatable:: resN(:), resA(:)
      real*8, allocatable:: Wreal( : )
C
c     --- set flags ---
      PrintInfo = .false.
      CheckDispPolynom = .true.  !(for testing)
      DoAnalyticEVal   = .false. !can be true only for TypeLinModel = 2 .and
                                 !DoIonDrag = .false. (for testing)
      if ( (TypeLinModel .eq. 1) .or. DoIonDrag )
     &      DoAnalyticEVal = .false.
      SimpleLinModel = 0 !0 - general model (default)
                         !1 - Vadas and Nicolls (2012) model (for testing)
                         !2 - Knight at al. (2024) model (for testing)
C
C     --- allocate ---
      M = 2*NModes
      allocate( Amat(M,M), resN(M), resA(M) )
      allocate( KRCN(M), EValN(M), EVctN(M,M) )
      allocate( MapJ(M) )
c
      kX  = 1.d-3 * kT  !kx in 1/m
      kX2 = kX * kX
C
      OmegaC  = dcmplx( Omega,  OmegaIm ) !Omega
      Omega0C = dcmplx( Omega0, 0.d0 )    !Omega0
C
      im   = ( 0.d0, 1.d0 )
      one  = ( 1.d0, 0.d0 )
      zero = ( 0.d0, 0.d0 )
C     -------------------------------------------------------------------------
C
C
C                             loop over layers
C
C
C     -------------------------------------------------------------------------
      do iT = 1, nzWT
c
c       --- midpoint of each isothermal region ---
        i = izT(iT) + IdzT / 2 !iT = 1 => izT(iT) = 1 => i = 2
                               !iT = 2 => izT(iT) = 3 => i = 4
                               !iT = 3 => izT(iT) = 5 => i = 6
c
        alpha = alphaT(iT)
        T0  = tnW0(i)  !K
        u0  = uT0(i)   !m/s
        ro0 = roW0(i)  !kg/m**3
        p0  = pW0(i)   !N/m**2
        Cv  = CvW0(i)  !J/(kg*K)
        g0  = gW0(i)   !m/s**2
        HDens   = HDensW0(i) !m
        Gamma   = GammaW0(i)
        Csound2 = CsoundW0(i)**2 !(m/s)**2
        Dro0 = DroW0(i)  ![ro]/m
        Du0  = DuT0(i)   ![v]/m
        DT0  = DtnW0(i)  !K/m
        DT0T = DT0 / T0
        HDH  = DHDensW0(i) / HDens
c
        nu0ni = nu0niW(i)
        nu0in = nu0inW(i)
        DiffCoef = DiffCoef0W(i)
        uD0   = DiffVelocityW0(i)
        deni  = deniW0(i)
        Ddeni = DdeniW0(i)
        DuD0  = DDiffVelocityW0(i)
C
C       --- dynamic viscosity mu0  ---
        mu0 = 3.34d-7 * tnW0(i)**0.71
C
C       --- kinematic viscosity  ---
        muK0 = mu0 / ro0
C
C       --- derivatives of viscosity Dmu0  ---
        Dmu0    = 0.71d0 * mu0 * DtnW0(i) / tnW0(i)
        Dmu0ro0 = Dmu0 / ro0
c
        if ( TypeLinModel .eq. 2 ) then
          u0      = 0.d0
          Du0     = 0.d0
          DT0T    = 0.d0
          HDH     = 0.d0 !redundant
          Dmu0ro0 = - muK0 / HDens
        end if
C
C       --- Complex Dopler frequency OmegaDop ---
        OmegaDop = OmegaC - dcmplx(kX * u0,0.D0)
C       -----------------------------------------------------------------------
c                    Ion drag terms for b = (-CosDip, 0.0, -SinDip )
c       -----------------------------------------------------------------------
        if ( DoIonDrag ) then
          call IonDragForcesAndHeating
     i       ( DimLessStVct, kX, u0, Du0, T0, DT0T, CosDip, SinDip,
     i         g0, Cv, nu0in, DiffCoef, uD0, HDens, deni, Ddeni, DuD0,
     i         OmegaC, Omega0C, OmegaDop,
     o         FX, FZ, P )
        else
          do k = 1, 6
            FX(k) = zero
            FZ(k) = zero
            P(k)  = zero
          end do
        end if
C       -----------------------------------------------------------------------
C                               System matrix
c       -----------------------------------------------------------------------
        if ( SimpleLinModel .eq. 0 ) then
          call SYSTEMMatrix
     i       ( M, DimLessStVct, DoNormEigVctZGEEV,
     i         kX, kX2, Prandtl,
     i         T0, ro0, Cv,
     i         HDens, Gamma, Csound2,
     i         Du0, DT0T, HDH,
     i         mu0, muK0, Dmu0ro0,
     i         Omega0C, OmegaDop,
     i         nu0ni, FX, FZ, P,
     o         EvalN, EVctN, Amat )
        else if ( SimpleLinModel .eq. 1 ) then
          call SYSTEMMatrixVN
     i       ( M, DimLessStVct, DoNormEigVctZGEEV,
     i         kX, kX2, Prandtl,
     i         HDens, Gamma, Csound2,
     i         Du0, DT0T, HDH,
     i         muK0, Dmu0ro0,
     i         Omega0C, OmegaDop,
     i         nu0ni, FX, FZ, P,
     o         EvalN, EVctN, Amat )
        else if ( SimpleLinModel .eq. 2 ) then
          call SYSTEMMatrixK
     i       ( M, DimLessStVct, DoNormEigVctZGEEV,
     i         kX, kX2, Prandtl,
     i         HDens, Gamma, Csound2,
     i         Du0, DT0T, HDH,
     i         muK0,
     i         Omega0C, OmegaDop,
     i         nu0ni, FX, FZ, P,
     o         EvalN, EVctN, Amat )
        end if
C       -----------------------------------------------------------------------
c           Check if the absolute value of the dispersion polynomial > 1.e-10
C       -----------------------------------------------------------------------
       if ( CheckDispPolynom ) then
         call CheckDispersionPolynomial
     i      ( iT, M, EValN, OmegaDop, Csound2, Gamma, muK0,
     i        kX, HDens, Prandtl )
       end if
C      -----------------------------------------------------------------------
c          Recompute analytical eigenvalues EvalN for TypeLinModel = 2
c                  (the eigenvectors EVctN remain unchanged)
C      -----------------------------------------------------------------------
       if ( (TypeLinModel .eq. 2) .and. .not. DoIonDrag .and.
     &       DoAnalyticEVal ) then
         call DispersionEquation
     i      ( M, alpha, Omega0, T0, p0, g0, kT, HDens,
     i        mu0, Gamma, Cv, Prandtl,
     o        EvalN )
       end if
C      -----------------------------------------------------------------------
c       Normalize eigenvectors by EVctN(qNEV,j) corresponding to qNEV
C      -----------------------------------------------------------------------
        if ( .not. DoNormEigVctZGEEV ) then
          do j = 1, M
            scale = 1.d0 / abs( EVctN(qNEV,j) )
            do i = 1, M
              EVctN(i,j) = scale * EVctN(i,j)
            end do
          end do
        end if
C
C       --- wavenumbers in 1/km ---
        do j = 1, 6
          KRCN(j) = im * ( EValN(j) - 0.5d0*alpha ) * kT
        end do
C       -----------------------------------------------------------------------
C                          Identify double eigenvalue
C       -----------------------------------------------------------------------
        DoubleEVal = .false.
        do i = 1, M - 1
          do j = i+1, M
            DEValR = abs( EValN(i) - EValN(j) )
     &             / max(abs(EValN(i)), abs(EValN(j)) )
            if ( DEValR .lt. 1.d-6 ) then
              DoubleEVal = .true.
              IDouble1   = i
              IDouble2   = j
              goto 1
            end if
          end do
        end do
1       continue
        if ( DoubleEVal ) then
          write(6,'(a)') 'Eigenvalues:'
          do i = 1, M
            write(6,'(I2,1X,"(",1PE11.4,",",1PE11.4,")")') i, EValN(i)
          end do
          write(6,'(I2,1X,I2,1X,"(",1PE11.4,",",1PE11.4,")")')
     &          iT, IDouble1, EValN(IDouble1)
          write(6,'(I2,1X,I2,1X,"(",1PE11.4,",",1PE11.4,")")')
     &          iT, IDouble2, EValN(IDouble2)
          write(6,*)
        end if
C       -----------------------------------------------------------------------
C                 Residual of the eigenvalue equation (for checking)
C       -----------------------------------------------------------------------
        if ( PrintInfo ) then
          do j = 1, M
            resN(j) = 0.d0
            do i = 1, M
              AVRi = zero
              do k = 1, M
                AVRi = AVRi + Amat(i,k)*EVctN(k,j)
              end do
              resN(j) = resN(j) + abs( AVRi - EValN(j) * EVctN(i,j) )**2
            end do
            normJ = 0.d0
            do i = 1, M
              normJ = normJ + abs(EVctN(i,j))**2
            end do
            resN(j) = sqrt( resN(j) / normJ )
          end do
C
C         --- write results to screen ---
          do j = 1, M
            WRITE(6,'(A,I5,A,I2,A,1PE11.4)')
     &     'Altitude index = ', iT,
     &     ', Vector component = ', j,
     &     ', Residual of Eigenvalue Equation = ', resN(j)
          end do
          write(6,*)
          if ( mod(iT,50) .eq. 0 ) then
            write(6,'(a)') 'Press Enter to continue'
            read(5,*)
          end if
        end if
C       -----------------------------------------------------------------------
C                     Order the solutions according to
C            - the real part of the eigenvalues EVal(6), or
C            - the imaginary part of the eigenvectors KR(6)
C
C       In the latter case, we have:
C         1. Ascending Acoustic-Gravity waves AAGW, with for example Im(KR) = -1,
C         2. Ascending Viscosity waves AVW, with for example Im(KR) = -5,
C         3. Ascending Thermal Conduction waves ATW, with for example Im(KR) = -10,
C
C         4. Descending Acoustic-Gravity waves DAGW, with for example Im(KR) = 1,
C         5. Descending Viscosity waves DVW, with for example Im(KR) = 5,
C         6. Descending Thermal Conduction waves DTW, with for example Im(KR) = 10,
C       ------------------------------------------------------------------------
        allocate( Wreal(M) )
        do i = 1, M
          Wreal(i) = dble( EvalN(i) )
        end do
        call SortEigenValEigenVct( M, M, Wreal, KRCN, EValN, EVctN )
C       -----------------------------------------------------------------------
C       After sorting, we have:
C         Im(KR(1)) < Im(KR(2)) < Im(KR(3)) < Im(KR(4)) < Im(KR(5)) < Im(KR(6))
C           AscTW       AsCvW       AscGW      DescGW      DesCvW      DescTW
C
C       We define the solution vector KRT by
C              KRT(1) = KR(3) = ascending  gravity waves  GW,
C              KRT(2) = KR(2) = ascending  viscosity waves VW,
C              KRT(3) = KR(1) = ascending  thermal TW,
C              KRT(4) = KR(4) = descending gravity waves GW,
C              KRT(5) = KR(5) = descending viscosity waves VW,
C              KRT(6) = KR(6) = descending thermal TW
C       Thus, the pairs of waves are
C                Wave          Ascending     Descending
C            --------------------------------------------
C             gravity GW        KRT(1)        KRT(4)
C             viscosity VW      KRT(2)        KRT(5)
C             thermal TW        KRT(3)        KRT(6)
C       For the solution vector. we have
C       Im(KRT(3)) < Im(KRT(2)) < Im(KRT(1)) < Im(KRT(4)) < Im(KRT(5)) < Im(KRT(6))
C          AscTW       AsCvW        AscGW        DescGW       DesCvW       DescTW
C       -----------------------------------------------------------------------
        MapJ(1) = 3
        MapJ(2) = 2
        MapJ(3) = 1

        MapJ(4) = 4
        MapJ(5) = 5
        MapJ(6) = 6
C
        do j = 1, M
          KRCT(j,iT)  = KRCN(MapJ(j))
          EValT(j,iT) = EValN(MapJ(j))
          do i = 1, M
            EVctT(i,j,iT) = EVctN(i,MapJ(j))
          end do
        end do
C
C       --- deallocate Lapack ---
        deallocate( Wreal )
      end do !end layers
C
C     --- deallocate ---
      deallocate( Amat, resN, resA )
      deallocate( KRCN, EValN, EVctN )
      deallocate( MapJ )
      end
c     *************************************************************************
      subroutine SYSTEMMatrix
     i         ( M, DimLessStVct, DoNormEigVctZGEEV,
     i           kX, kX2, Prandtl,
     i           T0, ro0, Cv,
     i           HDens, Gamma, Csound2,
     i           Du0, DT0T, HDH,
     i           mu0, muK0, Dmu0ro0,
     i           Omega0C, OmegaDop,
     i           nu0ni, FX, FZ, P,
     o           EvalN, EVctN, Amat )
      implicit none
c
C      -- arguments ---
      integer    M
      real*8     kX, kX2, Prandtl
      real*8     T0, ro0, Cv
      real*8     HDens, Gamma, Csound2
      real*8     Du0, DT0T, HDH
      real*8     mu0, muK0, Dmu0ro0, nu0ni
      complex*16 Omega0C, OmegaDop
      complex*16 FX(6), FZ(6), P(6)
      logical    DimLessStVct, DoNormEigVctZGEEV
c
      complex*16 EValN(M), EVctN(M,M)
      complex*16 Amat(6,6)
c
C     --- local scalars ---
      integer    i, j, k
      real*8     GammaPr, fctu2, fctR, FctDimLessStVct
      complex*16 im, one, zero
      complex*16 GmOmDop, GmOm0, GmOmDop2
      complex*16 fct, fct0
c
      integer    INFO, LDVL, LDVR, LWORK
      complex*16, allocatable:: AmatT(:,:)
      real*8, allocatable:: RWORK( : )
      complex*16, allocatable:: VL( :, : ), VR( :, : ),
     $                          W( : ), WORK( : )
      real*8, allocatable:: Wreal( : )

      im   = ( 0.0D0, 1.0D0 )
      one  = ( 1.0D0, 0.0D0 )
      zero = ( 0.0D0, 0.0D0 )
c
      GmOmDop  = Gamma * OmegaDop
      GmOmDop2 = GmOmDop * OmegaDop
      GmOm0    = Gamma * Omega0C
      GammaPr  = Gamma / Prandtl
c
      fctu2 = mu0 * Du0**2 / (Cv * T0 * ro0)
c
c     --- factor for obtaining a dimensionless system of equations ---
      if ( DimLessStVct ) then
        FctDimLessStVct = kX
      else
        FctDimLessStVct = 1.d0
      end if
C
C     --- initialize Amat ---
      Amat = zero
C
C     --- Equation 1 ---
      Amat(1,4) = FctDimLessStVct / kX
C
C     --- Equation 2 ---
      Amat(2,5) = FctDimLessStVct / KX
C
C     --- Equation 3 ---
      Amat(3,6) = FctDimLessStVct / kX
C     -----------------------------------------------------------------------
C                           Equation 4
c     -----------------------------------------------------------------------
      fct0 = one / (muK0 * kX * FctDimLessStVct)
      fct  = Dmu0ro0 * Du0 - im * kX * Csound2 / Gamma
      Amat(4,1) = im * OmegaDop + 4.d0 * kX2 * muK0 / 3.d0
     &          + kX * fct / OmegaDop

      Amat(4,2) = Du0 + im * kX * Dmu0ro0
     &          - im * fct / ( OmegaDop * HDens )

      Amat(4,3) = - ( im * Csound2 * kX2 /Gamma
     &          + 0.71d0 * kX * Du0 * Dmu0ro0 ) / Omega0C

      Amat(4,4) = - Dmu0ro0
      Amat(4,4) = FctDimLessStVct * Amat(4,4)

      Amat(4,5) = im * kX * muK0 / 3.d0 + im * fct / OmegaDop
      Amat(4,5) = FctDimLessStVct * Amat(4,5)

      Amat(4,6) = - 0.71d0 * kX * mu0 * Du0 / ( ro0 * Omega0C )
      Amat(4,6) = FctDimLessStVct * Amat(4,6)
      do k = 1, 6
        Amat(4,k) = fct0 * ( Amat(4,k) + nu0ni * FX(k) )
      end do
C     -----------------------------------------------------------------------
C                               Equation 5
c     -----------------------------------------------------------------------
      fct  = 4.d0 * muK0 / 3.d0 - im * Csound2 / ( Gamma * OmegaDop )
      fct0 = one / (fct * kX * FctDimLessStVct)
      Amat(5,1) = - 2.d0 * im * kX * Dmu0ro0 / 3.d0
     &          + Csound2 * kX2 * Du0 / GmOmDop2

      Amat(5,2) = im * OmegaDop + kX2 * muK0
     &          - im * Csound2 * ( kX * Du0 / OmegaDop - HDH )
     &          / (GmOmDop * HDens)

      Amat(5,3) = - Csound2 * kX * (1.d0 / HDens - DT0T)
     &          / GmOm0

      Amat(5,4) = kX * ( im * muK0 / 3.d0 + Csound2 / GmOmDop )
      Amat(5,4) = FctDimLessStVct * Amat(5,4)

      Amat(5,5) = im * Csound2 * ( kx * Du0/OmegaDop - 1.d0 / HDens )
     &          / GmOmDop - 4.d0 * Dmu0ro0 / 3.d0
      Amat(5,5) = FctDimLessStVct * Amat(5,5)

      Amat(5,6) = kX * Csound2 / GmOm0
      Amat(5,6) = FctDimLessStVct * Amat(5,6)
      do k = 1, 6
        Amat(5,k) = fct0 * ( Amat(5,k) + nu0ni * FZ(k) )
      end do
C     -----------------------------------------------------------------------
C                               Equation 6
c     -----------------------------------------------------------------------
      fctR = muK0 * GammaPr
      fct0 = one / (fctR * kX * FctDimLessStVct)
      Amat(6,1) = ( - im * (Gamma - 1.d0) + fctu2 / OmegaDop ) * Omega0C

      Amat(6,2) = Omega0C * ( DT0T + im * 2.d0 * muK0 * kX * Du0
     &          /(Cv * T0) - im * fctu2 / (OmegaDop * HDens ) ) / kX

      Amat(6,3) = im * OmegaDop + muK0 * GammaPr * kX2
     &          - 1.71d0 * GammaPr * Dmu0ro0 * DT0T - 0.71d0 * fctu2

      Amat(6,4) = - 2.d0 * muK0 * Omega0C * Du0 / (Cv * T0 * kX)
      Amat(6,4) = FctDimLessStVct * Amat(6,4)

      Amat(6,5) = ( (Gamma - 1.d0) + im * fctu2 / OmegaDop )
     &          * Omega0C / kX
      Amat(6,5) = FctDimLessStVct * Amat(6,5)

      Amat(6,6) = - GammaPr * ( Dmu0ro0 + 2.71d0 * muK0 * DT0T )
      Amat(6,6) = FctDimLessStVct * Amat(6,6)
      do k = 1, 6
        Amat(6,k) = fct0 * ( Amat(6,k) + nu0ni * P(k) )
      end do
C     -----------------------------------------------------------------------
C                        Eigenvalues and Eigenvectors
c     -----------------------------------------------------------------------
      allocate( AmatT(M,M) )
      do i = 1, M
        do j = 1, M
          AmatT(i,j) = Amat(i,j)
        end do
      end do
C
      LDVL  = 1
      LDVR  = M
      LWORK = 20*M !2*M
      allocate( W(M) )
      allocate( VL( LDVL, M ), VR( LDVR, M ) )
      allocate( RWORK(2*M) )
      allocate( WORK( LWORK ) )
      W  = zero
      VR = zero
      call ZGEEV0( 'NL', 'VR', DoNormEigVctZGEEV,
     &     M, AmatT, M, W, VL, LDVL, VR, LDVR,
     &     WORK, LWORK, RWORK, INFO )
C
C     --- sort eigenvalues with their realparts in ascending order ---
      allocate( Wreal(M) )
      do i = 1, M
        Wreal(i) = dble( W(i) )
      end do
      call SortEigenvalues ( M, M, Wreal, W, VR )
C
C     --- eigenvalues and eigenvectors ---
      do j = 1, M
        EValN(j) = W(j)
        do i = 1, M
          EVctN(i,j) = VR(i,j)
        end do
      end do
c
      deallocate( AmatT, W, VL, VR, RWORK, WORK, Wreal )
      end
c     *************************************************************************
      subroutine SYSTEMMatrixVN
     i         ( M, DimLessStVct, DoNormEigVctZGEEV,
     i           kX, kX2, Prandtl,
     i           HDens, Gamma, Csound2,
     i           Du0, DT0T, HDH,
     i           muK0, Dmu0ro0,
     i           Omega0C, OmegaDop,
     i           nu0ni, FX, FZ, P,
     o           EvalN, EVctN, Amat )
      implicit none
c
C      -- arguments ---
      integer    M
      real*8     kX, kX2, Prandtl
      real*8     HDens, Gamma, Csound2
      real*8     Du0, DT0T, HDH
      real*8     muK0, Dmu0ro0, nu0ni
      complex*16 Omega0C, OmegaDop
      complex*16 FX(6), FZ(6), P(6)
      logical    DimLessStVct, DoNormEigVctZGEEV
c
      complex*16 EValN(M), EVctN(M,M)
      complex*16 Amat(6,6)
c
C     --- local scalars ---
      integer    i, j, k
      real*8     GammaPr, fctR, FctDimLessStVct
      complex*16 im, one, zero
      complex*16 GmOmDop, GmOm0, GmOmDop2
      complex*16 fct, fct0
c
      integer    INFO, LDVL, LDVR, LWORK
      complex*16, allocatable:: AmatT(:,:)
      real*8, allocatable:: RWORK( : )
      complex*16, allocatable:: VL( :, : ), VR( :, : ),
     $                          W( : ), WORK( : )
      real*8, allocatable:: Wreal( : )

      im   = ( 0.0D0, 1.0D0 )
      one  = ( 1.0D0, 0.0D0 )
      zero = ( 0.0D0, 0.0D0 )
c
      GmOmDop  = Gamma * OmegaDop
      GmOmDop2 = GmOmDop * OmegaDop
      GmOm0    = Gamma * Omega0C
      GammaPr  = Gamma / Prandtl
c
c     --- factor for obtaining a dimensionless system of equations ---
      if ( DimLessStVct ) then
        FctDimLessStVct = kX
      else
        FctDimLessStVct = 1.d0
      end if
C
C     --- initialize Amat ---
      Amat = zero
C
C     --- Equation 1 ---
      Amat(1,4) = FctDimLessStVct / kX
C
C     --- Equation 2 ---
      Amat(2,5) = FctDimLessStVct / kX
C
C     --- Equation 3 ---
      Amat(3,6) = FctDimLessStVct / kX
C     -----------------------------------------------------------------------
C                           Equation 4
c     -----------------------------------------------------------------------
      fct0 = one / (muK0 * kX * FctDimLessStVct)
      fct  = - im * kX * Csound2 / Gamma
      Amat(4,1) = im * OmegaDop + 4.d0 * kX2 * muK0 / 3.d0
     &          + kX * fct / OmegaDop

      Amat(4,2) = Du0 + im * kX * Dmu0ro0
     &          - im * fct / ( OmegaDop * HDens )

      Amat(4,3) = - im * Csound2 * kX2 / GmOm0

      Amat(4,4) = - Dmu0ro0
      Amat(4,4) = FctDimLessStVct * Amat(4,4)

      Amat(4,5) = im * kX * muK0 / 3.d0 + im * fct / OmegaDop
      Amat(4,5) = FctDimLessStVct * Amat(4,5)

      Amat(4,6) = zero
      do k = 1, 6
        Amat(4,k) = fct0 * ( Amat(4,k) + nu0ni * FX(k) )
      end do
C     -----------------------------------------------------------------------
C                               Equation 5
c     -----------------------------------------------------------------------
      fct  = 4.d0 * muK0 / 3.d0 - im * Csound2 / ( Gamma * OmegaDop )
      fct0 = one / (fct * kX * FctDimLessStVct)
      Amat(5,1) = - 2.d0 * im * kX * Dmu0ro0 / 3.d0
     &          + Csound2 * kX2 * Du0 / GmOmDop2

      Amat(5,2) = im * OmegaDop + kX2 * muK0
     &          - im * Csound2 * ( kX * Du0 / OmegaDop - HDH )
     &          / (GmOmDop * HDens)

      Amat(5,3) = - Csound2 * kX * (1.d0 / HDens - DT0T)
     &          / GmOm0

      Amat(5,4) = kX * ( im * muK0 / 3.d0 + Csound2 / GmOmDop )
      Amat(5,4) = FctDimLessStVct * Amat(5,4)

      Amat(5,5) = im * Csound2 * ( kx * Du0/OmegaDop - 1.d0 / HDens )
     &          / GmOmDop - 4.d0 * Dmu0ro0 / 3.d0
      Amat(5,5) = FctDimLessStVct * Amat(5,5)

      Amat(5,6) = kX * Csound2 / GmOm0
      Amat(5,6) = FctDimLessStVct * Amat(5,6)
      do k = 1, 6
        Amat(5,k) = fct0 * ( Amat(5,k) + nu0ni * FZ(k) )
      end do
C     -----------------------------------------------------------------------
C                               Equation 6
c     -----------------------------------------------------------------------
      fctR = muK0 * GammaPr
      fct0 = one / (fctR * kX * FctDimLessStVct)
      Amat(6,1) = - im * (Gamma - 1.d0) * Omega0C

      Amat(6,2) = Omega0C * DT0T / kX

      Amat(6,3) = im * OmegaDop + muK0 * GammaPr * kX2
     &          - GammaPr * Dmu0ro0 * DT0T

      Amat(6,4) = zero

      Amat(6,5) = (Gamma - 1.d0) * Omega0C / kX
      Amat(6,5) = FctDimLessStVct * Amat(6,5)

      Amat(6,6) = - GammaPr * ( Dmu0ro0 + 2.d0 * muK0 * DT0T )
      Amat(6,6) = FctDimLessStVct * Amat(6,6)
      do k = 1, 6
        Amat(6,k) = fct0 * ( Amat(6,k) + nu0ni * P(k) )
      end do
C     -----------------------------------------------------------------------
C                        Eigenvalues and Eigenvectors
c     -----------------------------------------------------------------------
      allocate( AmatT(M,M) )
      do i = 1, M
        do j = 1, M
          AmatT(i,j) = Amat(i,j)
        end do
      end do
C
      LDVL  = 1
      LDVR  = M
      LWORK = 20*M !2*M
      allocate( W(M) )
      allocate( VL( LDVL, M ), VR( LDVR, M ) )
      allocate( RWORK(2*M) )
      allocate( WORK( LWORK ) )
      W  = zero
      VR = zero
      call ZGEEV0( 'NL', 'VR', DoNormEigVctZGEEV,
     &     M, AmatT, M, W, VL, LDVL, VR, LDVR,
     &     WORK, LWORK, RWORK, INFO )
C
C     --- sort eigenvalues with their realparts in ascending order ---
      allocate( Wreal(M) )
      do i = 1, M
        Wreal(i) = dble( W(i) )
      end do
      call SortEigenvalues ( M, M, Wreal, W, VR )
C
C     --- eigenvalues and eigenvectors ---
      do j = 1, M
        EValN(j) = W(j)
        do i = 1, M
          EVctN(i,j) = VR(i,j)
        end do
      end do
c
      deallocate( AmatT, W, VL, VR, RWORK, WORK, Wreal )
      end
c     *************************************************************************
      subroutine SYSTEMMatrixK
     i         ( M, DimLessStVct, DoNormEigVctZGEEV,
     i           kX, kX2, Prandtl,
     i           HDens, Gamma, Csound2,
     i           Du0, DT0T, HDH,
     i           muK0,
     i           Omega0C, OmegaDop,
     i           nu0ni, FX, FZ, P,
     o           EvalN, EVctN, Amat )
      implicit none
c
C      -- arguments ---
      integer    M
      real*8     kX, kX2, Prandtl
      real*8     HDens, Gamma, Csound2
      real*8     Du0, DT0T, HDH
      real*8     muK0, nu0ni
      complex*16 Omega0C, OmegaDop
      complex*16 FX(6), FZ(6), P(6)
      logical    DimLessStVct, DoNormEigVctZGEEV
c
      complex*16 EValN(M), EVctN(M,M)
      complex*16 Amat(6,6)
c
C     --- local scalars ---
      integer    i, j, k
      real*8     GammaPr, fctR, FctDimLessStVct
      complex*16 im, one, zero
      complex*16 GmOmDop, GmOm0, GmOmDop2
      complex*16 fct, fct0
c
      integer    INFO, LDVL, LDVR, LWORK
      complex*16, allocatable:: AmatT(:,:)
      real*8, allocatable:: RWORK( : )
      complex*16, allocatable:: VL( :, : ), VR( :, : ),
     $                          W( : ), WORK( : )
      real*8, allocatable:: Wreal( : )

      im   = ( 0.0D0, 1.0D0 )
      one  = ( 1.0D0, 0.0D0 )
      zero = ( 0.0D0, 0.0D0 )
c
      GmOmDop  = Gamma * OmegaDop
      GmOmDop2 = GmOmDop * OmegaDop
      GmOm0    = Gamma * Omega0C
      GammaPr  = Gamma / Prandtl
c
c     --- factor for obtaining a dimensionless system of equations ---
      if ( DimLessStVct ) then
        FctDimLessStVct = kX
      else
        FctDimLessStVct = 1.d0
      end if
C
C     --- initialize Amat ---
      Amat = zero
C
C     --- Equation 1 ---
      Amat(1,4) = FctDimLessStVct / kX
C
C     --- Equation 2 ---
      Amat(2,5) = FctDimLessStVct / kX
C
C     --- Equation 3 ---
      Amat(3,6) = FctDimLessStVct / kX
C     -----------------------------------------------------------------------
C                           Equation 4
c     -----------------------------------------------------------------------
      fct0 = one / (muK0 * kX * FctDimLessStVct)
      fct  = - im * kX * Csound2 / Gamma
      Amat(4,1) = im * OmegaDop + 4.d0 * kX2 * muK0 / 3.d0
     &          + kX * fct / OmegaDop

      Amat(4,2) = Du0 - im * fct / ( OmegaDop * HDens )

      Amat(4,3) = - im * Csound2 * kX2 / GmOm0

      Amat(4,4) = zero

      Amat(4,5) = im * kX * muK0 / 3.d0 + im * fct / OmegaDop
      Amat(4,5) = FctDimLessStVct * Amat(4,5)

      Amat(4,6) = zero
      do k = 1, 6
        Amat(4,k) = fct0 * ( Amat(4,k) + nu0ni * FX(k) )
      end do
C     -----------------------------------------------------------------------
C                               Equation 5
c     -----------------------------------------------------------------------
      fct  = 4.d0 * muK0 / 3.d0 - im * Csound2 / ( Gamma * OmegaDop )
      fct0 = one / (fct * kX * FctDimLessStVct)
      Amat(5,1) = Csound2 * kX2 * Du0 / GmOmDop2

      Amat(5,2) = im * OmegaDop + kX2 * muK0
     &          - im * Csound2 * ( kX * Du0 / OmegaDop - HDH )
     &          / (GmOmDop * HDens)

      Amat(5,3) = - Csound2 * kX * (1.d0 / HDens - DT0T)
     &          / GmOm0

      Amat(5,4) = kX * ( im * muK0 / 3.d0 + Csound2 / GmOmDop )
      Amat(5,4) = FctDimLessStVct * Amat(5,4)

      Amat(5,5) = im * Csound2 * ( kx * Du0/OmegaDop - 1.d0 / HDens )
     &          / GmOmDop
      Amat(5,5) = FctDimLessStVct * Amat(5,5)

      Amat(5,6) = kX * Csound2 / GmOm0
      Amat(5,6) = FctDimLessStVct * Amat(5,6)
      do k = 1, 6
        Amat(5,k) = fct0 * ( Amat(5,k) + nu0ni * FZ(k) )
      end do
C     -----------------------------------------------------------------------
C                               Equation 6
c     -----------------------------------------------------------------------
      fctR = muK0 * GammaPr
      fct0 = one / (fctR * kX * FctDimLessStVct)
      Amat(6,1) = - im * (Gamma - 1.d0) * Omega0C

      Amat(6,2) = Omega0C * DT0T / kX

      Amat(6,3) = im * OmegaDop + muK0 * GammaPr * kX2

      Amat(6,4) = zero

      Amat(6,5) = (Gamma - 1.d0) * Omega0C / kX
      Amat(6,5) = FctDimLessStVct * Amat(6,5)

      Amat(6,6) = zero
      do k = 1, 6
        Amat(6,k) = fct0 * ( Amat(6,k) + nu0ni * P(k) )
      end do
C     -----------------------------------------------------------------------
C                        Eigenvalues and Eigenvectors
c     -----------------------------------------------------------------------
      allocate( AmatT(M,M) )
      do i = 1, M
        do j = 1, M
          AmatT(i,j) = Amat(i,j)
        end do
      end do
C
      LDVL  = 1
      LDVR  = M
      LWORK = 20*M !2*M
      allocate( W(M) )
      allocate( VL( LDVL, M ), VR( LDVR, M ) )
      allocate( RWORK(2*M) )
      allocate( WORK( LWORK ) )
      W  = zero
      VR = zero
      call ZGEEV0( 'NL', 'VR', DoNormEigVctZGEEV,
     &     M, AmatT, M, W, VL, LDVL, VR, LDVR,
     &     WORK, LWORK, RWORK, INFO )
C
C     --- sort eigenvalues with their realparts in ascending order ---
      allocate( Wreal(M) )
      do i = 1, M
        Wreal(i) = dble( W(i) )
      end do
      call SortEigenvalues ( M, M, Wreal, W, VR )
C
C     --- eigenvalues and eigenvectors ---
      do j = 1, M
        EValN(j) = W(j)
        do i = 1, M
          EVctN(i,j) = VR(i,j)
        end do
      end do
c
      deallocate( AmatT, W, VL, VR, RWORK, WORK, Wreal )
      end
C     *************************************************************************
C
c
C                              Dispersion Equation
c      for an isothermal, homogeneous, and windless atmosphere without ion drag
C
C
C
C     ****************************************·*********************************
      subroutine DispersionEquation
     i         ( M, alpha, Omega0, T0, p0, g0, kT, HDens,
     i           mu0, Gamma, Cv, Prandtl,
     o           EValA )
      implicit none
      integer M
      real*8  alpha, Omega0, T0, p0, g0, kT, HDens
      real*8  mu0, Gamma, Cv, Prandtl
c
      complex*16 EValA(M)
c
c     --- local variables ---
      integer     j
      real*8      kX, kX2, beta, lam, H
      complex*16  im, etaC, niuC
      complex*16  C0C, C1C, C2C, C3C, R1C, R2C, R3C, kRC
      complex*16, allocatable:: KRCA(:)

      im = ( 0.d0, 1.d0 )
      allocate( KRCA(M) )
c
      kX   = 1.d-3 * kT
      kX2  = kX * kX
c
      beta = Omega0**2 / (g0 * KX2 * HDens)
c
      etaC = im * Omega0 * mu0 / ( 3.d0 * p0 )
c
      lam  = Gamma * Cv * mu0 / Prandtl
      niuC = im * T0 * kX2 * lam / ( Omega0 * p0 )
c
      C3C = - 3.d0 * etaC * niuC * (1.d0 + 4.d0*etaC)

      C2C = 3.d0*etaC * (1.d0 + 4.d0*etaC) / (gamma - 1.d0)
     &  + niuC * beta * (1.d0 + 7.d0*etaC) + 3.d0 * etaC

      C1C = - (beta**2 - 2.d0*etaC * alpha**2
     &  * (1.d0 + 3.d0*etaC))*niuC
     &  - beta * (1.d0 + 7.d0*etaC) / (gamma - 1.d0) - beta

      C0C = (beta**2 - 2.d0*etaC * alpha**2 * (1.d0 + 3.d0*etaC))
     &  / (gamma - 1.d0) + alpha**2 * (1.d0 + 3.d0*etaC)
c
      call SolveCubicEquationDouble
     i    ( C0C, C1C, C2C, C3C,
     o      R1C, R2C, R3C )
c
c     --- wavenumbers ---
      H = 1.d-3 * HDens !in km
      kRC = sqrt( (R1C - 1.d0) * kT**2
     &    - 1.d0 / (2.d0*H)**2 ) !in 1/km
      KRCA(1) = kRC
      KRCA(2) = -kRC

      kRC = sqrt( (R2C - 1.d0) * kT**2
     &    - 1.d0 / (2.d0*H)**2 ) !in 1/km
      KRCA(3) = kRC
      KRCA(4) = -kRC

      kRC = sqrt( (R3C - 1.d0) * kT**2
     &    - 1.d0 / (2.d0*H)**2 ) !in 1/km
      KRCA(5) = -kRC
      KRCA(6) = kRC
c
c     --- Eigenvalues ---
      do j = 1, 6
        EValA(j) = - im * KRCA(j) / kT + 0.5d0*alpha
      end do
c
      deallocate( KRCA )
      end
c     ********************************************************************
      subroutine SolveCubicEquationDouble
     i         ( C0C, C1C, C2C, C3C,
     o           R1C, R2C, R3C )
      implicit none
      complex*16    C0C, C1C, C2C, C3C
      complex*16    R1C, R2C, R3C
c
      real*8     pow13
      complex*16 im, p, q, r, a1, a2, a, b
      real*8     res1, res2, res3
      logical    DoResidual
c
      DoResidual = .false.
!
      im    = (0.d0, 1.d0)
      pow13 = 1.d0 / 3.d0
c
      p = C2C / C3C
      q = C1C / C3C
      r = C0C / C3C
c
      a1 = - 2.d0 * p**3 + 9.d0 * p * q - 27.d0 * r
      a2 = - p**2 * q**2 + 4.d0 * q**3 + 4.d0 * p**3 * r
     &     - 18.d0 * p * q * r + 27.d0 * r**2
      a  = ( a1 + 3.d0*sqrt(3.d0 * a2) )**pow13
     &   / ( 3.d0 * 2.d0**pow13 )
      b  = ( - p**2 + 3.d0 * q ) / ( 9.d0 * a )
c
      a1  = - 0.5d0 * ( 1.d0 + im * sqrt(3.d0) )
      a2  = - 0.5d0 * ( 1.d0 - im * sqrt(3.d0) )
      R1C = - p / 3.0 + a - b
      R2C = - p / 3.0 + a1 * a - a2 * b
      R3C = - p / 3.0 + a2 * a - a1 * b
c
      if ( DoResidual ) then
        res1 = abs(C3C*R1C**3 + C2C*R1C**2 + C1C*R1C + C0C )
        res2 = abs(C3C*R2C**3 + C2C*R2C**2 + C1C*R2C + C0C )
        res3 = abs(C3C*R3C**3 + C2C*R3C**2 + C1C*R3C + C0C )
        write(6,*) res1, res2, res3
      end if
      end
C     *************************************************************************
C
c
C                          Dispersion Polynomial
C
C
C     ****************************************·*********************************
      subroutine CheckDispersionPolynomial
     i         ( iT, M, EValN, OmegaDop, Csound2, Gamma, mu0,
     i           kX, HDens, Prandtl )
      implicit none
      integer    iT, M
      complex*16 EValN(M), OmegaDop
      real*8     Csound2, Gamma, mu0, kX, HDens, Prandtl
c
      integer    j
      complex*16 EVal, pLam
c
      do j = 1, M
        EVal = EvalN(j)
        call PolynomLambda
     i     ( EVal, OmegaDop, Csound2, Gamma, mu0,
     i       kX, HDens, Prandtl,
     o       pLam )
        if ( abs(pLam) .gt. 1.d-10 ) then
          write(6,'(a,1pe11.4,2(a,i4))')
     &    '- Dispersion polynomial = ', abs(pLam),
     &    ', in layer ', iT, ', for mode = ', j
        end if
      end do
      end
c     *************************************************************************
      subroutine PolynomLambda
     i         ( EVal, OmegaDop, Csound2, Gamma, mu0,
     i           kX, HDens, Prandtl,
     o           pLam )
      implicit none
      complex*16 EVal, OmegaDop
      real*8     Csound2, Gamma, mu0, kX, HDens, Prandtl
c
      complex*16 pLam
c
c     --- locals ---
      complex*16 im, one, lam2, oneMinusLam2
      complex*16 A1, A2, A3, A4
      complex*16 B, term1, term2
      real*8     kX2, term3
c
      im  = (0.0d0, 1.0d0)
      one = (1.0d0, 0.0d0)
c
      kX2 = kX*kX
c
      lam2         = EVal * EVal
      oneMinusLam2 = one - lam2
c
c     --- factors A1..A4 ---
      A1 = - OmegaDop
     &     + im * ( (Gamma*mu0/Prandtl)*kX2 ) * oneMinusLam2

      A2 = - OmegaDop
     &     + im * ( mu0*kX2 ) * oneMinusLam2

      A3 = - OmegaDop
     &     + im * ( (4.0d0*mu0/3.0d0)*kX2 ) * oneMinusLam2

      A4 = - OmegaDop
     &     + im * ( (mu0/Prandtl)*kX2 ) * oneMinusLam2
c
c     --- B = kX^2*(1-lambda^2) + (kX*lambda)/HDens ---
      B  = kX2 * oneMinusLam2 + (kX/HDens) * EVal
c
c     --- assemble p(lambda) ---
      term1 = (OmegaDop/Csound2) * A1 * A2 * A3
      term2 = A2 * A4 * B

      term3 = (kX2*Csound2/(Gamma*Gamma*HDens*HDens))
     &      * (Gamma - 1.0d0)

      pLam = term1 + term2 - term3
      end
C     *************************************************************************
C
C
c
C                                Ion Drag
C
C
C
C     ****************************************·*********************************
      subroutine IonDragForcesAndHeating
     i         ( DimLessStVct, kX, u0, Du0, T0, DT0T, cI, sI,
     i           g0, Cv, nu0in, DiffCoef, uD0, HDens, deni, Ddeni, DuD0,
     i           OmegaC, Omega0C, OmegaDop,
     o           FX, FZ, P )
C     -------------------------------------------------------------------------
C   Step 3: Compute F- and Q-coefficients for f_x, f_z and q.
C
C   A_x = u_d0*cos(I) + u_0*sin^2(I)
C   A_z = u_d0*sin(I) - u_0*cos(I)*sin(I)
C   B   = u_d0^2 + u_0^2*sin^2(I)
C
C   For f_x:
C     F_xu  = sin^2(I)*ω0/kx + A_x*N_u + cos(I)*U_u
C     F_xw  = -sin(I)cos(I)*ω0/kx + A_x*N_w + cos(I)*U_w
C     F_xT  = 0.37*A_x + A_x*N_T + cos(I)*U_T
C     F_xU  = A_x*N_U
C     F_xW  = A_x*N_W + cos(I)*U_W
C     F_xTc = A_x*N_Tc + cos(I)*U_Tc
C
C   For f_z:
C     F_zu  = -sin(I)cos(I)*ω0/kx + A_z*N_u + sin(I)*U_u
C     F_zw  =  cos^2(I)*ω0/kx    + A_z*N_w + sin(I)*U_w
C     F_zT  = 0.37*A_z + A_z*N_T + sin(I)*U_T
C     F_zU  = A_z*N_U
C     F_zW  = A_z*N_W + sin(I)*U_W
C     F_zTc = A_z*N_Tc + sin(I)*U_Tc
C
C   For q:
C     Q_u  = 2*sin^2(I)*ω0/kx*u_0 + B*N_u + 2*u_d0*U_u
C     Q_w  = -2*cos(I)sin(I)*ω0/kx*u_0 + B*N_w + 2*u_d0*U_w
C     Q_T  = 0.37*B + B*N_T + 2*u_d0*U_T
C     Q_U  = B*N_U
C     Q_W  = B*N_W + 2*u_d0*U_W
C     Q_Tc = B*N_Tc + 2*u_d0*U_Tc
C     -------------------------------------------------------------------------
      implicit none
      real*8     kX, u0, Du0, T0, DT0T, cI, sI, g0, Cv
      real*8     nu0in, DiffCoef, uD0, HDens, deni, Ddeni, DuD0
      complex*16 OmegaC, Omega0C, OmegaDop
      logical    DimLessStVct
c
      complex*16 FX(6), FZ(6), P(6)
C
C     --- local variables ---
      integer     i
      real*8      FctDimLessStVct
      real*8      Ax, Az, B
      complex*16  im, one, zero, ScaleF, ScaleP
      complex*16  AxC, AzC, BC
      complex*16  Uu, Uw, UcW, UT, UcT
      complex*16  Nu, Nw, NT, NUc, NWc, NTc
c
      im   = ( 0.d0, 1.d0 )
      one  = ( 1.d0, 0.d0 )
      zero = ( 0.d0, 0.d0 )
c
c     --- factor for obtaining a dimensionless system of equations ---
      if ( DimLessStVct ) then
        FctDimLessStVct = kX
      else
        FctDimLessStVct = 1.d0
      end if
c
      call ComputeUcoeff
     i   ( kX, DT0T, cI, sI, g0,
     i     nu0in, DiffCoef, HDens, deni, Ddeni,
     i     Omega0C, OmegaDop,
     o     Uu, Uw, UcW, UT, UcT )
c
      call ComputeNcoeff
     i   ( kX, u0, Du0, cI, sI,
     i     deni, Ddeni,
     i     OmegaC, Omega0C,
     i     ud0, DuD0,
     i     Uu, Uw, UcW, UT, UcT,
     o     Nu, Nw, NT, NUc, NWc, NTc )
C     -------------------------------------------------------------------------
C                     Abbreviations A_x, A_z, B (all real)
C     -------------------------------------------------------------------------
      Ax = ud0 * cI + u0 * sI * sI
      Az = ud0 * sI - u0 * cI * sI
      B  = ud0 * ud0 + u0 * u0 * sI * sI
C
      AxC = dcmplx( Ax, 0.0d0 )
      AzC = dcmplx( Az, 0.0d0 )
      BC  = dcmplx( B , 0.0d0 )
C     -------------------------------------------------------------------------
C                          Coefficients for f_x
C     -------------------------------------------------------------------------
c
C     --- F_xu = sin^2 I * ω0/kx + A_x*N_u + cosI*U_u
      FX(1) =  Omega0C * (sI*sI)/kX
     &    + AxC * Nu
     &    + cI * Uu
C
C     --- F_xw = -sinI cosI * ω0/kx + A_x*N_w + cosI*U_w
      FX(2) = - Omega0C * (sI*cI)/kX
     &    + AxC * Nw
     &    + cI * Uw
C
C     --- F_xT = 0.37*A_x + A_x*N_T + cosI*U_T
      FX(3) = dcmplx( 0.37d0*Ax, 0.0d0 )
     &    + AxC * NT
     &    + cI * UT
C
C     --- F_xU = A_x * N_U
      FX(4) = AxC * NUc
C
C     --- F_xW = A_x*N_W + cosI*U_W
      FX(5) = AxC * NWc
     &      + cI * UcW
C
C     --- F_xTc = A_x*N_Tc + cosI*U_Tc
      FX(6) = AxC * NTc
     &      + cI * UcT
c
      do i = 4, 6
        FX(i) = FctDimLessStVct * FX(i)
      end do
C     -------------------------------------------------------------------------
C                          Coefficients for f_z
C     -------------------------------------------------------------------------
c
C     --- F_zu = -sinI cosI * ω0/kx + A_z*N_u + sinI*U_u ---
      FZ(1) = - Omega0C * (sI*cI)/kX
     &      + AzC * Nu
     &      + sI * Uu
C
C     --- F_zw = cos^2I * ω0/kx + A_z*N_w + sinI*U_w ---
      FZ(2) =   Omega0C * (cI*cI)/kX
     &      + AzC * Nw
     &      + sI * Uw
C
C     --- F_zT = 0.37*A_z + A_z*N_T + sinI*U_T ---
      FZ(3) = dcmplx( 0.37d0*Az, 0.0d0 )
     &      + AzC * NT
     &      + sI * UT
C
C     --- F_zU = A_z * N_U ---
      FZ(4) = AzC * NUc
C
C     --- F_zW = A_z*N_W + sinI*U_W ---
      FZ(5) = AzC * NWc
     &      + sI * UcW
C
C     --- F_zTc = A_z*N_Tc + sinI*U_Tc ---
      FZ(6) = AzC * NTc
     &      + sI * UcT
c
      do i = 4, 6
        FZ(i) = FctDimLessStVct * FZ(i)
      end do
C     -------------------------------------------------------------------------
C                         Coefficients for q
C     -------------------------------------------------------------------------
c
C     --- Q_u = 2*sin^2I * ω0/kx * u_0 + B*N_u + 2*u_d0*U_u
      P(1) =  Omega0C * 2.d0*u0 * sI*sI/kX
     &     + BC * Nu
     &     + 2.d0*ud0 * Uu
C
C     --- Q_w = -2*cosI sinI * ω0/kx * u_0 + B*N_w + 2*u_d0*U_w
      P(2) = - Omega0C * 2.d0*u0 * cI*sI/kX
     &     + BC * Nw
     &     + 2.d0*ud0 * Uw
C
C     --- Q_T = 0.37*B + B*N_T + 2*u_d0*U_T
      P(3) = dcmplx( 0.37d0*B, 0.0d0 )
     &     + BC * NT
     &     + 2.d0*ud0 * UT
C
C     --- Q_U = B*N_U
      P(4) = BC * NUc
C
C     --- Q_W = B*N_W + 2*u_d0*U_W
      P(5) = BC * NWc
     &     + 2.d0*ud0 * UcW
C
C     --- Q_Tc = B*N_Tc + 2*u_d0*U_Tc
      P(6) = BC * NTc
     &     + 2.d0*ud0 * UcT
c
      do i = 4, 6
        P(i) = FctDimLessStVct * P(i)
      end do
c     -------------------------------------------------------------------------
c                 Apply scaling factors for the underlying equations
c     -------------------------------------------------------------------------
      ScaleF = kX / Omega0C
      ScaleP = one / (Cv * T0)
      do i = 1, 6
        FX(i) = ScaleF * FX(i)
        FZ(i) = ScaleF * FZ(i)
        P(i)  = ScaleP * P(i)
      end do
      end
c     *************************************************************************
      subroutine ComputeUcoeff
     i          ( kX, DT0T, cI, sI, g0,
     i           nu0in, DiffCoef, HDens, deni, Ddeni,
     i           Omega0C, OmegaDop,
     o           Uu, Uw, UcW, UT, UcT )
      implicit none
C
C     --- arguments ---
      real*8      cI, sI, kX
      complex*16  Omega0C, OmegaDop
      real*8      nu0in, deni, Ddeni, DiffCoef
      real*8      DT0T, g0, HDens
      complex*16  Uu, Uw, UcW, UT, UcT
C
C     --- local variables ---
      real*8      P, B, DB
      complex*16  im, ratio, DBc
c
      im   = ( 0.d0, 1.d0 )
C
C     --- P = sin(I) * [ (1/n_i0) dn_i0/dz + (1/T0) dT0/dz ]
      P = sI * ( Ddeni/deni + DT0T )
C
C     --- B = g * sin(I) / nu_in0
      B = g0 * sI / nu0in
C
C     --- DB = D_a0 * P + B   (real)
      DB  = DiffCoef * P + B
C
C     --- ratio = omega_0 / Omega  (complex)
      ratio = Omega0C / OmegaDop
C
C     --- DBc = (D_a0 * P + B) as complex
      DBc = dcmplx( DB, 0.0d0 )
C     -------------------------------------------------------------------------
C     U_u = - (omega_0 / Omega) * (D_a0 * P + B)
C     -------------------------------------------------------------------------
      Uu = - ratio * DBc
C     -------------------------------------------------------------------------
C     U_w = j * (omega_0 / Omega)
C    &      * (1 / (k_x * H_rho)) * (D_a0 * P + B)
C     -------------------------------------------------------------------------
      Uw = im * ratio * DBc / (kX * HDens)
C     -------------------------------------------------------------------------
C     U_W = -j * (omega_0 / Omega)
C    &       * (1 / k_x) * (D_a0 * P + B)
C     -------------------------------------------------------------------------
      UcW = - im * ratio * DBc / kX
C     -------------------------------------------------------------------------
C     U_T = -j * k_x * D_a0 * cos(I)
C    &       + 0.63 * D_a0 * P - 0.37 * B
C     -------------------------------------------------------------------------
      UT = - im * kX * DiffCoef * cI
     &     + 0.63d0 * DiffCoef * P - 0.37d0 * B
C     -------------------------------------------------------------------------
C     U_Tc = D_a0 * sin(I)
C     -------------------------------------------------------------------------
      UcT = DiffCoef * sI
      end
c     *************************************************************************
      subroutine ComputeNcoeff
     i          ( kX, u0, Du0, cI, sI,
     i           deni, Ddeni,
     i           OmegaC, Omega0C,
     i           ud0, DuD0,
     i           Uu, Uw, UcW, UT, UcT,
     o           Nu, Nw, NT, NUc, NWc, NTc )
      implicit none
C
C     --- arguments ---
      real*8      cI, sI, kX
      complex*16  Omega0C, OmegaC
      real*8      deni, Ddeni, u0
      real*8      Du0, ud0, DuD0
      complex*16  Uu, Uw, UcW, UT, UcT
      complex*16  Nu, Nw, NT, NUc, NWc, NTc
C
C     --- local variables ---
      real*8      A         ! A = (1/n_i0) dn_i0/dz * sin(I)
      complex*16  im, N0        ! N_0
      complex*16  term1, term2
      complex*16  fact_u, fact_w
c
      im = ( 0.d0, 1.d0 )
C
C     --- A = (1/n_i0) dn_i0/dz * sin(I)
      A = (Ddeni / deni) * sI
C
C     --- N_0 = j * omega - u_d0 * A
      N0 = im * OmegaC - ud0 * A + Du0 * cI * sI + u0 * A * cI
     &   - DuD0 * sI
C
C     --- Common factors for omega_0 / k_x:
      fact_u = Omega0C * cI / kX   ! cosI * omega_0 / k_x
      fact_w = Omega0C * sI / kX   ! sinI * omega_0 / k_x
C     -------------------------------------------------------------------------
C     N_u
C       term1 = (j k_x cosI - A) * cosI * (omega_0 / k_x)
C       term2 = A * U_u
C     -------------------------------------------------------------------------
      term1 = ( im * kX * cI - A ) * fact_u
      term2 = A * Uu
      Nu = ( term1 + term2 ) / N0
C     -------------------------------------------------------------------------
C     N_w
C       term1 = (j k_x cosI - A) * sinI * (omega_0 / k_x)
C       term2 = A * U_w
C     -------------------------------------------------------------------------
      term1 = ( im * kX * cI - A ) * fact_w
      term2 = A * Uw
      Nw = ( term1 + term2 ) / N0
C     -------------------------------------------------------------------------
C     N_T = (A * U_T) / N_0
C     -------------------------------------------------------------------------
      NT = (A * UT) / N0
C     -------------------------------------------------------------------------
C     N_U = - (cosI sinI * omega_0 / k_x) / N_0
C     -------------------------------------------------------------------------
      NUc = - (fact_u * sI) / N0
C     -------------------------------------------------------------------------
C     N_W = -1/N_0 * ( sin^2 I * (omega_0/k_x) - A * U_W )
C         = ( -sin^2I * omega_0/k_x + A * U_W ) / N_0
C     -------------------------------------------------------------------------
      term1 = - fact_w * sI
      term2 = A * UcW
      NWc   = ( term1 + term2 ) / N0
C     -------------------------------------------------------------------------
C     N_Tcal = (A * U_Tcal) / N_0
C     -------------------------------------------------------------------------
      NTc = ( A * UcT ) / N0
      end
C     *************************************************************************
C
c
C                          Imaginary Frequency Shift
C
C
C
c     *************************************************************************
      subroutine LayerwiseCausalityCondition
     i         ( NFFT, iT, OmegaK,
     i           NModes, nzW, nzWT, IdzT, izT,
     i           AbsOmegaImMin, AbsOmegaImMax, DOmegaIm,
     i           Omega0, kT, alphaT, Prandtl,
     i           TypeLinModel, DoIonDrag,
     i           DimLessStVct, DoNormEigVctZGEEV, qNEV,
C
     i           tnW0, roW0, uT0,
     i           CvW0, gW0, HDensW0, GammaW0, CsoundW0,
     i           DroW0, DtnW0, DuT0, DHDensW0,
     i           SinDip, CosDip, nu0niW, deniW0,
     i           nu0inW, DiffCoef0W, DiffVelocityW0,
     i           DdeniW0, DDiffVelocityW0,
C
     o           fail, OmegaIm, dmin, dmax )
C
C     --- input variables ---
      implicit none
      integer NFFT, iT
      integer NModes, nzW, nzWT, IdzT, TypeLinModel, qNEV
      integer izT(nzWT+1)
      real*8  OmegaK( NFFT )
      real*8  AbsOmegaImMin, AbsOmegaImMax, DOmegaIm
      real*8  Omega0, kT, Prandtl
      real*8  alphaT(nzWT)
C
      logical DoIonDrag
      logical DimLessStVct, DoNormEigVctZGEEV
C
      real*8 tnW0(nzW)
      real*8 roW0(nzW)
      real*8 uT0(nzW)
      real*8 CvW0(nzW)
      real*8 gW0(nzW)
      real*8 HDensW0(nzW)
      real*8 GammaW0(nzW)
      real*8 CsoundW0(nzW)
c
      real*8 DroW0(nzW)
      real*8 DtnW0(nzW)
      real*8 DuT0(nzW)
      real*8 DHDensW0(nzW)
c
      real*8 SinDip, CosDip
      real*8 nu0niW(nzW)
      real*8 deniW0(nzW)
      real*8 nu0inW(nzW)
      real*8 DiffCoef0W(nzW)
      real*8 DiffVelocityW0(nzW)
      real*8 DdeniW0(nzW)
      real*8 DDiffVelocityW0(nzW)
c
c     --- output ---
      real*8  OmegaIm, dmin, dmax
      logical fail
C
C     --- local variables ---
      integer    i, j, k, M, iOmega
      real*8     AbsOmegaIm, Omega, alpha, g0, Cv, kX, kX2
      real*8     mu0, Dmu0, muK0, ro0, DRo0, u0, Du0, T0, DT0
      real*8     HDens, Gamma, Csound2, HDH
      real*8     nu0ni, nu0in, DiffCoef, uD0, deni, Ddeni, DuD0
      real*8     Dmu0ro0, DT0T
      real*8     scale
      real*8     d, eps
      complex*16 im, one, zero
      complex*16 OmegaC, Omega0C, OmegaDop
      complex*16 FX(6), FZ(6), P(6)
      logical    more, cross
C
      integer, allocatable:: MapJ(:)
      complex*16, allocatable:: KRCN(:), EValN(:), EVctN(:,:)
C
      complex*16, allocatable:: Amat(:,:)
      real*8, allocatable:: Wreal( : )
      real*8, allocatable:: x(:), F1(:), F4(:)
C
C     --- allocate ---
      M = 2*NModes
      allocate( Amat(M,M) )
      allocate( KRCN(M), EValN(M), EVctN(M,M) )
      allocate( MapJ(M) )
      allocate( x(NFFT), F1(NFFT), F4(NFFT) )
c
      kX  = 1.d-3 * kT  !kx in 1/m
      kX2 = kX * kX
C
      Omega0C = dcmplx( Omega0, 0.d0 )    !Omega0
C
      im   = ( 0.d0, 1.d0 )
      one  = ( 1.d0, 0.d0 )
      zero = ( 0.d0, 0.d0 )
C     -------------------------------------------------------------------------
C                             Loop over OmegaIm
C     -------------------------------------------------------------------------
      more = .true.
      AbsOmegaIm = AbsOmegaImMin
      do while ( more )
C       -------------------------------------------------------------------------
C                             Loop over frequency
C       -------------------------------------------------------------------------
        do iOmega = 1, NFFT
          OmegaIm = - AbsOmegaIm
c
          Omega   = OmegaK(iOmega)
          OmegaC  = dcmplx( Omega,  OmegaIm )
c
          i = izT(iT) + IdzT / 2
c
          alpha = alphaT(iT)
          T0  = tnW0(i)  !K
          u0  = uT0(i)   !m/s
          ro0 = roW0(i)  !kg/m**3
          Cv  = CvW0(i)  !J/(kg*K)
          g0  = gW0(i)   !m/s**2
          HDens   = HDensW0(i) !m
          Gamma   = GammaW0(i)
          Csound2 = CsoundW0(i)**2 !(m/s)**2
          Dro0 = DroW0(i)  ![ro]/m
          Du0  = DuT0(i)   ![v]/m
          DT0  = DtnW0(i)  !K/m
          DT0T = DT0 / T0
          HDH  = DHDensW0(i) / HDens
c
          nu0ni = nu0niW(i)
          nu0in = nu0inW(i)
          DiffCoef = DiffCoef0W(i)
          uD0   = DiffVelocityW0(i)
          deni  = deniW0(i)
          Ddeni = DdeniW0(i)
          DuD0  = DDiffVelocityW0(i)
C
C         --- dynamic viscosity mu0  ---
          mu0 = 3.34d-7 * tnW0(i)**0.71
C
C         --- kinematic viscosity  ---
          muK0 = mu0 / ro0
C
C         --- derivatives of viscosity Dmu0  ---
          Dmu0    = 0.71d0 * mu0 * DtnW0(i) / tnW0(i)
          Dmu0ro0 = Dmu0 / ro0
c
          if ( TypeLinModel .eq. 2 ) then
            u0      = 0.d0
            Du0     = 0.d0
            DT0T    = 0.d0
            HDH     = 0.d0 !redundant
            Dmu0ro0 = - muK0 / HDens
          end if
C
C         --- Complex Dopler frequency OmegaDop ---
          OmegaDop = OmegaC - dcmplx(kX * u0,0.D0)
C         -----------------------------------------------------------------------
c                      Ion drag terms for b = (-CosDip, 0.0, -SinDip )
c         -----------------------------------------------------------------------
          if ( DoIonDrag ) then
            call IonDragForcesAndHeating
     i         ( DimLessStVct, kX, u0, Du0, T0, DT0T, CosDip, SinDip,
     i           g0, Cv, nu0in, DiffCoef, uD0, HDens, deni, Ddeni, DuD0,
     i           OmegaC, Omega0C, OmegaDop,
     o           FX, FZ, P )
          else
            do k = 1, 6
              FX(k) = zero
              FZ(k) = zero
              P(k)  = zero
            end do
          end if
C       -----------------------------------------------------------------------
C                                  System matrix
c       -----------------------------------------------------------------------
          call SYSTEMMatrix
     i       ( M, DimLessStVct, DoNormEigVctZGEEV,
     i         kX, kX2, Prandtl,
     i         T0, ro0, Cv,
     i         HDens, Gamma, Csound2,
     i         Du0, DT0T, HDH,
     i         mu0, muK0, Dmu0ro0,
     i         Omega0C, OmegaDop,
     i         nu0ni, FX, FZ, P,
     o         EvalN, EVctN, Amat )
c
c         --- normalize eigenvectors by EVctN(2,j) corresponding to qNEV ---
          if ( .not. DoNormEigVctZGEEV ) then
            do j = 1, M
              scale = 1.d0 / abs( EVctN(qNEV,j) )
              do i = 1, M
                EVctN(i,j) = scale * EVctN(i,j)
              end do
            end do
          end if
C
C         --- wavenumbers in 1/km ---
          do j = 1, 6
            KRCN(j) = im * ( EValN(j) - 0.5d0*alpha ) * kT
          end do
C         -----------------------------------------------------------------------
C                               Order solutions
C         ------------------------------------------------------------------------
          allocate( Wreal(M) )
          do i = 1, M
            Wreal(i) = dble( EvalN(i) )
          end do
          call SortEigenValEigenVct( M, M, Wreal, KRCN, EValN, EVctN )
c
          MapJ(1) = 3
          MapJ(2) = 2
          MapJ(3) = 1

          MapJ(4) = 4
          MapJ(5) = 5
          MapJ(6) = 6
C
          x(iOmega)  = dble(OmegaDop)
          F1(iOmega) = dble(EValN(MapJ(1)))
          F4(iOmega) = dble(EValN(MapJ(4)))
C
C         --- deallocate ---
          deallocate( Wreal )
        end do !end loop frequency
c
c       --- sort with respect to intrinsic frequency ---
c       call SortOneArrayPlusTwoArrays( NFFT, NFFT, x, F1, F4  )
C       -----------------------------------------------------------------------
C                                 Crossing curves
C       -----------------------------------------------------------------------
        eps   = 1.0d-10
        cross = .false.
c
        dmin =  1.0d300
        dmax = -1.0d300
        do i = 1, NFFT
          d    = F4(i) - F1(i)
          dmin = min(dmin, d)
          dmax = max(dmax, d)
        end do
        if (dmin .lt. -eps .and. dmax .gt. eps) cross = .true.
        if (abs(dmin) .le. eps .or. abs(dmax) .le. eps) cross = .true.
C
        if ( cross ) then
          if (AbsOmegaIm .lt. AbsOmegaImMax - 2.d0*DOmegaIm ) then
            AbsOmegaIm = AbsOmegaIm + DOmegaIm
          else
            more =.false.
            fail = .true.
          end if
        else
           more = .false.
        end if
      end do !end WHILE
C
C     --- deallocate ---
      deallocate( Amat )
      deallocate( KRCN, EValN, EVctN )
      deallocate( MapJ )
      deallocate( x, F1, F4 )
      end
C     *************************************************************************
      subroutine GlobalCausalityCondition
     i         ( NModes, nzW, nzWT, IdzT, izT,
     i           Omega, AbsOmegaImMin, AbsOmegaImMax, DOmegaIm,
     i           Omega0, kT, alphaT, Prandtl,
     i           TypeLinModel, DoIonDrag,
     i           DimLessStVct, DoNormEigVctZGEEV, qNEV,
C
     i           tnW0, roW0, uT0,
     i           CvW0, gW0, HDensW0, GammaW0, CsoundW0,
     i           DroW0, DtnW0, DuT0, DHDensW0,
     i           SinDip, CosDip, nu0niW, deniW0,
     i           nu0inW, DiffCoef0W, DiffVelocityW0,
     i           DdeniW0, DDiffVelocityW0,
C
     o           fail, OmegaIm )
C
C     --- input variables ---
      implicit none
      integer NModes, nzW, nzWT, IdzT, TypeLinModel, qNEV
      integer izT(nzWT+1)
      real*8  Omega, AbsOmegaImMin, AbsOmegaImMax, DOmegaIm
      real*8  Omega0, kT, Prandtl
      real*8  alphaT(nzWT)
C
      logical DoIonDrag
      logical DimLessStVct, DoNormEigVctZGEEV
C
      real*8 tnW0(nzW)
      real*8 roW0(nzW)
      real*8 uT0(nzW)
      real*8 CvW0(nzW)
      real*8 gW0(nzW)
      real*8 HDensW0(nzW)
      real*8 GammaW0(nzW)
      real*8 CsoundW0(nzW)
c
      real*8 DroW0(nzW)
      real*8 DtnW0(nzW)
      real*8 DuT0(nzW)
      real*8 DHDensW0(nzW)
c
      real*8 SinDip, CosDip
      real*8 nu0niW(nzW)
      real*8 deniW0(nzW)
      real*8 nu0inW(nzW)
      real*8 DiffCoef0W(nzW)
      real*8 DiffVelocityW0(nzW)
      real*8 DdeniW0(nzW)
      real*8 DDiffVelocityW0(nzW)
c
c     --- output ---
      real*8  OmegaIm
      logical fail
C
C     --- local variables ---
      integer    i, j, k, iT, M
      real*8     AbsOmegaIm, alpha, g0, Cv, kX, kX2
      real*8     mu0, Dmu0, muK0, ro0, DRo0, u0, Du0, T0, DT0
      real*8     HDens, Gamma, Csound2, HDH
      real*8     nu0ni, nu0in, DiffCoef, uD0, deni, Ddeni, DuD0
      real*8     Dmu0ro0, DT0T
      real*8     scale
      real*8     MaxEVal1, MinEval4
      complex*16 im, one, zero
      complex*16 OmegaC, Omega0C, OmegaDop
      complex*16 FX(6), FZ(6), P(6)
      logical    more, PrintInfo
C
      integer, allocatable:: MapJ(:)
      complex*16, allocatable:: KRCN(:), EValN(:), EVctN(:,:)
C
      complex*16, allocatable:: Amat(:,:)
      real*8, allocatable:: Wreal( : )
      real*8, allocatable:: F1(:), F4(:)
C
      PrintInfo = .false.
C
C     --- allocate ---
      M = 2*NModes
      allocate( Amat(M,M) )
      allocate( KRCN(M), EValN(M), EVctN(M,M) )
      allocate( MapJ(M) )
      allocate( F1(nzwT), F4(nzwT) )
c
      kX  = 1.d-3 * kT  !kx in 1/m
      kX2 = kX * kX
C
      Omega0C = dcmplx( Omega0, 0.d0 )    !Omega0
C
      im   = ( 0.d0, 1.d0 )
      one  = ( 1.d0, 0.d0 )
      zero = ( 0.d0, 0.d0 )
C     -------------------------------------------------------------------------
C                             Loop over OmegaIm
C     -------------------------------------------------------------------------
      more = .true.
      AbsOmegaIm = AbsOmegaImMin
      do while ( more )
        OmegaIm = - AbsOmegaIm
        OmegaC  = dcmplx( Omega,  OmegaIm ) !Omega
C       -------------------------------------------------------------------------
C                             Loop over layers
C       -------------------------------------------------------------------------
        do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
c
          alpha = alphaT(iT)
          T0  = tnW0(i)  !K
          u0  = uT0(i)   !m/s
          ro0 = roW0(i)  !kg/m**3
          Cv  = CvW0(i)  !J/(kg*K)
          g0  = gW0(i)   !m/s**2
          HDens   = HDensW0(i) !m
          Gamma   = GammaW0(i)
          Csound2 = CsoundW0(i)**2 !(m/s)**2
          Dro0 = DroW0(i)  ![ro]/m
          Du0  = DuT0(i)   ![v]/m
          DT0  = DtnW0(i)  !K/m
          DT0T = DT0 / T0
          HDH  = DHDensW0(i) / HDens
c
          nu0ni = nu0niW(i)
          nu0in = nu0inW(i)
          DiffCoef = DiffCoef0W(i)
          uD0   = DiffVelocityW0(i)
          deni  = deniW0(i)
          Ddeni = DdeniW0(i)
          DuD0  = DDiffVelocityW0(i)
C
C         --- dynamic viscosity mu0  ---
          mu0 = 3.34d-7 * tnW0(i)**0.71
C
C         --- kinematic viscosity  ---
          muK0 = mu0 / ro0
C
C         --- derivatives of viscosity Dmu0  ---
          Dmu0    = 0.71d0 * mu0 * DtnW0(i) / tnW0(i)
          Dmu0ro0 = Dmu0 / ro0
c
          if ( TypeLinModel .eq. 2 ) then
            u0      = 0.d0
            Du0     = 0.d0
            DT0T    = 0.d0
            HDH     = 0.d0 !redundant
            Dmu0ro0 = - muK0 / HDens
          end if
C
C         --- Complex Dopler frequency OmegaDop ---
          OmegaDop = OmegaC - dcmplx(kX * u0,0.D0)
C         -----------------------------------------------------------------------
c                        Ion drag terms for b = (-CosDip, 0.0, -SinDip )
c         -----------------------------------------------------------------------
          if ( DoIonDrag ) then
            call IonDragForcesAndHeating
     i         ( DimLessStVct, kX, u0, Du0, T0, DT0T, CosDip, SinDip,
     i           g0, Cv, nu0in, DiffCoef, uD0, HDens, deni, Ddeni, DuD0,
     i           OmegaC, Omega0C, OmegaDop,
     o           FX, FZ, P )
          else
            do k = 1, 6
              FX(k) = zero
              FZ(k) = zero
              P(k)  = zero
            end do
          end if
C       -----------------------------------------------------------------------
C                                 System matrix
c       -----------------------------------------------------------------------
          call SYSTEMMatrix
     i       ( M, DimLessStVct, DoNormEigVctZGEEV,
     i         kX, kX2, Prandtl,
     i         T0, ro0, Cv,
     i         HDens, Gamma, Csound2,
     i         Du0, DT0T, HDH,
     i         mu0, muK0, Dmu0ro0,
     i         Omega0C, OmegaDop,
     i         nu0ni, FX, FZ, P,
     o         EvalN, EVctN, Amat )
c
c         --- normalize eigenvectors by EVctN(2,j) corresponding to qNEV ---
          if ( .not. DoNormEigVctZGEEV ) then
            do j = 1, M
              scale = 1.d0 / abs( EVctN(qNEV,j) )
              do i = 1, M
                EVctN(i,j) = scale * EVctN(i,j)
              end do
            end do
          end if
C
C         --- wavenumbers in 1/km ---
          do j = 1, 6
            KRCN(j) = im * ( EValN(j) - 0.5d0*alpha ) * kT
          end do
C         -----------------------------------------------------------------------
C                                 Order solutions
C         ------------------------------------------------------------------------
          allocate( Wreal(M) )
          do i = 1, M
            Wreal(i) = dble( EvalN(i) )
          end do
          call SortEigenValEigenVct( M, M, Wreal, KRCN, EValN, EVctN )
c
          MapJ(1) = 3
          MapJ(2) = 2
          MapJ(3) = 1

          MapJ(4) = 4
          MapJ(5) = 5
          MapJ(6) = 6
c
          F1(iT) = dble(EValN(MapJ(1)))
          F4(iT) = dble(EValN(MapJ(4)))
C
C         --- deallocate ---
          deallocate( Wreal )
        end do !end layers
C       -------------------------------------------------------------------------
C           Max REALPART(EValT(1,iT)) and Min REALPART(EValT(4,iT))
C       -------------------------------------------------------------------------
        MaxEVal1 = -1.d+20
        MinEVal4 =  1.d+20
        do iT = 1, nzWT
           MaxEVal1 = max( MaxEVal1, F1(iT) )
           MinEVal4 = min( MinEVal4, F4(iT) )
        end do
c
        if ( PrintInfo ) then
          write(6,'(a,1pe11.4)')
     &   'Imaginary frequency shift = ', OmegaIm
          write(6,'(2(a,1pe11.4))') 'MaxEVal1 = ', MaxEVal1,
     &          '; MinEVal4 = ', MinEVal4
        end if
        if ( MaxEVal1 .lt. MinEVal4 ) then
          more = .false.
        else
          if (AbsOmegaIm .lt. AbsOmegaImMax - 2.d0*DOmegaIm ) then
            AbsOmegaIm = AbsOmegaIm + DOmegaIm
          else
            more =.false.
            fail = .true.
          end if
        end if
      end do !end WHILE
C
C     --- deallocate ---
      deallocate( Amat )
      deallocate( KRCN, EValN, EVctN )
      deallocate( MapJ )
      deallocate( F1, F4 )
      end
C     *************************************************************************
C
C
C
C
C                             TypeSolMet = 1
C                    Global-Matrix Method for Amplitudes (GMMA)
C
C
C
C
C     *************************************************************************
      subroutine WaveAmplitudesModelGMMA
     i         ( TypeSclLayEq, TypeBC, qBC, NModes, nzWT, kT, zT,
     i           EValT, EVctT,
     o           AmpT )
c ------------------------------------------------------------------------------
c
c             o zT(nzW+1)
c             |
c             |   layer iT = nzW, AmpT(*,nzW), EVal(*,nzW)
c             |
c BC = nzW-1  o zT(nzW):: ALay1(nzW-1)*AmpT(nzW)-ALay0(nzW-1)*AmpT(nzW-1) = 0
c             |
c             |  layer iT = nzW-1, AmpT(*,nzW-1), EVal(*,nzW-1)
c             |
c             o zT(nzW-1)
c
c            ...
c
c             o zT(iT+2)
c             |
c             |   layer iT + 1, AmpT(*,iT+1), EVal(*,iT+1)
c             |
c    BC = iT  o zT(iT+1): ALay1(iT)*AmpT(iT+1) - ALay0(iT)*AmpT(iT) = 0
c             |
c             |   layer iT, AmpT(*,iT), EVal(6,iT); dz = zT(iT+1)-zT(iT)
c             |
c             o zT(iT)
c
c            ...
c
c             o zT(3)
c             |
c             |   layer iT = 2, AmpT(*,2), EVal(*,2)
c             |
c      BC = 1 o zT(2):ALay1(1)*AmpT(2) - ALay0(1)*AmpT(1) = 0
c             |
c             |   layer iT = 1, AmpT(*,1), EVal(*,1)
c             |
c             o zT(1) = zmin
C
c   Unknowns
c          AT(1),...,AT(nzW) = AT(Nlevels)
c   Layers equations
c          iT = 1, ...,nzW-1
c
c      Calling routines:
c       - LayerMatricesModelGMMA
c       - ExtendLayerMatricesGMM
c       - BoundaryConditionComponent
c       - ZGBTRF, ZGBTRS
c ------------------------------------------------------------------------------
      implicit none
      integer    TypeSclLayEq, TypeBC, qBC, NModes, nzWT
      real*8     KT
      real*8     zT(nzWT+1)
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
!
      complex*16  AmpT(2*NModes,nzWT)
C
C       --- local variables ---
      integer      Nlevels, Nlayers, NMAX,
     &             NSUBDIAG, NSUPDIAG, BAND, LDAARB, INFO
      integer      layer, iT, I, COL, ROWMIN, ROWMAX, ROW,
     &             BANDROW, NL, NC, level, j
      real*8       dz
      complex*16   one, zero
      integer, allocatable:: INDX( : )
      complex*16, allocatable:: ARB( : , : ), BR( : , : )
!
      complex*16, allocatable:: ALayer1(:,:)
      complex*16, allocatable:: ALayer0(:,:)
!
      complex*16, allocatable:: b( : )
!
      real*8     res
      complex*16 resI
      logical    PrintInfo,  CheckBoundaryConditions
!
      PrintInfo = .false.
      CheckBoundaryConditions  = .false.
!
      one  = (1.d0,0.d0)
      zero = (0.d0,0.d0)
!
      Nlevels = nzWT
      Nlayers = nzWT - 1
      NMAX = 2*NModes * Nlevels
C
C     --- dimensions of band matrix ---
      NSUBDIAG = 3 * NModes - 1
      NSUPDIAG = 3 * NModes - 1
      BAND     = NSUBDIAG + NSUPDIAG + 1
C
C     --- allocatable dimension of band matrix ---
      LDAARB = 9 * NModes - 2
C
C     --- allocate band matrix ---
      allocate( ARB ( LDAARB, NMAX ),
     &          INDX( NMAX ),
     &          BR  ( NMAX, 1 ) )
      ARB = zero
      BR  = zero
c     -------------------------------------------------------------------------
C       Insert the boundary conditions at the TOA into the band matrix ARB
c     -------------------------------------------------------------------------
      do I = 1, NModes
        COL    = I + NModes
        ROWMIN = MAX0( 1,    COL - NSUPDIAG )
        ROWMAX = MIN0( NMAX, COL + NSUBDIAG )
        ROW    = I
        if( ( ROWMIN .LE. ROW ) .AND. ( ROW .LE. ROWMAX ) ) then
          BANDROW = ROW - COL + BAND
          ARB(BANDROW,COL) = one
        end if
      end do
c     -------------------------------------------------------------------------
C              Insert the layer matrices into the band matrix ARB
c     -------------------------------------------------------------------------
      allocate( ALayer1(2*NModes,2*NModes) )
      allocate( ALayer0(2*NModes,2*NModes) )
      do layer = 1, Nlayers
        iT = Nlayers - layer + 1  !iT = Nlayers,...,1
        dz = zT(iT+1) - zT(iT)
!       -----------------------------------------------------------------------
!                           Layer matrices-
!       -----------------------------------------------------------------------
        call LayerMatricesModelGMMA
     i     ( iT, TypeSclLayEq, NModes, nzWT, dz, kT,
     i       EValT, EVctT,
     o       ALayer1, ALayer0 )
!      ------------------------------------------------------------------------
!      Extend layer matrices
!            ALayer1(2*NModes,2*NModes) and ALayer0(2*NModes,2*NModes)
!      into ARB starting from
!      (NL,NC) = (NModes,0)...
!                (NModes+(Nlayers-1)*(2*NModes); (Nlayers-1)*(2*NModes)
!      ------------------------------------------------------------------------
        NL = NModes + ( layer - 1 ) * 2*NModes
        NC = ( layer - 1 ) * 2*NModes
        call ExtendLayerMatricesGMM
     I     ( NModes, NMAX, LDAARB, NSUPDIAG, NSUBDIAG,
     I       BAND, NL, NC, ALayer1, ALayer0,
     O       ARB )
      end do
c     -------------------------------------------------------------------------
C       Insert the boundary conditions at the BOTTOM into the band
C       matrix ARB and the vector BR
c     -------------------------------------------------------------------------
      NL = NModes + Nlayers * 2*NModes
      NC = Nlayers * 2*NModes
      do I = 1, NModes
        COL    = I + NC
        ROWMIN = MAX0( 1,    COL - NSUPDIAG )
        ROWMAX = MIN0( NMAX, COL + NSUBDIAG )
        ROW    = I + NL
        if( ( ROWMIN .LE. ROW ) .AND. ( ROW .LE. ROWMAX ) ) then
          BANDROW = ROW - COL + BAND
          ARB(BANDROW,COL) = ONE
        end if
      end do
!
!     --- BR ---
      allocate( b( NModes ) )
      if ( TypeBC .eq. 1 ) then
        do I = 1, NModes
          if ( I .EQ. 1 ) then
            BR(I+NL,1) = one
          else
            BR(I+NL,1) = zero
          end if
        end do
      else
        call BoundaryConditionComponent
     i     ( qBC, NModes, nzWT, kT, EValT, EVctT,
     o       b )
        do I = 1, NModes
          BR(I+NL,1) = b(I)
        end do
      end if
C     -------------------------------------------------------------------------
C                     Solve the system of equations
C     -------------------------------------------------------------------------
      call ZGBTRF( NMAX, NMAX, NSUBDIAG, NSUPDIAG, ARB,
     &     LDAARB, INDX, INFO )
      if ( INFO .NE. 0 ) then
        write(6,*) 'INFO = ', INFO
        STOP 'Error in ZGBTRF'
      end if
C
      call ZGBTRS( 'NO TRANSPOSE', NMAX, NSUBDIAG, NSUPDIAG,
     &     1, ARB, LDAARB, INDX, BR, NMAX, INFO )
      if ( INFO .NE. 0 ) STOP 'Error in ZGBTRS'
C     -------------------------------------------------------------------------
C                        Get amplitudes
C     -------------------------------------------------------------------------
      do level = 1, Nlevels
        iT = Nlevels - level + 1  !iT = Nlevels,...,1
        NL = (level - 1) * 2*NModes
        do i = 1, 2*NModes
          AmpT(i,iT) = BR(NL+i,1)
        end do
      end do
C     -------------------------------------------------------------------------
C               Write amplitudes at top and bottom boundaries
C     -------------------------------------------------------------------------
      if ( PrintInfo ) then
        write(6,'(a)') ' Amplitudes at bottom boundary'
        do i = 1, 2*NModes
          write(6,'(a,i2,a,2(1pe11.4,x))')
     &  ' Mode: ', i, ', Amplitude = ',  AmpT(i,1)
        end do
        write(6,'(a)') ' Amplitudes at top boundary'
        do i = 1, 2*NModes
          write(6,'(a,i2,a,2(1pe11.4,x))')
     &  ' Mode: ', i, ', Amplitude = ',  AmpT(i,Nlevels)
        end do
      end if
c     -------------------------------------------------------------------------
c                     Check boundary condition
c              ALay1(iT)*AmpT(iT+1) - Alay0(iT)*AmpT(iT) = 0
c     -------------------------------------------------------------------------
      if ( PrintInfo .and. CheckBoundaryConditions ) then
        do iT = 1, Nlayers
!         ---------------------------------------------------------------------
!                               Layer matrices
!         ---------------------------------------------------------------------
          call LayerMatricesModelGMMA
     i       ( iT, TypeSclLayEq, NModes, nzWT, dz, kT,
     i         EValT, EVctT,
     o         ALayer1, ALayer0 )
!         ---------------------------------------------------------------------
!                                   Residual
!         ---------------------------------------------------------------------
          res = 0.d0
          do i = 1, 2*NModes
            resI = zero
            do j = 1, 2*NModes
              resI = resI + ALayer1(i,j)*AmpT(j,iT+1)
     &             - ALayer0(i,j)*AmpT(j,iT)
            end do
            res = res + abs(resI)**2
          end do
          write(6,'(a,i6,a,1pe11.4)') 'Layer = ', iT,
     &          ', Residual = ', res
          if ( mod(iT,50).eq.0 ) then
            write(6,'(a)') 'Press Enter to continue'
            read(5,*)
          end if
        end do
      end if
!
      deallocate( ARB, INDX, BR )
      deallocate( ALayer1 )
      deallocate( ALayer0 )
      deallocate( b )
      end
C     ***********************************************************************
      subroutine ExtendLayerMatricesGMM
     I         ( NModes, NMAX, LDAARB, NSUPDIAG, NSUBDIAG,
     I           BAND, NL, NC, ALayer1, ALayer0,
     O           ARB )
C     -----------------------------------------------------------------------
C          Extend the layer matrices
C                   ALayer1, ALayer0
C          into the band matrix. Use the map:
C                   AB(i-j+band,j) = A(i,j)
C          for that rows i satisfying
C                max(1,j - nsupdiag) <= i <= min(n,j + nsubdiag)
C          with
C                band = nsubdiag + nsupdiag + 1
C     -----------------------------------------------------------------------
      implicit none
C
C     --- module arguments ---
      integer          NModes, NMAX, LDAARB, NSUPDIAG, NSUBDIAG,
     &                 BAND, NL, NC
!
      complex*16 ALayer1(2*NModes,2*NModes)
      complex*16 ALayer0(2*NModes,2*NModes)
C
      complex*16 ARB( LDAARB, NMAX )
C
C     --- local variables ---
      integer          J, COL, ROWMIN, ROWMAX, I, ROW, BANDROW
C
      do J = 1, 2*NModes
        COL   = J + NC
        ROWMIN = MAX0( 1,    COL - NSUPDIAG )
        ROWMAX = MIN0( NMAX, COL + NSUBDIAG )
        do I = 1, 2*NModes
          ROW = I + NL
          if ( ( ROWMIN .LE. ROW ) .AND. ( ROW .LE. ROWMAX ) )then
            BANDROW = ROW - COL + BAND
            ARB(BANDROW,COL) = ALayer1(i,j)
          end if
        end do
      end do
      NC = NC + 2*NModes
      do J = 1, 2*NModes
        COL    = J + NC
        ROWMIN = MAX0( 1,    COL - NSUPDIAG )
        ROWMAX = MIN0( NMAX, COL + NSUBDIAG )
        do I = 1, 2*NModes
          ROW = I + NL
          if ( ( ROWMIN .LE. ROW ) .AND. ( ROW .LE. ROWMAX ) )then
            BANDROW = ROW - COL + BAND
            ARB(BANDROW,COL) = - ALayer0(i,j)
          end if
        end do
      end do
      end
!     *************************************************************************
      subroutine BoundaryConditionComponent
     i         ( q, NModes, nzWT, kT, EValT, EVctT,
     o           b )
C     -------------------------------------------------------------------------
c      Calling routines:
C        - ZGETRF, ZGETRI
C     -------------------------------------------------------------------------
      implicit none
      integer    q, NModes, nzWT
      real*8     kT
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
!
      complex*16 b(NModes)
!
!     --- local Variables ---
      integer     m, k
      complex*16, allocatable:: A(:,:)
!
      integer  info
      integer, allocatable:: ipiv(:)
      complex*16, allocatable:: work(:)

      allocate( A(NModes,NModes) )
      do m = 1, NModes
        do k = 1, NModes
          A(m,k) = ( kT*EValT(m,1) )**(k-1) * EVctT(q,m,1)
        end do
      end do
!
      allocate( work( NModes ), ipiv( NModes ) )
      call ZGETRF( NModes, NModes, A, NModes, ipiv, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU factorization routine DGETRF'
        STOP
      end if
      call ZGETRI( NModes, A, NModes, ipiv, work, NModes, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU substitution routine DGETRI'
        STOP
      end if
!
      do m = 1, NModes
        b(m) = A(m,1)
      end do
!
      deallocate( A, work, ipiv )
      end
C     *************************************************************************
      subroutine LayerMatricesModelGMMA
     i         ( iT, TypeSclLayEq, NModes, nzWT, dz, kT,
     i           EValT, EVctT,
     o           ALayer1, ALayer0 )
C     -------------------------------------------------------------------------
c      Calling routines:
C        - ZGETRF, ZGETRI, ZGEMM
C     -------------------------------------------------------------------------
      implicit none
      integer    iT, TypeSclLayEq, NModes, nzWT
      real*8     dz, KT
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
!
      complex*16 ALayer1(2*NModes,2*NModes)
      complex*16 ALayer0(2*NModes,2*NModes)
!
!     --- local variables ---
      integer  Mrank, i, j
      real*8   lambdaR
      complex*16 one, zero
      complex*16, allocatable:: S0(:), S1(:)
      complex*16, allocatable:: VT(:,:), VT1(:,:), V1V(:,:)
      complex*16, allocatable:: VTtemp(:,:)
!
      integer  info
      integer, allocatable:: ipiv(:)
      complex*16, allocatable:: work(:)
!
      Mrank = 2*NModes
!
      one  = (1.d0,0.d0)
      zero = (0.d0,0.d0)
!
      ALayer1 = zero
      ALayer0 = zero
!     -------------------------------------------------------------------------
!                          Scaling matrices
!     -------------------------------------------------------------------------
      allocate( S0(Mrank), S1(Mrank) )
      if ( TypeSclLayEq .eq. 1 ) then
        do j = 1, Mrank
          lambdaR = dble( EValT(j,iT) )
          if ( lambdaR .gt. 0.d0 ) then
            S0(j) = one
            S1(j) = exp( -kT * dz * EValT(j,iT) )
          else
            S0(j) = exp( kT * dz * EValT(j,iT) )
            S1(j) = one
          end if
        end do
      else
        do j = 1, Mrank
          if ( j .le. NModes ) then
            S1(j) = one
          else
            S1(j) = exp( -kT * dz * EValT(j,iT) )
          end if
!
          if ( j .le. NModes ) then
            S0(j) = exp( kT * dz * EValT(j,iT) )
          else
            S0(j) = one
          end if
        end do
      end if
!     -------------------------------------------------------------------------
!                           Eigenvector matrices
!     -------------------------------------------------------------------------
      allocate( VT(Mrank,Mrank), VT1(Mrank,Mrank) )
      allocate( VTtemp(Mrank,Mrank) )
      do i = 1, Mrank
        do j = 1, Mrank
          VT(i,j)  = EVctT(i,j,iT)
          VT1(i,j) = EVctT(i,j,iT+1)
!
          VTtemp(i,j)  = EVctT(i,j,iT)
        end do
      end do
C     -------------------------------------------------------------------------
c                             VT^(-1)
c     -------------------------------------------------------------------------
      allocate( work( Mrank ), ipiv( Mrank ) )
      call ZGETRF( Mrank, Mrank, VT, Mrank, ipiv, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU factorization routine DGETRF'
        write( 6, * ) 'in layer ', iT, '/', nzWT
        do i = 1, Mrank
          write(6,'(100(1pe11.4))') ( Vttemp(i,j),j = 1,Mrank)
        end do
        STOP
      end if
      call ZGETRI( Mrank, VT, Mrank, ipiv, work, Mrank, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU substitution routine DGETRI'
        STOP
      end if
C     -------------------------------------------------------------------------
c                          VT^(-1)*VT1
c     -------------------------------------------------------------------------
      allocate( V1V(Mrank,Mrank) )
      call ZGEMM('N','N', Mrank, Mrank, Mrank, one, VT, Mrank,
     &     VT1, Mrank, zero, V1V, Mrank)
C     -------------------------------------------------------------------------
c                          Layer matrices
c     -------------------------------------------------------------------------
      do i = 1, Mrank
        do j = 1, Mrank
          ALayer1(i,j) = S1(i) * V1V(i,j)
        end do
        ALayer0(i,i) = S0(i)
      end do
!
      deallocate( S0, S1 )
      deallocate( VT, VT1, VTtemp )
      deallocate( work, ipiv )
      deallocate( V1V )
      end
c     *************************************************************************
      subroutine WaveParametersModelGMMAComplete
     i         ( NModes, nzWT, nzW, IdzT, izT,
     i           Omega, OmegaIm, Omega0, kT,
     i           EValT, EVctT, AmpT,
     i           TypeBC, DimLessStVct, DoScaling,
     i           t0BC, u0BC, w0BC, qBC,
     i           tnW0, roW0, pW0, uT0, HDensW0,
     o           UThat, UThat0, UThatPlus, UThatMinus,
     o           URhat, URhat0, URhatPlus, URhatMinus,
     o           TWhat, TWhat0, TWhatPlus, TWhatMinus,
     o           pWhat, roWhat )
      implicit none
      integer    NModes, nzWT, nzW, IdzT, TypeBC, qBC
      integer    izT(nzWT+1)
      real*8     Omega, OmegaIm, Omega0, KT
      real*8     t0BC, u0BC, w0BC
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
      complex*16 AmpT(2*NModes,nzWT)
      logical    DimLessStVct, DoScaling
c
      real*8 tnW0(nzW)
      real*8 roW0(nzW)
      real*8 pW0(nzW)
      real*8 uT0(nzW)
      real*8 HDensW0(nzW)
c
c     --- outputs ---
      complex*16 UThat(nzWT)
      complex*16 UThat0(nzWT)
      complex*16 UThatPlus(nzWT)
      complex*16 UThatMinus(nzWT)

      complex*16 URhat(nzWT)
      complex*16 URhat0(nzWT)
      complex*16 URhatPlus(nzWT)
      complex*16 URhatMinus(nzWT)

      complex*16 TWhat(nzWT)
      complex*16 TWhat0(nzWT)
      complex*16 TWhatPlus(nzWT)
      complex*16 TWhatMinus(nzWT)

      complex*16 pWhat(nzWT)
      complex*16 roWhat(nzWT)
c
c     --- local variables ---
      integer     Mrank, iT, i, j, k
      real*8      kX, FctDimLessStVct, Src0
      complex*16  im, zero, sum
      complex*16  U, W, T, Wcal
      complex*16  OmegaC, Omega0C, OmegaDop, fct
c
      complex*16, allocatable:: b( : )
      complex*16, allocatable:: VT(:,:)
      complex*16, allocatable:: eT0(:,:)
      complex*16, allocatable:: eTPlus(:,:)
      complex*16, allocatable:: eTMinus(:,:)
c
      complex*16, allocatable:: EVct(:)
C
      Mrank = 2*NModes
c
      kX  = 1.d-3*kT  !kx in 1/m
C
      OmegaC  = dcmplx( Omega,  OmegaIm ) !Omega
      Omega0C = dcmplx( Omega0, 0.d0 )    !Omega0
c
      im   = (0.d0,1.d0)
      zero = (0.d0,0.d0)
c
c     --- factor for obtaining a dimensionless system of equations ---
      if ( DimLessStVct ) then
        FctDimLessStVct = kX
      else
        FctDimLessStVct = 1.d0
      end if
c
c     --- eigenvector matrix ---
      allocate( b( NModes ) )
      allocate( VT(Mrank,Mrank) )
      allocate( eT0(2*NModes,nzWT) )
      allocate( eTPlus(2*NModes,nzWT) )
      allocate( eTMinus(2*NModes,nzWT) )
c
c     --- amplitude vector for BC on state vector component ---
      call BoundaryConditionComponent
     i     ( qBC, NModes, nzWT, kT, EValT, EVctT,
     o       b )
c
c    --- amplitudes at top and bottom from boundary conditions ---
      if ( TypeBC .eq. 1 ) then
        do i = 1, 2*NModes
          eTPlus(i,1)     = EVctT(i,1,1)
          eTMinus(i,nzWT) = zero
        end do
      else
        do i = 1, 2*NModes
          sum = zero
          do j = 1, NModes
            sum = sum + EVctT(i,j,1) * b(j)
          end do
          eTPlus(i,1) = sum
          eTMinus(i,nzWT) = zero
        end do
      end if
c
c     --- remaining amplitudes  ---
      do iT = 1, nzwT
        do i = 1, Mrank
          do j = 1, Mrank
            VT(i,j)  = EVctT(i,j,iT)
          end do
        end do
c
        if ( iT .ge. 2 ) then
          do i = 1,  Mrank
            sum = zero
            do j = 1, NModes
              sum = sum + VT(i,j) * AmpT(j,iT)
            end do
            eTPlus(i,iT) = sum
          end do
        end if
        if ( iT .le. nzWT - 1 ) then
          do i = 1,  Mrank
            sum = zero
            do j = 1, NModes
              sum = sum + VT(i,j+NModes) * AmpT(j+NModes,iT)
            end do
            eTMinus(i,iT) = sum
          end do
        end if
      end do
c
c     --- eT0(*,iT) = eTPlus(*,iT)+ eTMinus(*,iT)  ---
      do iT = 1, nzWT
        do i = 1, 2*NModes
          eT0(i,iT) = eTPlus(i,iT) + eTMinus(i,iT)
        end do
      end do
c     -------------------------------------------------------------------------
c                           Wave Amplitudes
c     -------------------------------------------------------------------------
      allocate( EVct( 2*NModes ) )
      do iT = 1, nzWT
c
c       --- midpoint of each isothermal region ---
        i = izT(iT) + IdzT / 2
        OmegaDop = OmegaC - dcmplx(kX * uT0(i),0.D0)
        fct = Omega0C / OmegaDop
c       -----------------------------------------------------------------------
c                e(2*NModes) = V(2*NModes,2*NModes) * a(2*NModes)
c       -----------------------------------------------------------------------
        do j = 1, 2*NModes
          sum = zero
          do k = 1, 2*NModes
            sum = sum + EVctT(j,k,iT) * AmpT(k,iT)
          end do
          EVct(j) = sum
        end do
c
        U  = EVct(1)
        W  = EVct(2)
        T  = EVct(3)
        Wcal = FctDimLessStVct * EVct(5)
c
        UThat(iT)  = Omega0C * U / kX !in m/s
        URhat(iT)  = Omega0C * W / kX !in m/s
        TWhat(iT)  = tnW0(i) * T       !in K
c
        roWhat(iT) = roW0(i) * ( fct * U
     &             - im * fct * W / (kX * HDensW0(i))
     &             + im * fct * Wcal / kX )  !in Kg/m**3
        pWhat(iT) = pW0(i) * ( TWhat(iT)/tnW0(i) + roWhat(iT)/roW0(i) ) !in N/m**2
c       -----------------------------------------------------------------------
c             e0(2*NModes) = ePlus(2*NModes) + eTMinus(2*NModes)
c       -----------------------------------------------------------------------
        do j = 1, 2*NModes
          EVct(j) = eT0(j,iT)
        end do
c
        U  = EVct(1)
        W  = EVct(2)
        T  = EVct(3)
!
        UThat0(iT)  = Omega0C * U / kX !in m/s
        URhat0(iT)  = Omega0C * W / kX !in m/s
        TWhat0(iT)  = tnW0(i) * T      !in K
c       -----------------------------------------------------------------------
c           ePlus(2*NModes) = VPlus(2*NModes,NModes) * aPlus(NModes)
c       -----------------------------------------------------------------------
        do j = 1, 2*NModes
          EVct(j) = eTPlus(j,iT)
        end do
c
        U  = EVct(1)
        W  = EVct(2)
        T  = EVct(3)
c
        UThatPlus(iT)  = Omega0C * U / kX !in m/s
        URhatPlus(iT)  = Omega0C * W / kX !in m/s
        TWhatPlus(iT)  = tnW0(i) * T      !in K
c       -----------------------------------------------------------------------
c         eMinus(2*NModes) = VMinus(2*NModes,NModes) * aMinus(NModes)
c       -----------------------------------------------------------------------
        do j = 1, 2*NModes
          EVct(j) = eTMinus(j,iT)
        end do
c
        U  = EVct(1)
        W  = EVct(2)
        T  = EVct(3)
c
        UThatMinus(iT)  = Omega0C * U / kX !in m/s
        URhatMinus(iT)  = Omega0C * W / kX !in m/s
        TWhatMinus(iT)  = tnW0(i) * T      !in K
      end do
c
c     --- scale quantities ---
      if ( DoScaling ) then
        if ( qBC .eq. 1 ) then
          Src0 = u0BC / abs( dble(UThat(1)) )
        else if ( qBC .eq. 2 ) then
          Src0 = w0BC / abs( dble(URhat(1)) )
        else
          Src0 = t0BC / abs( dble(TWhat(1)) )
        end if
        do iT = 1, nzWT
          UThat(iT)  = Src0 * UThat(iT)
          URhat(iT)  = Src0 * URhat(iT)
          TWhat(iT)  = Src0 * TWhat(iT)
          roWhat(iT) = Src0 * roWhat(iT)
          pWhat(iT)  = Src0 * pWhat(iT)
c
          UThat0(iT)  = Src0 * UThat0(iT)
          URhat0(iT)  = Src0 * URhat0(iT)
          TWhat0(iT)  = Src0 * TWhat0(iT)
c
          UThatPlus(iT)  = Src0 * UThatPlus(iT)
          URhatPlus(iT)  = Src0 * URhatPlus(iT)
          TWhatPlus(iT)  = Src0 * TWhatPlus(iT)
c
          UThatMinus(iT)  = Src0 * UThatMinus(iT)
          URhatMinus(iT)  = Src0 * URhatMinus(iT)
          TWhatMinus(iT)  = Src0 * TWhatMinus(iT)
        end do
      end if
c
      deallocate( b )
      deallocate( VT )
      deallocate( eT0 )
      deallocate( eTPlus )
      deallocate( eTMinus )
      deallocate( EVct )
      end
C     *************************************************************************
C
C
C
C
C                           TypeSolMet = 2
C                Global-Matrix Method for Nodal Values (GMMN)
C
C
C
C
C     *************************************************************************
      subroutine WaveAmplitudesModelGMMN
     i         ( TypeSclLayEq, TypeBC, qBC, NModes, nzWT, kT, zT,
     i           EValT, EVctT,
     o           eT, eT0, eTPlus, eTMinus )
c------------------------------------------------------------------------------
c
c             o zT(nzW+1), eT(nzW+1)
c             |
c             |   layer iT = nzW,
c             |
c             o zT(nzW), eT(nzW)
c             |
c BC = nzW-1  |  layer iT = nzW-1, A1(nzW-1)*eT(nzW) - A0(nzW-1)*eT(nzW-1) = 0
c             |
c             o zT(nzW-1), eT(nzW-1)
c
c            ...
c
c             o zT(iT+1), eT(iT+1)
c             |
c    BC = iT  |   layer iT, A1(iZ)*eT(iT+1) - A0(iT)*eT(iT) = 0
c             |
c             o zT(iT), eT(iT)
c
c            ...
c
c             o zT(3), eT(3)
c             |
c    BC = 2   |   layer iT = 2, A1(t=2)*eT(3) - A0(t=2)*eT(2) = 0
c             |
c             o zT(2), eT(2)
c             |
c    BC = 1   |   layer iT = 1, A1(t=1)*eT(2) - A0(t=1)*eT(1) = 0
c             |
c             o zT(1) = zmin, eT(1)
c
c   Unknowns
c          eT(1),...,eT(nzW) = eT(Nlevels)
c   Layers equations
c          iT = 1, ...,nzW-1
C
C   Calling routines
c     - LayerMatricesModelGMMN
c     - ExtendLayerMatricesGMM
c     - ZGBTRF, ZGBTRS
C     - ETPlusETMinusGMMN
c------------------------------------------------------------------------------
      implicit none
      integer    TypeSclLayEq, TypeBC, qBC, NModes, nzWT
      real*8     KT
      real*8     zT(nzWT+1)
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
!
      complex*16  eT(2*NModes,nzWT)
      complex*16  eT0(2*NModes,nzWT)
      complex*16  eTPlus(2*NModes,nzWT)
      complex*16  eTMinus(2*NModes,nzWT)
C
C       --- local variables ---
      integer      Nlevels, Nlayers, NMAX,
     &             NSUBDIAG, NSUPDIAG, BAND, LDAARB, INFO
      integer      layer, iT, I, COL, ROWMIN, ROWMAX, ROW,
     &             BANDROW, NL, NC, level, j
      real*8       dz
      complex*16   one, zero
      integer, allocatable:: INDX( : )
      complex*16, allocatable:: ARB( : , : ), BR( : , : )
!
      complex*16, allocatable:: VT(:,:)
      complex*16, allocatable:: VTI(:,:)
      complex*16, allocatable:: ALayer1(:,:)
      complex*16, allocatable:: ALayer0(:,:)
      complex*16, allocatable:: ATop(:,:)
      complex*16, allocatable:: ABot(:,:)
c
      complex*16, allocatable:: b( : )
!
      real*8     res
      complex*16 resI
      logical    PrintInfo, CheckBoundaryConditions
!
      PrintInfo = .false.
      CheckBoundaryConditions  = .false.
!
      one  = (1.d0,0.d0)
      zero = (0.d0,0.d0)
!
      Nlevels = nzWT
      Nlayers = nzWT - 1
      NMAX = 2*NModes * Nlevels
C
C     --- dimensions of band matrix ---
      NSUBDIAG = 3 * NModes - 1
      NSUPDIAG = 3 * NModes - 1
      BAND     = NSUBDIAG + NSUPDIAG + 1
C
C     --- allocatable dimension of band matrix ---
      LDAARB = 9 * NModes - 2
C
C     --- allocate band matrix ---
      allocate( ARB ( LDAARB, NMAX ),
     &          INDX( NMAX ),
     &          BR  ( NMAX, 1 ) )
      ARB = zero
      BR  = zero
c     -------------------------------------------------------------------------
C             Insert the layer matrices into the band matrix ARB
c     -------------------------------------------------------------------------
      allocate( VT(2*NModes,2*NModes) )
      allocate( VTI(2*NModes,2*NModes) )
      allocate( ALayer1(2*NModes,2*NModes) )
      allocate( ALayer0(2*NModes,2*NModes) )
      allocate( ATop(NModes,2*NModes) )
      allocate( ABot(NModes,2*NModes) )
!     -------------------------------------------------------------------------
!                    Extend ATop(NModes,2*NModes) into ARB
!     -------------------------------------------------------------------------
      dz = zT(nzWT+1) - zT(nzWT)
      call LayerMatricesModelGMMN
     i   ( nzWT, TypeSclLayEq, NModes, nzWT, dz, kT,
     i     EValT, EVctT,
     o     VT, VTI, ALayer1, ALayer0, ATop, ABot )
      DO J = 1, 2*NModes
        COL    = J
        ROWMIN = MAX0( 1,    COL - NSUPDIAG )
        ROWMAX = MIN0( NMAX, COL + NSUBDIAG )
        DO I = 1, NModes
          ROW = I
          IF ( ( ROWMIN .LE. ROW ) .AND. ( ROW .LE. ROWMAX ) ) THEN
            BANDROW = ROW - COL + BAND
            ARB(BANDROW,COL) = ATop(I,J)
          END IF
        END DO
      END DO
!     -------------------------------------------------------------------------
!                           Layer Equations
!     -------------------------------------------------------------------------
      do layer = 1, Nlayers
        iT = Nlayers - layer + 1  !iT = Nlayers,...,1
        dz = zT(iT+1) - zT(iT)
!       -----------------------------------------------------------------------
!                           layer matrices
!       -----------------------------------------------------------------------
        call LayerMatricesModelGMMN
     i     ( iT, TypeSclLayEq, NModes, nzWT, dz, kT,
     i       EValT, EVctT,
     o       VT, VTI, ALayer1, ALayer0, ATop, ABot )
!      ------------------------------------------------------------------------
!      Extend layer matrices
!            ALayer1(2*NModes,2*NModes) and ALayer0(2*NModes,2*NModes)
!      into ARB starting from
!      (NL,NC) = (NModes,0)...
!                (NModes+(Nlayers-1)*(2*NModes); (Nlayers-1)*(2*NModes)
!      ------------------------------------------------------------------------
        NL = NModes + ( layer - 1 ) * 2*NModes
        NC = ( layer - 1 ) * 2*NModes
        call ExtendLayerMatricesGMM
     I     ( NModes, NMAX, LDAARB, NSUPDIAG, NSUBDIAG,
     I       BAND, NL, NC, ALayer1, ALayer0,
     O       ARB )
!       -----------------------------------------------------------------------
!                 Extend ABot(NModes,2*NModes) into ARB
!       -----------------------------------------------------------------------
        if ( iT .eq. 1 ) then
          NL = NModes + Nlayers * 2*NModes
          NC = Nlayers * 2*NModes
          DO J = 1, 2*NModes
            COL    = J + NC
            ROWMIN = MAX0( 1,    COL - NSUPDIAG )
            ROWMAX = MIN0( NMAX, COL + NSUBDIAG )
            DO I = 1, NModes
              ROW = I + NL
              IF ( ( ROWMIN .LE. ROW ) .AND. ( ROW .LE. ROWMAX ) ) THEN
                BANDROW = ROW - COL + BAND
                ARB(BANDROW,COL) = ABot(I,J)
              END IF
            END DO
          END DO
        end if
!
      end do
c     -------------------------------------------------------------------------
C      Insert the boundary conditions at the BOTTOM into the vector BR
c     -------------------------------------------------------------------------
      NL = NModes + Nlayers * 2*NModes
      allocate( b( NModes ) )
      if ( TypeBC .eq. 1 ) then
        do I = 1, NModes
          if ( I .EQ. 1 ) then
            BR(I+NL,1) = one
          else
            BR(I+NL,1) = zero
          end if
        end do
      else
        call BoundaryConditionComponent
     i     ( qBC, NModes, nzWT, kT, EValT, EVctT,
     o       b )
        do I = 1, NModes
          BR(I+NL,1) = b(I)
        end do
      end if
C     -------------------------------------------------------------------------
C                     Solve the system of equations
C     -------------------------------------------------------------------------
      call ZGBTRF( NMAX, NMAX, NSUBDIAG, NSUPDIAG, ARB,
     &     LDAARB, INDX, INFO )
      if ( INFO .NE. 0 ) then
        write(6,*) 'INFO = ', INFO
        STOP 'Error in ZGBTRF'
      end if
C
      call ZGBTRS( 'NO TRANSPOSE', NMAX, NSUBDIAG, NSUPDIAG,
     &     1, ARB, LDAARB, INDX, BR, NMAX, INFO )
      if ( INFO .NE. 0 ) STOP 'Error in ZGBTRS'
C     -------------------------------------------------------------------------
C                        Get solutions
C     -------------------------------------------------------------------------
      do level = 1, Nlevels
        iT = Nlevels - level + 1  !iT = Nlevels,...,1
        NL = (level - 1) * 2*NModes
        do i = 1, 2*NModes
          eT(i,iT) = BR(NL+i,1)
        end do
      end do
C     -------------------------------------------------------------------------
C               Write amplitudes at top and bottom boundaries
C     -------------------------------------------------------------------------
      if ( PrintInfo ) then
        write(6,'(a)') ' Solutions at bottom boundary'
        do i = 1, 2*NModes
          write(6,'(a,i2,a,2(1pe11.4,x))')
     &  ' Mode: ', i, ', Solution = ',  eT(i,1)
        end do
        write(6,'(a)') ' Solutions at top boundary'
        do i = 1, 2*NModes
          write(6,'(a,i2,a,2(1pe11.4,x))')
     &  ' Mode: ', i, ', Solution = ',  eT(i,Nlevels)
        end do
      end if
C     -------------------------------------------------------------------------
C                        Ascending and Descending Waves
C     -------------------------------------------------------------------------
      call ETPlusETMinusGMMN
     i   ( NModes, nzWT, kT, zT,
     i     EValT, EVctT, eT,
     o     eT0, eTPlus, eTMinus )
c     -------------------------------------------------------------------------
c                     Check boundary condition
c         ALay1(iT)*AmpT(iT+1) - Alay0(iT)*AmpT(iT) = 0
c     -------------------------------------------------------------------------
      if ( PrintInfo .and. CheckBoundaryConditions ) then
        do iT = 1, Nlayers
!         ---------------------------------------------------------------------
!                              Layer matrices
!         ---------------------------------------------------------------------
          call LayerMatricesModelGMMN
     i       ( iT, TypeSclLayEq, NModes, nzWT, dz, kT,
     i         EValT, EVctT,
     o         VT, VTI, ALayer1, ALayer0, ATop, ABot )
!         ---------------------------------------------------------------------
!                                Residual
!         ---------------------------------------------------------------------
          res = 0.d0
          do i = 1, 2*NModes
            resI = zero
            do j = 1, 2*NModes
              resI = resI + ALayer1(i,j)*eT(j,iT+1)
     &             - ALayer0(i,j)*eT(j,iT)
            end do
            res = res + abs(resI)**2
          end do
          write(6,'(a,i6,a,1pe11.4)') 'Layer = ', iT,
     &          ', Residual = ', res
          if ( mod(iT,50).eq.0 ) then
            write(6,'(a)') 'Press Enter to continue'
            read(5,*)
          end if
        end do
      end if
!
      deallocate( ARB, INDX, BR )
      deallocate( VT, VTI )
      deallocate( ALayer1 )
      deallocate( ALayer0 )
      deallocate( ATop, ABot )
      deallocate( b )
      end
C     *************************************************************************
      subroutine LayerMatricesModelGMMN
     i         ( iT, TypeSclLayEq, NModes, nzWT, dz, kT,
     i           EValT, EVctT,
     o           VT, VTI, ALayer1, ALayer0, ATop, ABot )
c     -------------------------------------------------------------------------
c       Calling routines:
c           - ZGETRF, ZGETRI
c     -------------------------------------------------------------------------
      implicit none
      integer    iT, TypeSclLayEq, NModes, nzWT
      real*8     dz, KT
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
!
      complex*16 VT(2*NModes,2*NModes)
      complex*16 VTI(2*NModes,2*NModes)
!
      complex*16 ALayer1(2*NModes,2*NModes)
      complex*16 ALayer0(2*NModes,2*NModes)
      complex*16 ATop(NModes,2*NModes)
      complex*16 ABot(NModes,2*NModes)
!
!     --- local variables ---
      integer  Mrank, i, j
      real*8   lambdaR
      complex*16 one, zero
      complex*16, allocatable:: S0(:), S1(:)
!
      integer  info
      integer, allocatable:: ipiv(:)
      complex*16, allocatable:: work(:)
!
      Mrank = 2*NModes
!
      one  = (1.d0,0.d0)
      zero = (0.d0,0.d0)
!
      ALayer1 = zero
      ALayer0 = zero
      ATop    = zero
      ABot    = zero
!     -------------------------------------------------------------------------
!                              Scaling matrices
!     -------------------------------------------------------------------------
      allocate( S0(Mrank), S1(Mrank) )
      if ( TypeSclLayEq .eq. 1 ) then
        do j = 1, Mrank
          lambdaR = dble( EValT(j,iT) )
          if ( lambdaR .gt. 0.d0 ) then
            S0(j) = one
            S1(j) = exp( -kT * dz * EValT(j,iT) )
          else
            S0(j) = exp( kT * dz * EValT(j,iT) )
            S1(j) = one
          end if
        end do
      else
        do j = 1, Mrank
          if ( j .le. NModes ) then
            S1(j) = one
          else
            S1(j) = exp( -kT * dz * EValT(j,iT) )
          end if
!
          if ( j .le. NModes ) then
            S0(j) = exp( kT * dz * EValT(j,iT) )
          else
            S0(j) = one
          end if
        end do
      end if
!     -------------------------------------------------------------------------
!                          Eigenvector matrices
!     -------------------------------------------------------------------------
      do i = 1, Mrank
        do j = 1, Mrank
          VT(i,j)  = EVctT(i,j,iT)
          VTI(i,j) = VT(i,j)
        end do
      end do
C     -------------------------------------------------------------------------
c                            VTI = VT^(-1)
c     -------------------------------------------------------------------------
      allocate( work( Mrank ), ipiv( Mrank ) )
      call ZGETRF( Mrank, Mrank, VTI, Mrank, ipiv, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU factorization routine DGETRF'
        write( 6, * ) 'in layer ', iT, '/', nzWT
        do i = 1, Mrank
          write(6,'(100(1pe11.4))') ( VT(i,j),j = 1,Mrank)
        end do
        STOP
      end if
      call ZGETRI( Mrank, VTI, Mrank, ipiv, work, Mrank, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU substitution routine DGETRI'
        STOP
      end if
C     -------------------------------------------------------------------------
c                          Layer matrices
c        Alayer1 = S1 * VT^(-1), Alayer0 = S0 * VT^(-1)
c     -------------------------------------------------------------------------
      do i = 1, Mrank
        do j = 1, Mrank
          ALayer1(i,j) = S1(i) * VTI(i,j)
          ALayer0(i,j) = S0(i) * VTI(i,j)
        end do
      end do
C     -------------------------------------------------------------------------
c                       ATop = [0,I] * VT^(-1)
c     -------------------------------------------------------------------------
      if ( iT .eq. nzWT ) then
        do i = 1, NModes
          do j = 1, 2*NModes
            ATop(i,j) = VTI(i+NModes,j)
          end do
        end do
      end if
C     -------------------------------------------------------------------------
c                       ABot = [I,0] * VT^(-1)
c     -------------------------------------------------------------------------
      if ( iT .eq. 1 ) then
        do i = 1, NModes
          do j = 1, 2*NModes
            ABot(i,j) = VTI(i,j)
          end do
        end do
      end if
!
      deallocate( S0, S1 )
      deallocate( work, ipiv )
      end
!     *************************************************************************
      subroutine ETPlusETMinusGMMN
     i         ( NModes, nzWT, kT, zT,
     i           EValT, EVctT, eT,
     o           eT0, eTPlus, eTMinus )
c     -------------------------------------------------------------------------
c       Calling routines:
c           - TPlusTMinusGMMN
c     -------------------------------------------------------------------------
      implicit none
      integer    NModes, nzWT
      real*8     KT
      real*8     zT(nzWT+1)
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
      complex*16  eT(2*NModes,nzWT)
!
      complex*16  eT0(2*NModes,nzWT)
      complex*16  eTPlus(2*NModes,nzWT)
      complex*16  eTMinus(2*NModes,nzWT)
!
!     --- local variables ---
      integer    iT, i, j
      complex*16 one, zero, sum
      complex*16, allocatable:: TPlusT(:,:, :)
      complex*16, allocatable:: TMinusT(:,:, :)
!
      one  = (1.d0,0.d0)
      zero = (0.d0,0.d0)
!
      allocate( TPlusT(2*NModes,2*NModes, nzWT) )
      allocate( TMinusT(2*NModes,2*NModes, nzWT) )
!
      call TPlusTMinusGMMN
     i   ( NModes, nzWT, kT, zT, EValT, EVctT,
     o     TPlusT, TMinusT )
!
!    --- amplitudes at top and bottom from boundary conditions ---
      do i = 1, 2*NModes
        eTPlus(i,1)     = EVctT(i,1,1)
        eTMinus(i,nzWT) = zero
      end do
!
!     --- upward recurrnec ---
      do iT = 1, nzwT - 1
        do i = 1, 2*NModes
          sum = zero
          do j = 1, 2*NModes
            sum = sum + TPlusT(i,j,iT) * eT(j,iT)
          end do
          eTPlus(i,iT+1) = sum
        end do
      end do
!
!     --- downward recurrnec ---
      do iT = nzwT - 1, 1, -1
        do i = 1, 2*NModes
          sum = zero
          do j = 1, 2*NModes
            sum = sum + TMinusT(i,j,iT) * eT(j,iT+1)
          end do
          eTMinus(i,iT) = sum
        end do
      end do
!
!     --- eT (*,iT) = eTPlus(*,iT)+ eTMinus(*,iT)  ---
      do iT = 1, nzWT
        do i = 1, 2*NModes
          eT0(i,iT) = eTPlus(i,iT) + eTMinus(i,iT)
        end do
      end do
!
      deallocate( TPlusT, TMinusT )
      end
!     *************************************************************************
      subroutine TPlusTMinusGMMN
     i         ( NModes, nzWT, kT, zT, EValT, EVctT,
     o           TPlusT, TMinusT )
c     -------------------------------------------------------------------------
c       Calling routines:
c           - ZGETRF, ZGETRI, ZGEMM
c     -------------------------------------------------------------------------
      implicit none
      integer    NModes, nzWT
      real*8     kT
      real*8     zT(nzWT+1)
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
!
      complex*16 TPlusT(2*NModes,2*NModes, nzWT)
      complex*16 TMinusT(2*NModes,2*NModes, nzWT)
!
!     --- local variables ---
      integer     iT, i, j, Mrank
      real*8      dz, arg
      complex*16  one, zero
      complex*16, allocatable:: VT(:,:), VT1(:,:)
      complex*16, allocatable:: lambdaPlus(:), lambdaMinus(:)
      complex*16, allocatable:: x(:,:), TPlus(:,:), TMinus(:,:)
!
      integer  info
      integer, allocatable:: ipiv(:)
      complex*16, allocatable:: work(:)
!
      Mrank = 2*NModes
!
      one  = (1.d0,0.d0)
      zero = (0.d0,0.d0)
!
      allocate( VT(Mrank,Mrank), VT1(Mrank,Mrank) )
      allocate( work( Mrank ), ipiv( Mrank ) )
      allocate( lambdaPlus(Mrank) )
      allocate( lambdaMinus(Mrank) )
!
      allocate( x(Mrank,Mrank) )
      allocate( TPlus(Mrank,Mrank) )
      allocate( TMinus(Mrank,Mrank) )
!
      do iT = 1, nzWT
        dz  = zT(iT+1) - zT(iT)
        arg = kT * dz
!       ------------------------------------------------------------------------
!                             Eigenvector matrices
!       ------------------------------------------------------------------------
        do i = 1, Mrank
          do j = 1, Mrank
            VT(i,j)  = EVctT(i,j,iT)
            VT1(i,j) = EVctT(i,j,iT)
          end do
        end do
!       ------------------------------------------------------------------------
!                              Scaling Matrices
!       ------------------------------------------------------------------------
        lambdaPlus  = zero
        lambdaMinus = zero
        do j = 1, NModes
          lambdaPlus(j)         = exp( arg * EValT(j,iT) )
          lambdaMinus(j+NModes) = exp( - arg * EValT(j+NModes,iT) )
        end do
C       -----------------------------------------------------------------------
c                             VT1 = VT^(-1)
c       -----------------------------------------------------------------------
        call ZGETRF( Mrank, Mrank, VT1, Mrank, ipiv, info)
        if ( info .ne. 0 ) then
          write( 6, * ) 'Error in LU factorization routine DGETRF:VT'
          STOP
        end if
        call ZGETRI( Mrank, VT1, Mrank, ipiv, work, Mrank, info)
        if ( info .ne. 0 ) then
          write( 6, * ) 'Error in LU substitution routine DGETRI'
          STOP
        end if
!       ------------------------------------------------------------------------
!                        x = lambdaPlus * VT^(-1)  for TPlus
!       ------------------------------------------------------------------------
        do i = 1, Mrank
          do j = 1, Mrank
            x(i,j) = lambdaPlus(i) * VT1(i,j)
          end do
        end do
C       ----------------------------------------------------------------------
c                             TPlus = VT * x
c       ----------------------------------------------------------------------
        call ZGEMM('N','N', Mrank, Mrank, Mrank, one, VT, Mrank,
     &       x, Mrank, zero, TPlus, Mrank)
        do i = 1, Mrank
          do j = 1, Mrank
            TPlusT(i,j,iT) = TPlus(i,j)
          end do
        end do
!       ----------------------------------------------------------------------
!                        x = lambdaMinus * VT^(-1) for TMinus
!       ----------------------------------------------------------------------
        do i = 1, Mrank
          do j = 1, Mrank
            x(i,j) = lambdaMinus(i) * VT1(i,j)
          end do
        end do
C       -----------------------------------------------------------------------
c                            TMinus = VT * x
c       -----------------------------------------------------------------------
        call ZGEMM('N','N', Mrank, Mrank, Mrank, one, VT, Mrank,
     &       x, Mrank, zero, TMinus, Mrank)
        do i = 1, Mrank
          do j = 1, Mrank
            TMinusT(i,j,iT) = TMinus(i,j)
          end do
        end do
      end do
!
      deallocate( VT, VT1, lambdaPlus, lambdaMinus, x, TPlus, TMinus )
      deallocate( work, ipiv )
      end
!     *************************************************************************
      subroutine WaveParametersModelGMMN
     i         ( NModes, nzWT, nzW, IdzT, izT,
     i           Omega,  OmegaIm, Omega0, kT,
     i           eT, DimLessStVct, DoScaling,
     i           t0BC, u0BC, w0BC, qBC,
     i           tnW0, roW0, pW0, uT0, HDensW0,
     o           UThat, URhat, TWhat, pWhat, roWhat )
      implicit none
      integer    NModes, nzWT, nzW, IdzT, qBC
      integer    izT(nzWT+1)
      real*8     Omega,  OmegaIm, Omega0, KT
      real*8     t0BC, u0BC, w0BC
      complex*16 eT(2*NModes,nzWT)
      logical    DimLessStVct, DoScaling
!
      real*8 tnW0(nzW)
      real*8 roW0(nzW)
      real*8 pW0(nzW)
      real*8 uT0(nzW)
      real*8 HDensW0(nzW)
c
c     --- outputs ---
      complex*16 UThat(nzWT)
      complex*16 URhat(nzWT)
      complex*16 TWhat(nzWT)
      complex*16 pWhat(nzWT)
      complex*16 roWhat(nzWT)
c
c     --- local variables ---
      integer    iT, i, j
      real*8     kX, FctDimLessStVct, Src0
      complex*16 zero, im
      complex*16 U, W, T, Wcal
      complex*16 OmegaC, Omega0C, OmegaDop, fct
      complex*16, allocatable:: EVct(:)
c
      kX  = 1.d-3*kT  !kx in 1/m
C
      OmegaC  = dcmplx( Omega,  OmegaIm ) !Omega
      Omega0C = dcmplx( Omega0, 0.d0 )    !Omega0
c
      im   = (0.d0,1.d0)
      zero = (0.d0,0.d0)
c
c     --- factor for obtaining a dimensionless system of equations ---
      if ( DimLessStVct ) then
        FctDimLessStVct = kX
      else
        FctDimLessStVct = 1.d0
      end if
c
      allocate( EVct( 2*NModes ) )
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        OmegaDop = OmegaC - dcmplx(kX * uT0(i),0.D0)
        fct = Omega0C / OmegaDop
!       -----------------------------------------------------------------------
!                                  e(2*NModes)
!       -----------------------------------------------------------------------
        do j = 1, 2*NModes
          EVct(j) = eT(j,iT)
        end do
!
        U  = EVct(1)
        W  = EVct(2)
        T  = EVct(3)
        Wcal = FctDimLessStVct * EVct(5)
!
        UThat(iT)  = Omega0C * U / kX !in m/s
        URhat(iT)  = Omega0C * W / kX !in m/s
        TWhat(iT)  = tnW0(i) * T      !in K
c
        roWhat(iT) = roW0(i) * ( fct * U
     &             - im * fct * W / (kX * HDensW0(i))
     &             + im * fct * Wcal / kX ) !in Kg/m**3
        pWhat(iT) = pW0(i) * ( TWhat(iT)/tnW0(i) + roWhat(iT)/roW0(i) ) !in N/m**2
      end do
!
!     --- scale quantities ---
      if ( DoScaling ) then
        if ( qBC .eq. 1 ) then
          Src0 = u0BC / abs( dble(UThat(1)) )
        else if ( qBC .eq. 2 ) then
          Src0 = w0BC / abs( dble(URhat(1)) )
        else
          Src0 = t0BC / abs( dble(TWhat(1)) )
        end if
        do iT = 1, nzWT
          UThat(iT)  = Src0 * UThat(iT)
          URhat(iT)  = Src0 * URhat(iT)
          TWhat(iT)  = Src0 * TWhat(iT)
          roWhat(iT) = Src0 * roWhat(iT)
          pWhat(iT)  = Src0 * pWhat(iT)
        end do
      end if
!
      deallocate( EVct )
      end
c     *************************************************************************
c
C
C
c
c                             TypeSolMet = 3
c               Scattering-Matrix Method for Amplitudes (SMMA)
c
c
C
C
c     *************************************************************************
      subroutine WaveAmplitudesModelSMMA
     i         ( TypeSclLayEq, TypeBC, qBC, NModes, nzWT, kT, zT,
     i           EValT, EVctT,
     o           AmpT )
c------------------------------------------------------------------------------
c
c             o zT(nzW+1), et(nzW+1)
c             |
c             |   layer iT = nzW, AmpT(*,nzW), EVal(*,nzW)
c             |
c BC = nzW-1  o zT(nzW), et(nzW): ALay1(nzW-1)*AmpT(nzW)-ALay0(nzW-1)*AmpT(nzW-1) = 0
c             |
c             |  layer iT = nzW-1, AmpT(*,nzW-1), EVal(*,nzW-1)
c             |
c             o zT(nzW-1), et(nzW-1)
c
c            ...
c
c             o zT(iT+2), et(iT+2)
c             |
c             |   layer iT + 1, AT(*,iT+1), EVal(*,iT+1)
c             |
c    BC = iT  o zT(iT+1), et(iT+1): ALay1(iT)*AmpT(iT+1) - ALay0(iT)*AmpT(iT) = 0
c             |
c             |   layer iT, AmpT(*,iT), EVal(*,iT); dz = zT(iT+1)-zT(iT)
c             |
c             o zT(iT), et(IT)
c
c            ...
c
c             o zT(3), et(3)
c             |
c             |   layer iT = 2, AT(*,2), EVal(*,2)
c             |
c      BC = 1 o zT(2), et(2):ALay1(1)*AmpT(2) - ALay0(1)*AmpT(1) = 0
c             |
c             |   layer iT = 1, AmpT(*,1), EVal(*,1)
c             |
c             o zT(1) = zmin, et(1)
c
c
c     Unknowns
c          AT(1),...,AT(nzW) = AT(Nlevels)
c     Layers equations
c          iT = 1, ...,nzW-1
c
c     Calling routines
c        - RrondTrondMatricesUp_ModelSMMA
c        - RrondTrondMatricesDown_ModelSMMA
c        - ZGETRF, ZGETRI
c------------------------------------------------------------------------------
      implicit none
      integer    TypeSclLayEq, TypeBC, qBC, NModes, nzWT
      real*8     KT
      real*8     zT(nzWT+1)
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
!
      complex*16  AmpT(2*NModes,nzWT)
!
!     --- local variables ---
      integer    iT, i, j, k
      complex*16 one, zero, sum
!
      complex*16, allocatable:: RrondUp0(:,:)
      complex*16, allocatable:: RrondUp1(:,:)
      complex*16, allocatable:: TrondUp0(:,:)
      complex*16, allocatable:: TrondUp1(:,:)
!
      complex*16, allocatable:: RrondDown0(:,:)
      complex*16, allocatable:: RrondDown1(:,:)
      complex*16, allocatable:: TrondDown0(:,:)
      complex*16, allocatable:: TrondDown1(:,:)
      complex*16, allocatable:: Ti(:,:)
!
      complex*16, allocatable:: aTPlus(:,:)
      complex*16, allocatable:: aTMinus(:,:)
      complex*16, allocatable:: b( : )
!
      integer  info
      logical  PrintInfo
      integer, allocatable:: ipiv(:)
      complex*16, allocatable:: work(:)
c
      PrintInfo = .false.
!
      allocate( RrondUp0(NModes,NModes) )
      allocate( RrondUp1(NModes,NModes) )
      allocate( TrondUp0(NModes,NModes) )
      allocate( TrondUp1(NModes,NModes) )
!
      allocate( RrondDown0(NModes,NModes) )
      allocate( RrondDown1(NModes,NModes) )
      allocate( TrondDown0(NModes,NModes) )
      allocate( TrondDown1(NModes,NModes) )

      allocate( TI(NModes,NModes) )
!
      allocate( aTPlus(NModes,nzWT) )
      allocate( aTMinus(NModes,nzWT) )
      allocate( b( NModes ) )
!
      allocate( work( NModes ), ipiv( NModes ) )
!
      one  = (1.d0,0.d0)
      zero = (0.d0,0.d0)
!
      call RrondTrondMatricesUp_ModelSMMA
     i   ( nzWT, TypeSclLayEq, NModes, nzWT, kT, zT, EValT, EVctT,
     o     RrondUp0, RrondUp1, TrondUp0, TrondUp1 )
!
!    --- amplitudes at top and bottom from boundary conditions ---
      if ( TypeBC .eq. 1 ) then
        b = zero
        b(1) = one
      else
        call BoundaryConditionComponent
     i     ( qBC, NModes, nzWT, kT, EValT, EVctT,
     o       b )
      end if
      do i = 1, NModes
        aTPlus(i,1) = b(i)
        aTMinus(i,nzWT) = zero
      end do
!
!     --- rest of amplitudes at top and bottom ---
      do i = 1, NModes
        sum = zero
        do j = 1, NModes
          sum = sum + RrondUp0(i,j) * aTPlus(j,1)
        end do
        aTMinus(i,1) = sum
c
        sum = zero
        do j = 1, NModes
          sum = sum + TrondUp1(i,j) * aTPlus(j,1)
        end do
        aTPlus(i,nzWT)  = sum
      end do
!
!     --- internal amplitudes ---
      do iT = 2, nzWT - 1
        call RrondTrondMatricesUp_ModelSMMA
     i     ( iT, TypeSclLayEq, NModes, nzWT, kT, zT, EValT, EVctT,
     o       RrondUp0, RrondUp1, TrondUp0, TrondUp1 )
!
        call RrondTrondMatricesDown_ModelSMMA
     i     ( iT, TypeSclLayEq, NModes, nzWT, kT, zT, EValT, EVctT,
     o       RrondDown0, RrondDown1, TrondDown0, TrondDown1 )
!
!       --- TI = I - RrondUp1*RrondDown0 ---
        do i = 1, NModes
          do j = 1, NModes
            sum = zero
            do k = 1, NModes
              sum = sum + RrondUp1(i,k) * RrondDown0(k,j)
            end do
            TI(i,j) = - sum
          end do
          TI(i,i) = TI(i,i) + one
        end do
        call ZGETRF( NModes, NModes, TI, NModes, ipiv, info)
        if ( info .ne. 0 ) then
          write( 6, * ) 'Error in LU factorization routine DGETRF'
          STOP
        end if
        call ZGETRI( NModes, TI, NModes, ipiv, work, NModes, info)
        if ( info .ne. 0 ) then
          write( 6, * ) 'Error in LU substitution routine DGETRI'
          STOP
        end if
!
!       ---aTPlus(i,iT) = TI(i,j) * TrondUp1(j,k) * aTPlus(k,1)) )
        do i = 1, NModes
          do j = 1, NModes
            sum = zero
            do k = 1, NModes
              sum = sum + TrondUp1(j,k) * aTPlus(k,1)
            end do
            b(j) = sum
          end do
          sum = zero
          do j = 1, NModes
            sum = sum + TI(i,j) * b(j)
          end do
          aTPlus(i,iT) = sum
        end do
!
!       ---aTMinus(i,iT) = RrondDown0(i,j)*aTPlus(j,iT)  ---
        do i = 1, NModes
          sum = zero
          do j = 1, NModes
            sum = sum + RrondDown0(i,j) * aTPlus(j,iT)
          end do
          aTMinus(i,iT) = sum
        end do
      end do
!
!     --- AmpT (*,iT) = [aTPlus(*,iT), aTMinus(*,iT)]^T ---
      do iT = 1, nzWT
        do i = 1, NModes
          AmpT(i,iT)        = aTPlus(i,iT)
          AmpT(i+NModes,iT) = aTMinus(i,iT)
        end do
      end do
!
      if ( PrintInfo ) then
        write(6,'(a)') ' Amplitudes at bottom boundary'
        do i = 1, 2*NModes
          write(6,'(a,i2,a,2(1pe11.4,x))')
     &  ' Mode: ', i, ', Amplitude = ',  AmpT(i,1)
        end do
        write(6,'(a)') ' Amplitudes at top boundary'
        do i = 1, 2*NModes
          write(6,'(a,i2,a,2(1pe11.4,x))')
     &  ' Mode: ', i, ', Amplitude = ',  AmpT(i,nzWT)
        end do
      end if
!
      deallocate( RrondUp0 )
      deallocate( RrondUp1 )
      deallocate( TrondUp0 )
      deallocate( TrondUp1 )

      deallocate( RrondDown0 )
      deallocate( RrondDown1 )
      deallocate( TrondDown0 )
      deallocate( TrondDown1 )
      deallocate( TI )
!
      deallocate( aTPlus )
      deallocate( aTMinus )
      deallocate( b )
!
      deallocate( work, ipiv )
      end
c     *************************************************************************
      subroutine WaveAmplitudesModelSMMA1
     i         ( TypeSclLayEq, NModes, nzWT, kT, zT,
     i           EValT, EVctT,
     o           AmpT )
c------------------------------------------------------------------------------
c
c             o zT(nzW+1), et(nzW+1)
c             |
c             |   layer iT = nzW, AmpT(*,nzW), EVal(*,nzW)
c             |
c BC = nzW-1  o zT(nzW), et(nzW): ALay1(nzW-1)*AmpT(nzW)-ALay0(nzW-1)*AmpT(nzW-1) = 0
c             |
c             |  layer iT = nzW-1, AmpT(*,nzW-1), EVal(*,nzW-1)
c             |
c             o zT(nzW-1), et(nzW-1)
c
c            ...
c
c             o zT(iT+2), et(iT+2)
c             |
c             |   layer iT + 1, AT(*,iT+1), EVal(*,iT+1)
c             |
c    BC = iT  o zT(iT+1), et(iT+1): ALay1(iT)*AmpT(iT+1) - ALay0(iT)*AmpT(iT) = 0
c             |
c             |   layer iT, AmpT(*,iT), EVal(*,iT); dz = zT(iT+1)-zT(iT)
c             |
c             o zT(iT), et(IT)
c
c            ...
c
c             o zT(3), et(3)
c             |
c             |   layer iT = 2, AT(*,2), EVal(*,2)
c             |
c      BC = 1 o zT(2), et(2):ALay1(1)*AmpT(2) - ALay0(1)*AmpT(1) = 0
c             |
c             |   layer iT = 1, AmpT(*,1), EVal(*,1)
c             |
c             o zT(1) = zmin, et(1)
c
c
c     Unknowns
c          AT(1),...,AT(nzW) = AT(Nlevels)
c     Layers equations
c          iT = 1, ...,nzW-1
c
c     Calling routines
c        - RrondTrondMatricesUp_ModelSMMA
c        - RrondTrondMatricesDown_ModelSMMA
c        - ZGETRF, ZGETRI
c------------------------------------------------------------------------------
      implicit none
      integer    TypeSclLayEq, NModes, nzWT
      real*8     KT
      real*8     zT(nzWT+1)
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
!
      complex*16  AmpT(2*NModes,nzWT)
!
!     --- local variables ---
      integer    iT, i, j, k
      complex*16 one, zero, sum
!
      complex*16, allocatable:: RrondUp0(:,:)
      complex*16, allocatable:: RrondUp1(:,:)
      complex*16, allocatable:: TrondUp0(:,:)
      complex*16, allocatable:: TrondUp1(:,:)
!
      complex*16, allocatable:: RrondDown0(:,:)
      complex*16, allocatable:: RrondDown1(:,:)
      complex*16, allocatable:: TrondDown0(:,:)
      complex*16, allocatable:: TrondDown1(:,:)
      complex*16, allocatable:: Ti(:,:)
!
      complex*16, allocatable:: aTPlus(:,:)
      complex*16, allocatable:: aTMinus(:,:)
!
      integer  info
      logical  PrintInfo
      integer, allocatable:: ipiv(:)
      complex*16, allocatable:: work(:)
c
      PrintInfo = .false.
!
      allocate( RrondUp0(NModes,NModes) )
      allocate( RrondUp1(NModes,NModes) )
      allocate( TrondUp0(NModes,NModes) )
      allocate( TrondUp1(NModes,NModes) )
!
      allocate( RrondDown0(NModes,NModes) )
      allocate( RrondDown1(NModes,NModes) )
      allocate( TrondDown0(NModes,NModes) )
      allocate( TrondDown1(NModes,NModes) )

      allocate( TI(NModes,NModes) )
!
      allocate( aTPlus(NModes,nzWT) )
      allocate( aTMinus(NModes,nzWT) )
!
      allocate( work( NModes ), ipiv( NModes ) )
!
      one  = (1.d0,0.d0)
      zero = (0.d0,0.d0)
!
      call RrondTrondMatricesUp_ModelSMMA
     i   ( nzWT, TypeSclLayEq, NModes, nzWT, kT, zT, EValT, EVctT,
     o     RrondUp0, RrondUp1, TrondUp0, TrondUp1 )
!
!    --- amplitudes at top and bottom from boundary conditions ---
      do i = 1, NModes
        if ( i .eq. 1 ) then
          aTPlus(i,1) = one
        else
          aTPlus(i,1) = zero
        end if
        aTMinus(i,nzWT) = zero
      end do
!
!     --- rest of amplitudes at top and bottom ---
      do i = 1, NModes
        aTMinus(i,1)    = RrondUp0(i,1)
        aTPlus(i,nzWT)  = TrondUp1(i,1)
      end do
!
!     --- internal amplitudes ---
      do iT = 2, nzWT - 1
        call RrondTrondMatricesUp_ModelSMMA
     i     ( iT, TypeSclLayEq, NModes, nzWT, kT, zT, EValT, EVctT,
     o       RrondUp0, RrondUp1, TrondUp0, TrondUp1 )
!
        call RrondTrondMatricesDown_ModelSMMA
     i     ( iT, TypeSclLayEq, NModes, nzWT, kT, zT, EValT, EVctT,
     o       RrondDown0, RrondDown1, TrondDown0, TrondDown1 )
!
!       --- TI = I - RrondUp1*RrondDown0 ---
        do i = 1, NModes
          do j = 1, NModes
            sum = zero
            do k = 1, NModes
              sum = sum + RrondUp1(i,k) * RrondDown0(k,j)
            end do
            TI(i,j) = - sum
          end do
          TI(i,i) = TI(i,i) + one
        end do
        call ZGETRF( NModes, NModes, TI, NModes, ipiv, info)
        if ( info .ne. 0 ) then
          write( 6, * ) 'Error in LU factorization routine DGETRF'
          STOP
        end if
        call ZGETRI( NModes, TI, NModes, ipiv, work, NModes, info)
        if ( info .ne. 0 ) then
          write( 6, * ) 'Error in LU substitution routine DGETRI'
          STOP
        end if
!
!       ---aTPlus(*,iT) = TI * TrondUp1 * aTPlus(*,1)) )
        do i = 1, NModes
          sum = zero
          do j = 1, NModes
            sum = sum + TI(i,j) * TrondUp1(j,1)
          end do
          aTPlus(i,iT) = sum
        end do
!
!       ---aTMinus(*,iT) = RrondDown0*aTPlus(*,iT)  ---
        do i = 1, NModes
          sum = zero
          do j = 1, NModes
            sum = sum + RrondDown0(i,j) * aTPlus(j,iT)
          end do
          aTMinus(i,iT) = sum
        end do
      end do
!
!     --- AmpT (*,iT) = [aTPlus(*,iT), aTMinus(*,iT)]^T ---
      do iT = 1, nzWT
        do i = 1, NModes
          AmpT(i,iT)        = aTPlus(i,iT)
          AmpT(i+NModes,iT) = aTMinus(i,iT)
        end do
      end do
!
      if ( PrintInfo ) then
        write(6,'(a)') ' Amplitudes at bottom boundary'
        do i = 1, 2*NModes
          write(6,'(a,i2,a,2(1pe11.4,x))')
     &  ' Mode: ', i, ', Amplitude = ',  AmpT(i,1)
        end do
        write(6,'(a)') ' Amplitudes at top boundary'
        do i = 1, 2*NModes
          write(6,'(a,i2,a,2(1pe11.4,x))')
     &  ' Mode: ', i, ', Amplitude = ',  AmpT(i,nzWT)
        end do
      end if
!
      deallocate( RrondUp0 )
      deallocate( RrondUp1 )
      deallocate( TrondUp0 )
      deallocate( TrondUp1 )

      deallocate( RrondDown0 )
      deallocate( RrondDown1 )
      deallocate( TrondDown0 )
      deallocate( TrondDown1 )
      deallocate( TI )
!
      deallocate( aTPlus )
      deallocate( aTMinus )
!
      deallocate( work, ipiv )
      end
!     *************************************************************************
      subroutine RrondTrondMatricesUp_ModelSMMA
     i         ( iT0, TypeSclLayEq, NModes, nzWT, kT, zT, EValT, EVctT,
     o           Rrond0, Rrond1, Trond0, Trond1 )
C     -------------------------------------------------------------------------
C      Calling routines
C        - LayerMatricesModelGMMA
c        - RTmatrices_ModelSMMA
C        - InteractionPrincipleModelSMMA
C     -------------------------------------------------------------------------
      implicit none
      integer    iT0, TypeSclLayEq, NModes, nzWT
      real*8     KT
      real*8     zT(nzWT+1)
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
!
      complex*16 Rrond0(NModes,NModes)
      complex*16 Rrond1(NModes,NModes)
      complex*16 Trond0(NModes,NModes)
      complex*16 Trond1(NModes,NModes)
!
!     --- local variables ---
      integer    iT, i, j
      real*8     dz
!
      complex*16, allocatable:: ALayer1(:,:)
      complex*16, allocatable:: ALayer0(:,:)
!
      complex*16, allocatable:: R0(:,:)
      complex*16, allocatable:: R1(:,:)
      complex*16, allocatable:: T0(:,:)
      complex*16, allocatable:: T1(:,:)

      allocate( ALayer1(2*NModes,2*NModes) )
      allocate( ALayer0(2*NModes,2*NModes) )
!
      allocate( R0(NModes,NModes) )
      allocate( R1(NModes,NModes) )
      allocate( T0(NModes,NModes) )
      allocate( T1(NModes,NModes) )
!
      do iT = 1, iT0 - 1
        dz = zT(iT+1) - zT(iT)
!       -----------------------------------------------------------------------
!                           Layer matrices
!       -----------------------------------------------------------------------
        call LayerMatricesModelGMMA
     i     ( iT, TypeSclLayEq, NModes, nzWT, dz, kT,
     i       EValT, EVctT,
     o       ALayer1, ALayer0 )
!       .......................................................................
!        Reflection and transmission matrices R0, R1, T0, T1 for this layer
!       .......................................................................
        call RTmatrices_ModelSMMA
     i     ( NModes, ALayer0, ALayer1,
     o       R0, R1, T0, T1 )
!       .......................................................................
!                Initialize Rrond0, Rrond1, Trond0, Trond1 for iT = 1
!       .......................................................................
        if ( iT .eq. 1 ) then
          do i = 1, NModes
            do j = 1, NModes
              Rrond0(i,j) = R0(i,j)
              Rrond1(i,j) = R1(i,j)
!
              Trond0(i,j) = T0(i,j)
              Trond1(i,j) = T1(i,j)
            end do
          end do
!       ........................................................................
!                       Recurrence for iT = 2,...,iT0 - 1
!       ........................................................................
        else
          call InteractionPrincipleModelSMMA
     i       ( NModes, R0, R1, T0, T1,
     o         Rrond0, Rrond1, Trond0, Trond1 )
        end if
      end do
!
      deallocate( ALayer1 )
      deallocate( ALayer0 )
!
      deallocate( R0 )
      deallocate( R1 )
      deallocate( T0 )
      deallocate( T1 )
      end
!     *************************************************************************
      subroutine RrondTrondMatricesDown_ModelSMMA
     i         ( iT0, TypeSclLayEq, NModes, nzWT, kT, zT, EValT, EVctT,
     o           Rrond0, Rrond1, Trond0, Trond1 )
C     -------------------------------------------------------------------------
C      Calling routines
C        - LayerMatricesModelGMMA
c        - RTmatrices_ModelSMMA
C        - InteractionPrincipleModelSMMA
C     -------------------------------------------------------------------------
      implicit none
      integer    iT0, TypeSclLayEq, NModes, nzWT
      real*8     KT
      real*8     zT(nzWT+1)
      complex*16 EValT(2*NModes, nzWT)
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
!
      complex*16 Rrond0(NModes,NModes)
      complex*16 Rrond1(NModes,NModes)
      complex*16 Trond0(NModes,NModes)
      complex*16 Trond1(NModes,NModes)
!
!     --- local variables ---
      integer    iT, i, j
      real*8     dz
!
      complex*16, allocatable:: ALayer1(:,:)
      complex*16, allocatable:: ALayer0(:,:)
!
      complex*16, allocatable:: R0(:,:)
      complex*16, allocatable:: R1(:,:)
      complex*16, allocatable:: T0(:,:)
      complex*16, allocatable:: T1(:,:)

      allocate( ALayer1(2*NModes,2*NModes) )
      allocate( ALayer0(2*NModes,2*NModes) )
!
      allocate( R0(NModes,NModes) )
      allocate( R1(NModes,NModes) )
      allocate( T0(NModes,NModes) )
      allocate( T1(NModes,NModes) )
!
      do iT = iT0, nzWT - 1
        dz = zT(iT+1) - zT(iT)
!       -----------------------------------------------------------------------
!                           Layer matrices
!       -----------------------------------------------------------------------
        call LayerMatricesModelGMMA
     i     ( iT, TypeSclLayEq, NModes, nzWT, dz, kT,
     i       EValT, EVctT,
     o       ALayer1, ALayer0 )
!       .......................................................................
!        Reflection and transmission matrices R0, R1, T0, T1 for this layer
!       .......................................................................
        call RTmatrices_ModelSMMA
     i     ( NModes, ALayer0, ALayer1,
     o       R0, R1, T0, T1 )
!       .......................................................................
!            Initialize Rrond0, Rrond1, Trond0, Trond1 for iT = 1
!       .......................................................................
        if ( iT .eq. iT0 ) then
          do i = 1, NModes
            do j = 1, NModes
              Rrond0(i,j) = R0(i,j)
              Rrond1(i,j) = R1(i,j)
!
              Trond0(i,j) = T0(i,j)
              Trond1(i,j) = T1(i,j)
            end do
          end do
!       .......................................................................
!                 Recurrence for iT = iT0+1,...,nzWT - 1
!       .......................................................................
        else
          call InteractionPrincipleModelSMMA
     i       ( NModes, R0, R1, T0, T1,
     o         Rrond0, Rrond1, Trond0, Trond1 )
        end if
      end do
!
      deallocate( ALayer1 )
      deallocate( ALayer0 )
!
      deallocate( R0 )
      deallocate( R1 )
      deallocate( T0 )
      deallocate( T1 )
      end
!     **************************************************************************
      subroutine RTmatrices_ModelSMMA
     i         ( NModes, ALayer0, ALayer1,
     o           R0, R1, T0, T1 )
C     -------------------------------------------------------------------------
C      Calling routines
C        - ZGETRF
C        - ZGETRI
C     -------------------------------------------------------------------------
      implicit none
      integer    NModes
      complex*16 ALayer1(2*NModes,2*NModes)
      complex*16 ALayer0(2*NModes,2*NModes)
!
      complex*16 R0(NModes,NModes)
      complex*16 R1(NModes,NModes)
      complex*16 T0(NModes,NModes)
      complex*16 T1(NModes,NModes)
!
!     --- local variables ---
      integer     i, j, k
      complex*16  zero, sum
      complex*16, allocatable:: x(:,:), y(:,:), z(:,:)
!
      integer  info
      integer, allocatable:: ipiv(:)
      complex*16, allocatable:: work(:)
!
      zero = (0.d0,0.d0)
!     -------------------------------------------------------------------------
!                    Compute x(2*NModes,2*NModes) and y(2*NModes,2*NModes)
!     -------------------------------------------------------------------------
      allocate( x(2*NModes,2*NModes), y(2*NModes,2*NModes) )
      do i = 1, NModes
        do j = 1, NModes
          x(i,j)               =   Alayer0(i,j+NModes)
          x(i,j+NModes)        = - Alayer1(i,j)
          x(i+NModes,j)        =   Alayer0(i+NModes,j+NModes)
          x(i+NModes,j+NModes) = - Alayer1(i+NModes,j)
!
          y(i,j)               = - Alayer0(i,j)
          y(i,j+NModes)        =   Alayer1(i,j+NModes)
          y(i+NModes,j)        = - Alayer0(i+NModes,j)
          y(i+NModes,j+NModes) =   Alayer1(i+NModes,j+NModes)
        end do
      end do
!     ------------------------------------------------------------------------
!                            Compute X^(-1)
!     ------------------------------------------------------------------------
      allocate( work( 2*NModes ), ipiv( 2*NModes ) )
      call ZGETRF( 2*NModes, 2*NModes, x, 2*NModes, ipiv, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU factorization routine DGETRF'
        STOP
      end if
      call ZGETRI( 2*NModes, x, 2*NModes, ipiv, work, 2*NModes, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU substitution routine DGETRI'
        STOP
      end if
C     ------------------------------------------------------------------------
C                               Compute Z = X^(-1)*Y
C     ------------------------------------------------------------------------
      allocate( z(2*NModes,2*NModes) )
      do i = 1, 2*NModes
        do j = 1, 2*NModes
          sum = zero
          do k = 1, 2*NModes
            sum = sum + x(i,k)*y(k,j)
          end do
          z(i,j) = sum
        end do
      end do
!
      do i = 1, NModes
        do j = 1, NModes
          R0(i,j) = z(i,j)
          R1(i,j) = z(i+NModes,j+NModes)
          T0(i,j) = z(i,j+NModes)
          T1(i,j) = z(i+NModes,j)
        end do
      end do
!
      deallocate( x, y, z )
      deallocate( work, ipiv )
      end
!     **************************************************************************
      subroutine InteractionPrincipleModelSMMA
     i         ( NModes, R0, R1, T0, T1,
     o           Rrond0, Rrond1, Trond0, Trond1 )
C     -------------------------------------------------------------------------
C      Calling routines
C        - ZGETRF, ZGETRI
C     -------------------------------------------------------------------------
      implicit none
      integer   NModes
      complex*16 R0(NModes,NModes)
      complex*16 R1(NModes,NModes)
      complex*16 T0(NModes,NModes)
      complex*16 T1(NModes,NModes)

!
      complex*16 Rrond0(NModes,NModes)
      complex*16 Rrond1(NModes,NModes)
      complex*16 Trond0(NModes,NModes)
      complex*16 Trond1(NModes,NModes)
!
!     --- local variables ---
      integer     i, j, k
      complex*16  one, zero, sum
      complex*16, allocatable:: R0Rrond1(:,:), Rrond1R0(:,:)
      complex*16, allocatable:: a(:,:), b(:,:)
!
      integer  info
      integer, allocatable:: ipiv(:)
      complex*16, allocatable:: work(:)
!
      one  = (1.d0,0.d0)
      zero = (0.d0,0.d0)
!
      allocate( R0Rrond1(NModes,NModes) )
      allocate( Rrond1R0(NModes,NModes) )
      allocate( a(NModes, NModes) )
      allocate( b(NModes, NModes) )
!
      allocate( work( NModes ), ipiv( NModes ) )
!     -----------------------------------------------------------------------
!                            R0Rrond1 = I - R0 * Rrond1
!     -----------------------------------------------------------------------
      do i = 1, NModes
        do j = 1, NModes
          sum = zero
          do k = 1, NModes
            sum = sum + R0(i,k) * Rrond1(k,j)
          end do
          R0Rrond1(i,j) = - sum
        end do
        R0Rrond1(i,i) = R0Rrond1(i,i) + one
      end do
!     -----------------------------------------------------------------------
!                            Rrond1R0 = I - Rrond1 * R0
!     -----------------------------------------------------------------------
      do i = 1, NModes
        do j = 1, NModes
          sum = zero
          do k = 1, NModes
            sum = sum + Rrond1(i,k) * R0(k,j)
          end do
          Rrond1R0(i,j) = - sum
        end do
        Rrond1R0(i,i) = Rrond1R0(i,i) + one
      end do
!     -----------------------------------------------------------------------
!                       R0Rrond1^(-1) = (I - R0 * Rrond1)^(-1)
!     -----------------------------------------------------------------------
      call ZGETRF( NModes, NModes, R0Rrond1, NModes, ipiv, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU factorization routine DGETRF:IP1'
        STOP
      end if
      call ZGETRI( NModes, R0Rrond1, NModes, ipiv, work, NModes, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU substitution routine DGETRI'
        STOP
      end if
!     -----------------------------------------------------------------------
!                       Rrond1R0^(-1) = (I - Rrond1 * R0)^(-1)
!     -----------------------------------------------------------------------
      ipiv = 0
      work = zero
      call ZGETRF( NModes, NModes, Rrond1R0, NModes, ipiv, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU factorization routine DGETRF:IP2'
        STOP
      end if
      call ZGETRI( NModes, Rrond1R0, NModes, ipiv, work, NModes, info)
      if ( info .ne. 0 ) then
        write( 6, * ) 'Error in LU substitution routine DGETRI'
        STOP
      end if
!     -----------------------------------------------------------------------
!                         A = Trond0 * R0Rrond1^(-1)
!     -----------------------------------------------------------------------
      do i = 1, NModes
        do j = 1, NModes
          sum = zero
          do k = 1, NModes
            sum = sum + Trond0(i,k) * R0Rrond1(k,j)
          end do
          a(i,j) = sum
        end do
      end do
!     -----------------------------------------------------------------------
!                            B = R0 * Trond1
!     -----------------------------------------------------------------------
      do i = 1, NModes
        do j = 1, NModes
          sum = zero
          do k = 1, NModes
            sum = sum + R0(i,k) * Trond1(k,j)
          end do
          b(i,j) = sum
        end do
      end do
!     -------------------------------------------------------------------------
!    UPDATE Rrond0 = Rrond0 + Trond0*R0Rrond1^(-1) * R0*Trond1 = Rrond0 + A * B
!     -------------------------------------------------------------------------
      do i = 1, NModes
        do j = 1, NModes
          sum = zero
          do k = 1, NModes
            sum = sum + A(i,k) * B(k,j)
          end do
          Rrond0(i,j) = Rrond0(i,j) + sum
        end do
      end do
!     -----------------------------------------------------------------------
!            UPDATE Trond0 = Trond0*R0Rrond1^(-1) * T0 = A * T0
!     -----------------------------------------------------------------------
      do i = 1, NModes
        do j = 1, NModes
          sum = zero
          do k = 1, NModes
            sum = sum + A(i,k) * T0(k,j)
          end do
          Trond0(i,j) = sum
        end do
      end do
!     -----------------------------------------------------------------------
!                         A = T1 * Rrond1R0^(-1)
!     -----------------------------------------------------------------------
      do i = 1, NModes
        do j = 1, NModes
          sum = zero
          do k = 1, NModes
            sum = sum + T1(i,k) * Rrond1R0(k,j)
          end do
          a(i,j) = sum
        end do
      end do
!     -----------------------------------------------------------------------
!                            B = Rrond1 * T0
!     -----------------------------------------------------------------------
      do i = 1, NModes
        do j = 1, NModes
          sum = zero
          do k = 1, NModes
            sum = sum + Rrond1(i,k) * T0(k,j)
          end do
          b(i,j) = sum
        end do
      end do
!     -----------------------------------------------------------------------
!     UPDATE Rrond1 = R1 + T1*Rrond1R0^(-1) * Rrond1*T0 = R1 + A * B
!     -----------------------------------------------------------------------
      do i = 1, NModes
        do j = 1, NModes
          sum = zero
          do k = 1, NModes
            sum = sum + A(i,k) * B(k,j)
          end do
          Rrond1(i,j) = R1(i,j) + sum
        end do
      end do
!     -----------------------------------------------------------------------
!     UPDATE Trond1 = T1*Rrond1R0^(-1) * Trond1 = A * Trond1
!     -----------------------------------------------------------------------
      do i = 1, NModes
        do j = 1, NModes
          sum = zero
          do k = 1, NModes
            sum = sum + A(i,k) * Trond1(k,j)
          end do
          b(i,j) = sum
        end do
      end do
!
      do i = 1, NModes
        do j = 1, NModes
          Trond1(i,j) = b(i,j)
        end do
      end do
!
      deallocate( R0Rrond1 )
      deallocate( Rrond1R0 )
      deallocate( a )
      deallocate( b )
!
      deallocate( work, ipiv )
      end
!     *************************************************************************
      subroutine WaveParametersModelSMMA
     i         ( NModes, nzWT, nzW, IdzT, izT,
     i           Omega, OmegaIm, Omega0, kT,
     i           EVctT, AmpT,
     i           DimLessStVct, DoScaling,
     i           t0BC, u0BC, w0BC, qBC,
     i           tnW0, roW0, pW0,  uT0, HDensW0,
     o           UThat, URhat, TWhat, pWhat, roWhat )
      implicit none
      integer    NModes, nzWT, nzW, IdzT, qBC
      integer    izT(nzWT+1)
      real*8     Omega, OmegaIm, Omega0, KT
      real*8     t0BC, u0BC, w0BC
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
      complex*16 AmpT(2*NModes,nzWT)
      logical    DimLessStVct, DoScaling
!
      real*8 tnW0(nzW)
      real*8 roW0(nzW)
      real*8 pW0(nzW)
      real*8 uT0(nzW)
      real*8 HDensW0(nzW)
c
c     --- outputs ---
      complex*16 UThat(nzWT)
      complex*16 URhat(nzWT)
      complex*16 TWhat(nzWT)
      complex*16 pWhat(nzWT)
      complex*16 roWhat(nzWT)
c
c     --- local variables ---
      integer    iT, i, j, k
      real*8     kX, FctDimLessStVct, Src0
      complex*16 zero, sum, im
      complex*16 U, W, T, Wcal
      complex*16 OmegaC, Omega0C, OmegaDop, fct
      complex*16, allocatable:: EVct(:)
c
      kX  = 1.d-3*kT  !kx in 1/m
C
      OmegaC  = dcmplx( Omega,  OmegaIm ) !Omega
      Omega0C = dcmplx( Omega0, 0.d0 )    !Omega0
c
      im   = (0.d0,1.d0)
      zero = (0.d0,0.d0)
c
c     --- factor for obtaining a dimensionless system of equations ---
      if ( DimLessStVct ) then
        FctDimLessStVct = kX
      else
        FctDimLessStVct = 1.d0
      end if
c
      allocate( EVct( 2*NModes ) )
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        OmegaDop = OmegaC - dcmplx(kX * uT0(i),0.D0)
        fct = Omega0C / OmegaDop
c
        do j = 1, 2*NModes
          sum = zero
          do k = 1, 2*NModes
            sum = sum + EVctT(j,k,iT) * AmpT(k,iT)
          end do
          EVct(j) = sum
        end do
c
        U  = EVct(1)
        W  = EVct(2)
        T  = EVct(3)
        Wcal = FctDimLessStVct * EVct(5)
!
        UThat(iT)  = Omega0C * U / kX !in m/s
        URhat(iT)  = Omega0C * W / kX !in m/s
        TWhat(iT)  = tnW0(i) * T      !in K
c
        roWhat(iT) = roW0(i) * ( fct * U
     &             - im * fct * W / (kX * HDensW0(i))
     &             + im * fct * Wcal / kX ) !in Kg/m**3
        pWhat(iT) = pW0(i) * ( TWhat(iT)/tnW0(i) + roWhat(iT)/roW0(i) ) !in N/m**2
      end do
!
!     --- scale quantities ---
      if ( DoScaling ) then
        if ( qBC .eq. 1 ) then
          Src0 = u0BC / abs( dble(UThat(1)) )
        else if ( qBC .eq. 2 ) then
          Src0 = w0BC / abs( dble(URhat(1)) )
        else
          Src0 = t0BC / abs( dble(TWhat(1)) )
        end if
        do iT = 1, nzWT
          UThat(iT)  = Src0 * UThat(iT)
          URhat(iT)  = Src0 * URhat(iT)
          TWhat(iT)  = Src0 * TWhat(iT)
          roWhat(iT) = Src0 * roWhat(iT)
          pWhat(iT)  = Src0 * pWhat(iT)
        end do
      end if
c
      deallocate( EVct )
      end
C     *************************************************************************
c
c
c
c                            Time-dependent wave packet
c
c
c
C     *************************************************************************
      subroutine WaveParametersModelGMMA
     i         ( NModes, nzWT, nzW, IdzT, izT,
     i           Omega, OmegaIm, Omega0, kT,
     i           EVctT, AmpT,
     i           DimLessStVct, DoScaling,
     i           t0BC, u0BC, w0BC, qBC,
     i           tnW0, roW0, pW0, uT0, HDensW0,
     o           UThat, URhat, TWhat, pWhat, roWhat )
      implicit none
      integer    NModes, nzWT, nzW, IdzT, qBC
      integer    izT(nzWT+1)
      real*8     Omega, OmegaIm, Omega0, KT
      real*8     t0BC, u0BC, w0BC
      complex*16 EVctT(2*NModes,2*NModes, nzWT)
      complex*16 AmpT(2*NModes,nzWT)
      logical    DimLessStVct, DoScaling
c
      real*8 tnW0(nzW)
      real*8 roW0(nzW)
      real*8 pW0(nzW)
      real*8 uT0(nzW)
      real*8 HDensW0(nzW)
c
c     --- outputs ---
      complex*16 UThat(nzWT)
      complex*16 URhat(nzWT)
      complex*16 TWhat(nzWT)
      complex*16 pWhat(nzWT)
      complex*16 roWhat(nzWT)
c
c     --- local variables ---
      integer     Mrank, iT, i, j, k
      real*8      kX, FctDimLessStVct
      complex*16  im, zero, sum
      complex*16  U, W, T, Wcal
      complex*16  OmegaC, Omega0C, OmegaDop, fct, Src0
      complex*16, allocatable:: EVct(:)
C
      Mrank = 2*NModes
c
      kX  = 1.d-3*kT  !kx in 1/m
C
      OmegaC  = dcmplx( Omega,  OmegaIm ) !Omega
      Omega0C = dcmplx( Omega0, 0.d0 )    !Omega0
c
      im   = (0.d0,1.d0)
      zero = (0.d0,0.d0)
c
c     --- factor for obtaining a dimensionless system of equations ---
      if ( DimLessStVct ) then
        FctDimLessStVct = kX
      else
        FctDimLessStVct = 1.d0
      end if
c     -------------------------------------------------------------------------
c                           Wave Amplitudes
c     -------------------------------------------------------------------------
      allocate( EVct( 2*NModes ) )
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
c
        OmegaDop = OmegaC - dcmplx(kX * uT0(i),0.D0)
        fct = Omega0C / OmegaDop
c       -----------------------------------------------------------------------
c                e(2*NModes) = V(2*NModes,2*NModes) * a(2*NModes)
c       -----------------------------------------------------------------------
        do j = 1, 2*NModes
          sum = zero
          do k = 1, 2*NModes
            sum = sum + EVctT(j,k,iT) * AmpT(k,iT)
          end do
          EVct(j) = sum
        end do
c
        U  = EVct(1)
        W  = EVct(2)
        T  = EVct(3)
        Wcal = FctDimLessStVct * EVct(5)
c
        UThat(iT)  = Omega0C * U / kX !in m/s
        URhat(iT)  = Omega0C * W / kX !in m/s
        TWhat(iT)  = tnW0(i) * T       !in K
c
        roWhat(iT) = roW0(i) * ( fct * U
     &             - im * fct * W / (kX * HDensW0(i))
     &             + im * fct * Wcal / kX )  !in Kg/m**3
        pWhat(iT) = pW0(i) * ( TWhat(iT)/tnW0(i) + roWhat(iT)/roW0(i) ) !in N/m**2
      end do
c
c     --- scale quantities ---
      if ( DoScaling ) then
        if ( qBC .eq. 1 ) then
          Src0 = u0BC / abs( UThat(1) )
        else if ( qBC .eq. 2 ) then
          Src0 = w0BC / abs( URhat(1) )
        else
          Src0 = t0BC / abs( TWhat(1) )
        end if
        do iT = 1, nzWT
          UThat(iT)  = Src0 * UThat(iT)
          URhat(iT)  = Src0 * URhat(iT)
          TWhat(iT)  = Src0 * TWhat(iT)
          roWhat(iT) = Src0 * roWhat(iT)
          pWhat(iT)  = Src0 * pWhat(iT)
        end do
      end if
c
      deallocate( EVct )
      end
c     *************************************************************************
      subroutine TimeFTParameters1
     i         ( NFFT, Omega0, RatioSigmaOmega, NPeriod, TimeMin,
     o           TimeShift, OmegaK, TimeK, SigmaOmega, DOmega )
c     -------------------------------------------------------------------------
c      Example for input: NFFT  = 512, RatioSigmaOmega = 30.0D0
c     -------------------------------------------------------------------------
      implicit none
      integer NFFT, NPeriod
      real*8  Omega0
      real*8  RatioSigmaOmega, TimeMin
c
      real*8  TimeShift, SigmaOmega, DOmega
      real*8  OmegaK( NFFT ), TimeK( NFFT )
c
      integer i, NFFT2
      real*8  pi, OmegaMin, OmegaMax, DeltaOmega
      real*8  SigmaTime, PeriodTime, DTime, DeltaTime
      real*8  TimeMax, KSigmaTime
c
      Pi = acos(-1.d0)
c
C     -- frequency band centered on Omega0 ---
      DOmega = 1.36d-4
      DeltaOmega = DOmega * (NFFT - 1)
      SigmaOmega = min( Omega0 / 3.d0, 2.d0 * DOmega )
      OmegaMin   = Omega0 - 3.d0 * SigmaOmega
      OmegaMax   = Omega0 + 3.d0 * SigmaOmega
      RatioSigmaOmega = Omega0/SigmaOmega
c
C     --- Time interval ---
      SigmaTime  = 1.d0 / SigmaOmega
      PeriodTime = 2.d0*Pi / Omega0
      DeltaTime  = NPeriod * PeriodTime
      DTime      = DeltaTime / ( NFFT - 1 )

      DTime = 2.d0 * Pi / (DOmega*NFFT)
      DeltaTime = DTime * (NFFT - 1)
c
      TimeMin    = TimeMin * 3600.d0 !in seconds
      TimeMax    = TimeMin + DeltaTime
c
      do i = 1, NFFT
        OmegaK(i) = OmegaMin + (i - 1) * DOmega
      end do
      do i = 1, NFFT
        TimeK(i)  = TimeMin + (i - 1) * DTime
      end do
      NFFT2 = NFFT / 2
      TimeShift = (TimeMax - TimeMin) / 2.d0  !TimeK(NFFT2)
      KSigmaTime = DeltaTime / SigmaTime
c
      write(6,*)
      write(6,'(a)')
     & 'Fourier-transform parameters in frequency domain:'
      write(6,'(a,i6)') '- Number of discretization points = ', NFFT
      write(6,'(a,1pe11.4)') '- Reference frequency [1/s] = ', Omega0
      write(6,'(a,1pe11.4)') '- Frequency standard deviation [1/s] = ',
     &      SigmaOmega
      write(6,'(a,1pe11.4)') '- Frequency minimum [1/s]  = ', OmegaMin
      write(6,'(a,1pe11.4)') '- Frequency maximum [1/s]  = ', OmegaMax
      write(6,'(a,1pe11.4)') '- Frequency interval [1/s] = ', DeltaOmega
      write(6,'(a,1pe11.4)') '- Frequency discretization step [1/s] = ',
     &      DOmega
      write(6,'(a)')
     & 'Fourier-transform parameters in time domain:'
      write(6,'(a,1pe11.4)') '- Time period [min] = ',
     &      PeriodTime / 60.0
      write(6,'(a,1pe11.4)') '- Time standard deviation [hr] = ',
     &      SigmaTime / 3600.0
      write(6,'(a,1pe11.4)') '- Time minimum [hr]  = ',
     &      TimeMin / 3600.0
      write(6,'(a,1pe11.4)') '- Time maximum [hr]  = ',
     &      TimeMax / 3600.0
      write(6,'(a,1pe11.4)') '- Time interval [hr] = ',
     &      DeltaTime / 3600.0
      write(6,'(a,1pe11.4)') '- Time discretization step [min] = ',
     &      DTime / 60.0
      write(6,'(a,1pe11.4)')
     &'- Number of time standard deviations per '//
     & 'time interval KSigmaTime ( > 6) = ', KSigmaTime
      if ( KSigmaTime .lt. 6 ) then
        write(6,'(a)')
     &'- Warning: KSigmaTime < 6: increase NPeriod or '//
     & 'reduce RatioSigmaOmega'
      end if
      stop

c
c       write(6,'(a,i6)')
c     &'- Number of time samples in the time period = ',
c     & int( PeriodTime / Dtime)
c
      write(6,*)
C
      end
c     *************************************************************************
      subroutine TimeFTParameters
     i         ( NFFT, Omega0, RatioSigmaOmega, NPeriod, TimeMin,
     o           TimeShift, OmegaK, TimeK, SigmaOmega, DOmega )
c     -------------------------------------------------------------------------
c      Example for input: NFFT  = 512, RatioSigmaOmega = 30.0D0
c     -------------------------------------------------------------------------
      implicit none
      integer NFFT, NPeriod
      real*8  Omega0
      real*8  RatioSigmaOmega, TimeMin
c
      real*8  TimeShift, SigmaOmega, DOmega
      real*8  OmegaK( NFFT ), TimeK( NFFT )
c
      integer i
      real*8  pi, OmegaMin, OmegaMax, DeltaOmega, KSigmaOmega
      real*8  SigmaTime, PeriodTime, DTime, DeltaTime
      real*8  TimeMax, KSigmaTime, Nyquist

      real*8  ks
c
      Pi = acos(-1.d0)
c
C     -- frequency band centered on Omega0 ---
      ks = 3.d0
      SigmaOmega  = Omega0/RatioSigmaOmega
      OmegaMin    = Omega0 - ks * SigmaOmega
      OmegaMax    = Omega0 + ks * SigmaOmega
      DeltaOmega  = OmegaMax - OmegaMin
      DOmega      = DeltaOmega / ( NFFT - 1 )
      KSigmaOmega = DeltaOmega / SigmaOmega
c
C     --- Time interval ---
      SigmaTime  = 1.d0 / SigmaOmega
      PeriodTime = 2.d0*Pi / Omega0
      DeltaTime  = NPeriod * PeriodTime
      DTime      = DeltaTime / ( NFFT - 1 )
c
      TimeMin    = TimeMin * 3600.d0 !in seconds
      TimeMax    = TimeMin + DeltaTime
c
      do i = 1, NFFT
        OmegaK(i) = OmegaMin + (i - 1) * DOmega
      end do
      do i = 1, NFFT
        TimeK(i)  = TimeMin + (i - 1) * DTime
      end do
      TimeShift  = (TimeMax - TimeMin) / 2.d0
      KSigmaTime = DeltaTime / SigmaTime
c
      write(6,*)
      write(6,'(a)')
     & 'Fourier-transform parameters in frequency domain:'
      write(6,'(a,i6)') '- Number of discretization points = ', NFFT
      write(6,'(a,1pe11.4)') '- Reference frequency [1/s] = ', Omega0
      write(6,'(a,1pe11.4)') '- Frequency standard deviation [1/s] = ',
     &      SigmaOmega
      write(6,'(a,1pe11.4)') '- Frequency minimum [1/s]  = ', OmegaMin
      write(6,'(a,1pe11.4)') '- Frequency maximum [1/s]  = ', OmegaMax
      write(6,'(a,1pe11.4)') '- Frequency interval [1/s] = ', DeltaOmega
      write(6,'(a,1pe11.4)') '- Frequency discretization step [1/s] = ',
     &      DOmega
      write(6,'(a,1pe11.4)')
     &'- Number of frequency standard deviations per '//
     & 'frequency interval KSigmaOmega = ', KSigmaOmega
      write(6,'(a)')
     & 'Fourier-transform parameters in time domain:'
      write(6,'(a,1pe11.4)') '- Time period [min] = ',
     &      PeriodTime / 60.0
      write(6,'(a,1pe11.4)') '- Time standard deviation [hr] = ',
     &      SigmaTime / 3600.0
      write(6,'(a,1pe11.4)') '- Time minimum [hr]  = ',
     &      TimeMin / 3600.0
      write(6,'(a,1pe11.4)') '- Time maximum [hr]  = ',
     &      TimeMax / 3600.0
      write(6,'(a,1pe11.4)') '- Time interval [hr] = ',
     &      DeltaTime / 3600.0
      write(6,'(a,1pe11.4)') '- Time discretization step [min] = ',
     &      DTime / 60.0
c
      write(6,'(a,1pe11.4)')
     &'- Number of time standard deviations per '//
     & 'time interval KSigmaTime ( > 6) = ', KSigmaTime
      if ( KSigmaTime .lt. 6 ) then
        write(6,'(a)')
     &'- Warning: KSigmaTime < 6: increase NPeriod or '//
     & 'reduce RatioSigmaOmega'
      end if
      Nyquist = Pi/ OmegaMax
      if ( DTime .lt. Nyquist ) then
        write(6,'(2(a,1pe11.4))')
     &'- Nyquist condition is satisfied: ', DTime, ' < ', Nyquist
      else
        write(6,'(a)')  'Nyquist condition is not satisfied'
      end if
c
      Nyquist = 2.d0*Pi/DeltaTime
      if ( DOmega .lt. Nyquist ) then
        write(6,'(2(a,1pe11.4))')
     &'- Nyquist dual condition is satisfied: ', DOmega, ' < ', Nyquist
      else
        write(6,'(a)')  'Nyquist dual condition is not satisfied'
      end if
      write(6,*)
C
      end
C     *************************************************************************
      subroutine WaveParametersInTime
     i         ( NFFT, TimeShift, nzWT,
     i           Omega0, OmegaIm, RatioSigmaOmega, DOmega,
     i           OmegaK, TimeK, UTFFT, URFFT, TFFT,
     i           DoScaling, t0BC, u0BC, w0BC, qBC,
     o           UTTime, URTime, TTime )
      implicit none
      integer   NFFT, nzWT, qBC
      real*8    TimeShift, Omega0, OmegaIm, RatioSigmaOmega, DOmega
      real*8    OmegaK( NFFT ), TimeK( NFFT )
      real*8    t0BC, u0BC, w0BC
      complex*16 UTFFT(nzWT, NFFT)
      complex*16 URFFT(nzWT, NFFT)
      complex*16 TFFT(nzWT, NFFT)
      logical    DoScaling
c
c     --- outputs ---
      complex*16 UTTime(nzWT, NFFT)
      complex*16 URTime(nzWT, NFFT)
      complex*16 TTime(nzWT, NFFT)
!
!     --- local variables ---
      integer    iT, i, k, NFFT2
      real*8     pi, SigmaOmega, Time, Omega, fct, Src0
      complex*16 im, Z, zero, sum
      complex*16 SrcFrequencyAmplit
c
      im   = (0.d0,1.d0)
      zero = (0.d0,0.d0)
      pi   = acos(-1.d0)
c
c     --- omega standard deviation ---
      SigmaOmega = Omega0 / RatioSigmaOmega
c
      do i = 1, NFFT
        Time = TimeK(i)
        fct  = DOmega / (2.d0*pi)
        do iT = 1, nzWT
C
C         --- horizontal velocity ---
          sum  = zero
          do k = 1, NFFT
            Omega = OmegaK(k)
            Z = fct * UTFFT(iT,k) * exp( im*Omega*(Time-TimeShift) )
     &        * SrcFrequencyAmplit(Omega,OmegaIm,Omega0,SigmaOmega)
            sum = sum + Z
          end do
          UTTime(iT,i) = sum
          UTTime(iT,i) = exp(-OmegaIm*(Time-TimeShift)) * UTTime(iT,i)
C
C         --- vertical velocity ---
          sum  = zero
          do k = 1, NFFT
            Omega = OmegaK(k)
            Z = fct * URFFT(iT,k) * exp( im*Omega*(Time-TimeShift) )
     &        * SrcFrequencyAmplit(Omega,OmegaIm,Omega0,SigmaOmega)
            sum = sum + Z
          end do
          URTime(iT,i) = sum
          URTime(iT,i) = exp(-OmegaIm*(Time-TimeShift)) * URTime(iT,i)
c
c         --- temperature ---
          sum  = zero
          do k = 1, NFFT
            Omega = OmegaK(k)
            Z = fct * TFFT(iT,k) * exp( im*Omega*(Time-TimeShift) )
     &        * SrcFrequencyAmplit(Omega,OmegaIm,Omega0,SigmaOmega)
            sum = sum + Z
          end do
          TTime(iT,i) = sum
          TTime(iT,i) = exp(-OmegaIm*(Time-TimeShift)) * TTime(iT,i)
        end do
      end do
c
c     --- scale quantities ---
      if ( DoScaling ) then
        NFFT2 = NFFT / 2
        if ( qBC .eq. 1 ) then
          Src0 = u0BC / abs( dble(UTTime(1,NFFT2)) )
        else if ( qBC .eq. 2 ) then
          Src0 = w0BC / abs( dble(URTime(1,NFFT2)) )
        else
          Src0 = t0BC / abs( dble(TTime(1,NFFT2)) )
        end if
        do i = 1, NFFT
          do iT = 1, nzWT
            UTTime(iT,i) = Src0 * UTTime(iT,i)
            URTime(iT,i) = Src0 * URTime(iT,i)
            TTime(iT,i)  = Src0 * TTime(iT,i)
          end do
        end do
      end if
      end
c     *************************************************************************
      subroutine MaxWaveParametersInTime
     i         ( NFFT, nzWT, izT, zmin, dz,
     i           TimeK, UTTime, URTime, TTime,
     o           UMax, timeUMax, zUMax,
     o           WMax, timeWMax, zWMax,
     o           TMax, timeTMax, zTMax )
      implicit none
      integer NFFT, nzWT
      integer izT(nzWT+1)
      real*8  zmin, dz
      real*8  TimeK( NFFT )
c
      complex*16 UTTime(nzWT, NFFT)
      complex*16 URTime(nzWT, NFFT)
      complex*16 TTime(nzWT, NFFT)
c
      real*8 UMax, timeUMax, zUMax
      real*8 WMax, timeWMax, zWMax
      real*8 TMax, timeTMax, zTMax
C
C     --- locals ---
      integer  j, iT, i
      real*8   time, z
c
      UMax = -1.d+20
      do j = 1, NFFT
        time = TimeK(j)/60.0 !in min
        do iT = 2, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz !in km
          if ( dble(UTTime(iT, j)) .gt. UMax ) then
            UMax  = dble(UTTime(iT, j))
            timeUMax = time !in min
            zUMax = z !in km
          end if
        end do
      end do
c
      WMax = -1.d+20
      do j = 1, NFFT
        time = TimeK(j)/60.0 !in min
        do iT = 2, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz !in km
          if ( dble(URTime(iT, j)) .gt. WMax ) then
            WMax  = dble(URTime(iT, j))
            timeWMax = time !in min
            zWMax = z !in km
          end if
        end do
      end do
C
      TMax = -1.d+20
      do j = 1, NFFT
        time = TimeK(j)/60.0 !in min
        do iT = 2, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz !in km
          if ( dble(TTime(iT, j)) .gt. TMax ) then
            TMax  = dble(TTime(iT, j))
            timeTMax = time !in min
            zTMax = z !in km
          end if
        end do
      end do
      end
c     *************************************************************************
      subroutine SelectSolution
     i         ( nzWT, NFFT, NOmegaIm, NSol,
     i           UMaxSol, WMaxSol, TMaxSol, OmegaImSol,
     i           UTTimeSol, URTimeSol, TTimeSol,
     o           UMax, WMax, TMax, OmegaIm,
     o           UTTime, URTime, TTime )
      implicit none
      integer    nzWT, NFFT, NOmegaIm, NSol
      real*8     UMaxSol(NOmegaIm), WMaxSol(NOmegaIm)
      real*8     TMaxSol(NOmegaIm), OmegaImSol(NOmegaIm)
      complex*16 UTTimeSol(nzWT, NFFT, NOmegaIm)
      complex*16 URTimeSol(nzWT, NFFT, NOmegaIm)
      complex*16 TTimeSol(nzWT, NFFT, NOmegaIm)
c
      real*8     UMax, WMax, TMax, OmegaIm
      complex*16 UTTime(nzWT, NFFT)
      complex*16 URTime(nzWT, NFFT)
      complex*16 TTime(nzWT, NFFT)
c
c     --- local variables ---
      integer  isol, iT, j
      real*8   UC, WC, TC, distMin, dist
c
      if ( NSol .eq. 1 ) then
        isol = 1
        UMax = UMaxSol(isol)
        WMax = WMaxSol(isol)
        TMax = TMaxSol(isol)
        OmegaIm = OmegaImSol(isol)
        do j = 1, NFFT
          do iT = 2, nzWT
            UTTime(iT,j) = UTTimeSol(iT,j,isol)
            URTime(iT,j) = URTimeSol(iT,j,isol)
            TTime(iT,j)  = TTimeSol(iT,j,isol)
          end do
        end do
      else if ( NSol .eq. 2 ) then
        if ( abs(OmegaImSol(1)) .lt. abs(OmegaImSol(2)) ) then
          isol = 1
        else
          isol = 2
        end if
        UMax = UMaxSol(isol)
        WMax = WMaxSol(isol)
        TMax = TMaxSol(isol)
        OmegaIm = OmegaImSol(isol)
        do j = 1, NFFT
          do iT = 2, nzWT
            UTTime(iT,j) = UTTimeSol(iT,j,isol)
            URTime(iT,j) = URTimeSol(iT,j,isol)
            TTime(iT,j)  = TTimeSol(iT,j,isol)
          end do
        end do
      else
c
c       --- center of mass in (Umax,WMax,TMax) space ---
        UC = 0.d0
        WC = 0.d0
        TC = 0.d0
        do isol = 1, NSol
          UC = UC + UMaxSol(isol)
          WC = WC + WMaxSol(isol)
          TC = TC + TMaxSol(isol)
        end do
        UC = UC / dble(NSol)
        WC = WC / dble(NSol)
        TC = TC / dble(NSol)
c
        distMin = 1.d+20
        do isol = 1, NSol
          dist = (UMaxSol(isol) - UC)**2
     &         + (WMaxSol(isol) - WC)**2
     &         + (TMaxSol(isol) - TC)**2
          if ( dist .lt. distMin ) then
            distMin = dist
            UMax = UMaxSol(isol)
            WMax = WMaxSol(isol)
            TMax = TMaxSol(isol)
            OmegaIm = OmegaImSol(isol)
            do j = 1, NFFT
              do iT = 2, nzWT
                UTTime(iT,j) = UTTimeSol(iT,j,isol)
                URTime(iT,j) = URTimeSol(iT,j,isol)
                TTime(iT,j)  = TTimeSol(iT,j,isol)
              end do
            end do
          end if
        end do
      end if
      end
c     *************************************************************************
c
c
c                             Source function
c
c
c
c     *************************************************************************
      subroutine SourceFunctionReconstruction
     i         ( NFFT, TimeShift,
     i           Omega0, OmegaIm, RatioSigmaOmega, DOmega,
     i           OmegaK, TimeK,
     i           rms )
      implicit none
      integer NFFT
      real*8  Omega0, OmegaIm
      real*8  TimeShift, RatioSigmaOmega, DOmega
      real*8  OmegaK( NFFT ), TimeK( NFFT )
c
      real*8  rms
c
c     --- local variables ---
      integer    i, k
      real*8     pi, SigmaOmega, Time, Omega, fct
      real*8     SrcMax1, SrcMax2,  eps, norm
      complex*16 im, Z, zero, sum
      complex*16 SrcFrequencyAmplit, SrcTimeAmplit
      character*72 path, filename
      complex*16, allocatable:: Src1(:), Src2(:)
C
      path  = 'WaveResultsWP/'
c
      im   = (0.d0,1.d0)
      zero = (0.d0,0.d0)
      pi   = acos(-1.d0)
c
c     --- omega standard deviation ---
      SigmaOmega = Omega0 / RatioSigmaOmega

      allocate( Src1(NFFT), Src2(NFFT) )
      filename = TRIM(path)//'TimeSrcFctReconstruct.dat'
      open(unit=44,file=trim(filename),status='unknown')
      eps  = 0.d0
      norm = 0.d0
      fct  = DOmega / (2.0*pi)
      SrcMax1 = -1.d+20
      SrcMax2 = -1.d+20
      do i = 1, NFFT
        Time = TimeK(i)
        sum  = zero
        do k = 1, NFFT
          Omega = OmegaK(k)
          Z = fct * exp( im*Omega*(Time-TimeShift) )
     &      * SrcFrequencyAmplit(Omega,OmegaIm,Omega0,SigmaOmega)
          sum = sum + Z
        end do
        Src1(i) = sum * exp(-OmegaIm*(Time-TimeShift))
        Src2(i) = SrcTimeAmplit( Time, TimeShift,
     &            Omega0, SigmaOmega )
c
        SrcMax1 = max( dble(Src1(i)), SrcMax1 )
        SrcMax2 = max( dble(Src2(i)), SrcMax2 )
      end do
      do i = 1, NFFT
        Src1(i) = Src1(i) / SrcMax1
        Src2(i) = Src2(i) / SrcMax2
c
        eps = eps + abs( Src1(i) - Src2(i) )**2
        norm = norm + abs(Src2(i))**2
        write(44,'(3(2x,1pe11.4))')
     &        Time/3600.0, dble(Src1(i)), dble(Src2(i))
      end do
      rms = sqrt(eps/norm)
      close(unit=44)
c
      deallocate( Src1, Src2 )
      end
c     *************************************************************************
      complex*16 function SrcFrequencyAmplit
     &         ( Omega, OmegaIm, Omega0, SigmaOmega )
c     -------------------------------------------------------------------------
c        Source function excluding the factor
c                  exp( 0.5d0 * OmegaIm**2 / SigmaOmega**2)
c     -------------------------------------------------------------------------
      implicit none
      real*8   Omega, OmegaIm, Omega0, SigmaOmega
C
C     --- Locals ---
      real*8      pi, arg
      complex*16  im, zero, argC
C
      pi   = acos(-1.D0)
      im   = dcmplx(0.d0, 1.d0)
      zero = dcmplx(0.d0, 0.d0)
c
C     --- arguments of exponentials ---
      arg  = ( (Omega - Omega0) / SigmaOmega )**2
      argC = im * OmegaIm * (Omega - Omega0) / SigmaOmega**2
c
      if ( arg .lt. 100.0D0 ) then
         SrcFrequencyAmpliT = sqrt(2.0D0*pi)
     &       * exp( -0.5d0*arg ) * exp(-argC) / SigmaOmega
      else
         SrcFrequencyAmpliT = zero
      end if
c
c      argC = ( dcmplx( Omega - Omega0, OmegaIm) / SigmaOmega )**2
c      if ( arg .lt. 100.0D0 ) then
c         SrcFrequencyAmpliT = sqrt(2.0D0*pi) * exp( -0.5d0*argC )
c     &                      / SigmaOmega
c      else
c         SrcFrequencyAmpliT = zero
c      end if
c
      end
c     *************************************************************************
      complex*16 function SrcTimeAmplit( Time, TimeShift,
     &           Omega0, SigmaOmega )
      implicit none
      real*8  Time, TimeShift
      real*8  Omega0, SigmaOmega
C
C     --- Locals ---
      real*8      arg, arg2, gauss, phase
      complex*16  zero, phase_factor
c
      zero = dcmplx(0.0D0, 0.0D0)

      arg  = (Time - TimeShift) * SigmaOmega
      arg2 = arg * arg
c
C     --- Underflow/overflow safety for the Gaussian ---
      if ( 0.5D0*arg2 .GT. 700.0D0 ) then
         SrcTimeAmpliT = zero
         return
      end if
      gauss = EXP( -0.5D0 * arg2 )
c
C     --- Pure phase term: exp(i*Omega0*(Time-TimeShift)) = cos + i sin ---
      phase = Omega0 * (Time - TimeShift)
      phase_factor = DCMPLX( DCOS(phase), DSIN(phase) )

      SrcTimeAmpliT = gauss * phase_factor
      end
c     *************************************************************************
c
c
c                         RMS errors of solution methods
c
c
c
c     *************************************************************************
      subroutine RMSErrors
     I         ( TypeLinModel, nzWT,
     o           ErrU12, ErrW12, ErrT12,
     o           ErrU13, ErrW13, ErrT13,
     o           ErrU23, ErrW23, ErrT23)
      implicit none
      integer  TypeLinModel, nzWT
C
c     --- outputs ---
      real*8  ErrU12, ErrW12, ErrT12
      real*8  ErrU13, ErrW13, ErrT13
      real*8  ErrU23, ErrW23, ErrT23
c
c     --- locals ---
      integer  iT
      real*8   z, sumU, sumW, sumT
      real*8   NormU1, NormU2, NormU3
      real*8   NormW1, NormW2, NormW3
      real*8   NormT1, NormT2, NormT3
      character*1  charM
      character*72 path, path0, filename
      real*8, allocatable:: U1(:), U2(:), U3(:)
      real*8, allocatable:: W1(:), W2(:), W3(:)
      real*8, allocatable:: T1(:), T2(:), T3(:)
c
      write(charM,'(I1)') TypeLinModel
C
      allocate( U1(nzWT), U2(nzWT), U3(nzWT) )
      allocate( W1(nzWT), W2(nzWT), W3(nzWT) )
      allocate( T1(nzWT), T2(nzWT), T3(nzWT) )
c
      U1 = 0.d0
      U2 = 0.d0
      U3 = 0.d0
      W1 = 0.d0
      W2 = 0.d0
      W3 = 0.d0
      T1 = 0.d0
      T2 = 0.d0
      T3 = 0.d0
C
      path  = 'WaveResultsSF/'
c     -------------------------------------------------------------------------
c                         Solution Method GMMA
c     -------------------------------------------------------------------------
      path0 = 'SolMetGMMALinModel'//TRIM(charM)//'.dat'
c
c     --- U1 ---
      filename = TRIM(path)//'UVelocity'//TRIM(path0)
      call CheckIfFileExists( filename )
      open(unit=13,file=trim(filename),status='unknown')
      read(13,*)
      do iT = 2, nzWT
        read(13,*) U1(iT), z
      end do
      close(13)
c
c     --- W1 ---
      filename = TRIM(path)//'WVelocity'//TRIM(path0)
      call CheckIfFileExists( filename )
      open(unit=13,file=trim(filename),status='unknown')
      read(13,*)
      do iT = 2, nzWT
        read(13,*) W1(iT), z
      end do
      close(13)
c
c     --- T1 ---
      filename = TRIM(path)//'Temperature'//TRIM(path0)
      call CheckIfFileExists( filename )
      open(unit=13,file=trim(filename),status='unknown')
      read(13,*)
      do iT = 2, nzWT
        read(13,*) T1(iT), z
      end do
      close(13)
c     -------------------------------------------------------------------------
c                         Solution Method GMMN
c     -------------------------------------------------------------------------
      path0 = 'SolMetGMMNLinModel'//TRIM(charM)//'.dat'
c
c     --- U2 ---
      filename = TRIM(path)//'UVelocity'//TRIM(path0)
      call CheckIfFileExists( filename )
      open(unit=13,file=trim(filename),status='unknown')
      read(13,*)
      do iT = 2, nzWT
        read(13,*) U2(iT), z
      end do
      close(13)
c
c     --- W2 ---
      filename = TRIM(path)//'WVelocity'//TRIM(path0)
      call CheckIfFileExists( filename )
      open(unit=13,file=trim(filename),status='unknown')
      read(13,*)
      do iT = 2, nzWT
        read(13,*) W2(iT), z
      end do
      close(13)
c
c     --- T2 ---
      filename = TRIM(path)//'Temperature'//TRIM(path0)
      call CheckIfFileExists( filename )
      open(unit=13,file=trim(filename),status='unknown')
      read(13,*)
      do iT = 2, nzWT
        read(13,*) T2(iT), z
      end do
      close(13)
c     -------------------------------------------------------------------------
c                         Solution Method SMMA
c     -------------------------------------------------------------------------
      path0 = 'SolMetSMMALinModel'//TRIM(charM)//'.dat'
c
c     --- U3 ---
      filename = TRIM(path)//'UVelocity'//TRIM(path0)
      call CheckIfFileExists( filename )
      open(unit=13,file=trim(filename),status='unknown')
      read(13,*)
      do iT = 2, nzWT
        read(13,*) U3(iT), z
      end do
      close(13)
c
c     --- W3 ---
      filename = TRIM(path)//'WVelocity'//TRIM(path0)
      call CheckIfFileExists( filename )
      open(unit=13,file=trim(filename),status='unknown')
      read(13,*)
      do iT = 2, nzWT
        read(13,*) W3(iT), z
      end do
      close(13)
c
c     --- T3 ---
      filename = TRIM(path)//'Temperature'//TRIM(path0)
      call CheckIfFileExists( filename )
      open(unit=13,file=trim(filename),status='unknown')
      read(13,*)
      do iT = 2, nzWT
        read(13,*) T3(iT), z
      end do
      close(13)
c     -------------------------------------------------------------------------
c                               Norms
c     -------------------------------------------------------------------------
      sumU = 0.d0
      sumW = 0.d0
      sumT = 0.d0
      do iT = 2, nzWT
        sumU = sumU + U1(iT)**2
        sumW = sumW + W1(iT)**2
        sumT = sumT + T1(iT)**2
      end do
      NormU1 = sqrt(sumU)
      NormW1 = sqrt(sumW)
      NormT1 = sqrt(sumT)
C
      sumU = 0.d0
      sumW = 0.d0
      sumT = 0.d0
      do iT = 2, nzWT
        sumU = sumU + U2(iT)**2
        sumW = sumW + W2(iT)**2
        sumT = sumT + T2(iT)**2
      end do
      NormU2 = sqrt(sumU)
      NormW2 = sqrt(sumW)
      NormT2 = sqrt(sumT)
C
      sumU = 0.d0
      sumW = 0.d0
      sumT = 0.d0
      do iT = 2, nzWT
        sumU = sumU + U3(iT)**2
        sumW = sumW + W3(iT)**2
        sumT = sumT + T3(iT)**2
      end do
      NormU3 = sqrt(sumU)
      NormW3 = sqrt(sumW)
      NormT3 = sqrt(sumT)
c     -------------------------------------------------------------------------
c                             Errors 1-2
c     -------------------------------------------------------------------------
      sumU = 0.d0
      sumW = 0.d0
      sumT = 0.d0
      do iT = 2, nzWT
        sumU = sumU + ( U1(iT) - U2(iT) )**2
        sumW = sumW + ( W1(iT) - W2(iT) )**2
        sumT = sumT + ( T1(iT) - T2(iT) )**2
      end do
      ErrU12 = sqrt(sumU) / NormU1
      ErrW12 = sqrt(sumW) / NormW1
      ErrT12 = sqrt(sumT) / NormT1
c     -------------------------------------------------------------------------
c                             Errors 1-3
c     -------------------------------------------------------------------------
      sumU = 0.d0
      sumW = 0.d0
      sumT = 0.d0
      do iT = 2, nzWT
        sumU = sumU + ( U1(iT) - U3(iT) )**2
        sumW = sumW + ( W1(iT) - W3(iT) )**2
        sumT = sumT + ( T1(iT) - T3(iT) )**2
      end do
      ErrU13 = sqrt(sumU) / NormU1
      ErrW13 = sqrt(sumW) / NormW1
      ErrT13 = sqrt(sumT) / NormT1
c     -------------------------------------------------------------------------
c                             Errors 2-3
c     -------------------------------------------------------------------------
      sumU = 0.d0
      sumW = 0.d0
      sumT = 0.d0
      do iT = 2, nzWT
        sumU = sumU + ( U3(iT) - U2(iT) )**2
        sumW = sumW + ( W3(iT) - W2(iT) )**2
        sumT = sumT + ( T3(iT) - T2(iT) )**2
      end do
      ErrU23 = sqrt(sumU) / NormU2
      ErrW23 = sqrt(sumW) / NormW2
      ErrT23 = sqrt(sumT) / NormT2
C
      deallocate( U1, U2, U3 )
      deallocate( W1, W2, W3 )
      deallocate( T1, T2, T3 )
      end
c     ***********************************************************************
      subroutine CheckIfFileExists
     I           ( FileName )
      implicit none
      character*72 FileName

      logical FileExist

      INQUIRE(FILE=TRIM(FileName), EXIST= FileExist )
      if ( .not. FileExist ) then
          write(6,'(A)')
     &   'The file '// TRIM(FileName)//' does not exist'
          stop
      end if
      end
C     *************************************************************************
C
C
C
c
C                             Spatial Fourier Transform
C
C
C
C     *************************************************************************
      subroutine SpatialFourierTransform
     i         ( TypeSolMet, TypeLinModel,
     i           nzWT, IdzT, izT, dz,
     i           UThat, URhat, TWhat )
      implicit none
      integer TypeSolMet, TypeLinModel, nzWT, IdzT
      integer izT(nzWT+1)
      real*8  dz
C
      complex*16 UThat(nzWT)
      complex*16 URhat(nzWT)
      complex*16 TWhat(nzWT)
c
c     --- local variables ---
      integer  nfft, iT, i, IMin, IMax
      real*8   pi, LregSFT,  Omegazmin, OmegaZMax
      real*8   Lambdazmin, LambdaZMax
      real*8   dzF, DOmegaZ, ErrRMS, ErrMax, AMax
      logical  DoSFTSmoothing
      character*1   charM
      character*72  path, path0, filename
      real*8, allocatable:: zWT( : ), OmegaZ( : ), Fct( : )
      real*8, allocatable:: AK( : ), AKS( : )
c
c     --- parameters ---
      pi   = acos(-1.d0)
      nfft = 4096
      DoSFTSmoothing = .false.
      LregSFT    = 1.d+2
      Lambdazmin = 50.d0
c
      write(charM,'(I1)') TypeLinModel
c
      path  = 'FourierTransformSF/'
      if (TypeSolMet .eq. 1 ) then
        path0 = 'SolMetGMMALinModel'//TRIM(charM)//'.dat'
      else if (TypeSolMet .eq. 2 ) then
        path0 = 'SolMetGMMNLinModel'//TRIM(charM)//'.dat'
      else
        path0 = 'SolMetSMMALinModel'//TRIM(charM)//'.dat'
      end if
c
      allocate( zWT(nzWT) )
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        zWT(iT) = (i-1)*dz
      end do
      dzF = zWT(2) - zWT(1)
c
      DOmegaZ = 2.d0 * pi / (nfft * dzF)
c
      allocate( OmegaZ( nfft ), Fct( nfft ) )
      allocate( AK( nfft ), AKS( nfft ) )
      do i = 1, nfft
        OmegaZ(i) = (i - 1) * DOmegaZ
      end do
      Omegazmin  = 10.d0 * DOmegaZ
      LambdaZMax = 2.d0 * pi / Omegazmin
c
      OmegaZMax  = 2.d0 * pi / Lambdazmin
      IMax = int( OmegaZMax / DOmegaZ )
c
      IMin = 3 !cut low frequency
c     ...................................................................
c                          FT of Temperature
c     ...................................................................
      do iT = 1, nzWT
        Fct(iT) = dble(TWhat(iT))
      end do
      call DiscreteCosineTransform
     i   ( nzWT, nfft, Fct,
     o     AK, ErrRMS, ErrMax )
      do i = 1, nfft
        AK(i) = abs(AK(i))
      end do
c
      AMax = -1.d+20
      do I = IMin, nfft
        if ( AK(I) .gt. AMax ) then
          AMax = AK(I)
          OmegaZMax = OmegaZ(I)
        end if
      end do
c
c     --- smooth spectrum for printing ---
      if ( DoSFTSmoothing ) then
          call csplineD
     i      ( .true., nfft, AK, LregSFT,
     o         AKS )
          do i = IMin, nfft
            AK(i) = AKS(i)
          end do
      end if
c
      filename = TRIM(path)//'Temperature'//TRIM(path0)
      open(unit=31,file=TRIM(filename),status='UNKNOWN')
      write(31,'(a)') '# Spatial Fourier Transform of Temperature'
      write(31,'(a,1pe13.4)')
     & '# Maximum wavelength LambdaZ = ', 2.0 * pi / OmegaZMax
      write(31,'(a,1pe13.4)')
     & '# Maximum wavenumber kZ  = ', OmegaZMax

      do i = IMin, IMax
        write(31,'(1PE13.4,10(2X,1PE13.4))')
     &  OmegaZ(I), AK(I) / AMax
      end do
      close(unit=31)
c
      write(6,*)
      write(6,*) 'Spatial Fourier Transform'
      write(6,'(a,1pe13.4)')
     & '- Minimum wavelength = ', Lambdazmin
      write(6,'(a,1pe13.4)')
     & '- Maximum Wavelength = ', LambdaZMax
c
      write(6,'(2(a,1pe13.4))')
     & '- Wavelength for maximum temperature          = ',
     & 2.0 * pi / OmegaZMax,
     & ', RMS error in Fourier transform = ', ErrRMS
c     ...................................................................
c                          FT of vertical velocity
c     ...................................................................
      do iT = 1, nzWT
        Fct(iT) = dble(URhat(iT))
      end do
      call DiscreteCosineTransform
     i   ( nzWT, nfft, Fct,
     o     AK, ErrRMS, ErrMax )
      do i = 1, nfft
        AK(i) = abs(AK(i))
      end do
c
      AMax = -1.d+20
      do I = IMin, nfft
        if ( AK(I) .gt. AMax ) then
          AMax = AK(I)
          OmegaZMax = OmegaZ(I)
        end if
      end do
c
c     --- smooth spectrum for printing ---
      if ( DoSFTSmoothing ) then
          call csplineD
     i      ( .true., nfft, AK, LregSFT,
     o         AKS )
          do i = IMin, nfft
            AK(i) = AKS(i)
          end do
      end if
c
      filename = TRIM(path)//'WVelocity'//TRIM(path0)
      open(unit=31,file=TRIM(filename),status='UNKNOWN')
      write(31,'(a)') '# Spatial Fourier Transform of Vertical Velocity'
      write(31,'(a,1pe13.4)')
     & '# Maximum wavelength LambdaZ = ', 2.0 * pi / OmegaZMax
      write(31,'(a,1pe13.4)')
     & '# Maximum wavenumber kZ  = ', OmegaZMax
      do i = IMin, IMax
        write(31,'(1PE13.4,10(2X,1PE13.4))')
     &  OmegaZ(I), AK(I) / AMax
      end do
      close(unit=31)
c
      write(6,'(2(a,1pe13.4))')
     & '- Wavelength for maximum vertical velocity    = ',
     &  2.0 * pi / OmegaZMax,
     & ', RMS error in Fourier transform = ', ErrRMS
c     ...................................................................
c                         FT of horizontal velocity
c     ...................................................................
      do iT = 1, nzWT
        Fct(iT) = dble(UThat(iT))
      end do
      call DiscreteCosineTransform
     i   ( nzWT, nfft, Fct,
     o     AK, ErrRMS, ErrMax )
      do i = 1, nfft
        AK(i) = abs(AK(i))
      end do
c
      AMax = -1.d+20
      do I = IMin, nfft
        if ( AK(I) .gt. AMax ) then
          AMax = AK(I)
          OmegaZMax = OmegaZ(I)
        end if
      end do
c
c     --- smooth spectrum for printing ---
      if ( DoSFTSmoothing ) then
          call csplineD
     i      ( .true., nfft, AK, LregSFT,
     o         AKS )
          do i = IMin, nfft
            AK(i) = AKS(i)
          end do
      end if
c
      filename = TRIM(path)//'UVelocity'//TRIM(path0)
      open(unit=31,file=TRIM(filename),status='UNKNOWN')
      write(31,'(a)')
     &     '# Spatial Fourier Transform of Horizontal Velocity'
      write(31,'(a,1pe13.4)')
     & '# Maximum wavelength LambdaZ = ', 2.0 * pi / OmegaZMax
      write(31,'(a,1pe13.4)')
     & '# Maximum wavenumber kZ  = ', OmegaZMax
      do i = IMin, IMax
        write(31,'(1PE13.4,10(2X,1PE13.4))')
     &  OmegaZ(I), AK(I) / AMax
      end do
      close(unit=31)
c
      write(6,'(2(a,1pe13.4))')
     & '- Wavelength for maximum horizontal velocity  = ',
     &  2.0 * pi / OmegaZMax,
     & ', RMS error in Fourier transform = ', ErrRMS
c
      deallocate( zWT, OmegaZ, Fct, AK, AKS )
      end
c     *************************************************************************
      subroutine DiscreteCosineTransform
     i         ( nzWT, nfft, Fct,
     o           Ck, err_rms, err_max )
      implicit none
c
c     --- inputs ---
      integer nzWT, nfft
      real*8  Fct(nzWT)
c
c     --- outputs ---
      real*8  Ck(nfft)
      real*8  err_rms, err_max
c
c     --- locals ---
      integer i, j, iplen
      real*8  sum2, diff
c
c     --- allocatable work arrays ---
      integer, allocatable :: ip(:)
      real*8,  allocatable :: w(:)
      real*8,  allocatable :: a(:)
c
c     --- checks ---
      if (nzWT .lt. 2) then
         write(6,*) 'DCT_ALTITUDE_DDCT: nzWT must be >= 2'
         stop
      end if

      if (nfft .lt. nzWT) then
         write(6,*) 'DCT_ALTITUDE_DDCT: nfft must be >= nzWT'
         stop
      end if
c
c     --- check nfft power of 2 ---
      j = 1
 10   continue
      if (j .lt. nfft) then
         j = j*2
         goto 10
      end if
      if (j .ne. nfft) then
         write(6,*) 'DCT_ALTITUDE_DDCT: nfft must be a power of 2'
         stop
      end if
c
c     --- allocate work arrays ---
c     ip length: DDCT needs >= 2+sqrt(n/2). Use safe upper bound 2+nfft/4.
      iplen = 2 + nfft/4
      allocate(ip(0:iplen-1))
      allocate(w(0:(5*nfft)/4 - 1))
      allocate(a(0:nfft-1))
c
c     --- pack data and zero-pad ---
      do i = 0, nfft-1
         a(i) = 0.0d0
      end do
      do i = 1, nzWT
         a(i-1) = Fct(i)
      end do
c
c     --- forward DCT: isgn = -1 ---
      ip(0) = 0
      call ddct(nfft, -1, a, ip, w)
c
c     --- copy coefficients ---
      do i = 1, nfft
         Ck(i) = a(i-1)
      end do
c     -------------------------------------------------------------------------
c     Inverse check:
c              a(0)=a(0)/2
c              call ddct(n, 1, a, ip, w)
c              a(j)=a(j)*2/n
c     -------------------------------------------------------------------------
      a(0) = a(0) / 2.0d0
      call ddct(nfft, 1, a, ip, w)
      do i = 0, nfft-1
         a(i) = a(i) * 2.0d0 / dble(nfft)
      end do
c
c     --- error measures on the original nzWT points ---
      sum2    = 0.0d0
      err_max = 0.0d0
      do i = 1, nzWT
         diff = a(i-1) - Fct(i)
         sum2 = sum2 + diff*diff
         if (dabs(diff) .gt. err_max) err_max = dabs(diff)
      end do
      err_rms = dsqrt( sum2 / dble(nzWT) )
c
      deallocate(a)
      deallocate(w)
      deallocate(ip)
      return
      end

