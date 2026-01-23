C     *************************************************************************
C
C
C
C                               Printing routines
C
C
C
C     *************************************************************************
      subroutine PrintAtmosphericParameters
     i          ( nzW, nzWT, IdzT, izT, zmin, zminDiffVelocity, dz,
     i            tnW0, roW0, pW0, uT0, HScaleW0,
     i            HDensW0, CvW0, GammaW0, CsoundW0,
     i            deniW0S, nu0niW, nu0inW, uD0,
     i            DroW0, DpW0, DtnW0, DuT0, DHDensW0,
     i            DdeniW0, DuD0 )
      implicit none
      integer nzW, nzWT, IdzT
      integer izT(nzWT+1)
      real*8  zmin, zminDiffVelocity, dz
c
      real*8 tnW0(nzW)
      real*8 roW0(nzW)
      real*8 pW0(nzW)
      real*8 uT0(nzW)
      real*8 HScaleW0(nzW)
      real*8 HDensW0(nzW)
      real*8 CvW0(nzW)
      real*8 GammaW0(nzW)
      real*8 CsoundW0(nzW)
c
      real*8 DroW0(nzW)
      real*8 DpW0(nzW)
      real*8 DtnW0(nzW)
      real*8 DuT0(nzW)
      real*8 DHDensW0(nzW)
c
      real*8 deniW0S(nzW)
      real*8 nu0niW(nzW)
      real*8 nu0inW(nzW)
      real*8 uD0(nzW)
      real*8 DdeniW0(nzW)
      real*8 DuD0(nzW)
c
c     --- local ---
      integer iT, i
      real*8  z
      character*72 filename
c
      filename = 'AtmosPars/Temperature.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Temperature versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    tnW0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/MassDensity.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Mass density versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    roW0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/Pressure.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Pressure versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    pW0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/UVelocity.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Horizontal velocity versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    uT0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/ScaleHeightAtmos.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Scale height versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    HScaleW0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/ScaleHeightDensity.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Density scale height versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    HDensW0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/HeatConstantVolume.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Specific heat constant versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    CvW0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/RatioSpecificHeats.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Gamma versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    GammaW0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/SoundSpeed.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Sound speed versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    CsoundW0(i), z
      end do
      close(unit=44)
c     -------------------------------------------------------------------------
c                            Derivatives
c     -------------------------------------------------------------------------
      filename = 'AtmosPars/TemperatureDerivative.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Temperature derivative versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    DtnW0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/UVelocityDerivative.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Horizontal-velocity derivative versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    DuT0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/MassDensityDerivative.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Mass-density derivative versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    DroW0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/PressureDerivative.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Pressure derivative versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    DpW0(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/HScaleDensDerivative.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)')
     & '# Derivative of density-scale height versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    DHDensW0(i), z
      end do
      close(unit=44)
c     -------------------------------------------------------------------------
c                          Ion parameters
c     -------------------------------------------------------------------------
      filename = 'AtmosPars/IonNumberDensity.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Ion number density versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    deniW0S(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/FrequencyNeutralIon.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)')
     &  '# Neutral-ion collision frequency versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    nu0niW(i), z
      end do
      close(unit=44)
c
      filename = 'AtmosPars/FrequencyIonNeutral.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)')
     &  '# Neutral-ion collision frequency versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          if ( z .gt. 100.d0 ) then
            write(44,'(1pe13.4,2x,1pe13.4)')
     &      nu0inW(i), z
          end if
      end do
      close(unit=44)
c
      filename = 'AtmosPars/DiffusionVelocity.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)') '# Diffusion velocity versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          if ( z .gt. zminDiffVelocity ) then
            write(44,'(1pe13.4,2x,1pe13.4)')
     &      uD0(i), z
          end if
      end do
      close(unit=44)
c
      filename = 'AtmosPars/DiffusionVelocityDerivative.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)')
     &'# Diffusion velocity derivative versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          if ( z .gt. zminDiffVelocity ) then
            write(44,'(1pe13.4,2x,1pe13.4)')
     &      DuD0(i), z
          end if
      end do
      close(unit=44)
c
      filename = 'AtmosPars/IonNumberDensityDerivative.dat'
      open(unit=44,file=trim(filename),status='unknown')
      write(44,'(a)')
     & '# Derivative of ion number density versus altitude'
      do iT = 1, nzWT
          i = izT(iT) + IdzT / 2
          z = zmin + (i-1)*dz
          write(44,'(1pe13.4,2x,1pe13.4)')
     &    DdeniW0(i), z
      end do
      close(unit=44)
c
      write(6,'(a)')
     &   'PRINT MESSAGE: Atmospheric parameters '//
     &   'have been written to the directory AtmosPars'
      write(6,'(15X,A)')
     &   'in the files: Temperature.dat, MassDensity.dat, etc.'
      end
C     *************************************************************************
      subroutine PrintWAvenumbersAndEigenvaluesSF
     i         ( TypeLinModel, NModes, nzWT, IdzT, izT, zmin, dz,
     i            KRCT, EvalT )
      implicit none
      integer    TypeLinModel, NModes, nzWT, IdzT
      integer    izT(nzWT+1)
      real*8     zmin, dz
      complex*16 KRCT(2*NModes, nzWT)
      complex*16 EValT(2*NModes, nzWT)
c
      integer    i, iT, j, M
      real*8     z
      character*1  charM
      character*72 filename, path
c
      M = 2*NModes
c
      write(charM,'(I1)') TypeLinModel

      path = 'EigenValWaveNumSF/'
c     -------------------------------------------------------------------------
c                                   Eigenvalues
c     -------------------------------------------------------------------------
      filename =
     &  TRIM(path)//'EigenValRealLinModel'//TRIM(charM)//'.dat'
      open(unit=70,file=trim(filename),status='unknown')
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        z = zmin + (i-1)*dz
        write(70,'(1pe11.4,3x,7(1pe13.4,2x))')
     &        z, ( dble(EValT(j,iT)),j = 1,M ), 0.0
      end do
      close(unit=70)
C
      filename =
     &  TRIM(path)//'EigenValRealAscDescGWLinModel'//
     &  TRIM(charM)//'.dat'
      open(unit=70,file=trim(filename),status='unknown')
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        z = zmin + (i-1)*dz
        write(70,'(1pe11.4,3x,2(1pe13.4,2x))')
     &        z, dble(EValT(1,iT)), dble(EValT(4,iT))
      end do
      close(unit=70)
C
      filename =
     &  TRIM(path)//'EigenValImagLinModel'//TRIM(charM)//'.dat'
      open(unit=70,file=trim(filename),status='unknown')
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        z = zmin + (i-1)*dz
        write(70,'(1pe11.4,3x,7(1pe13.4,2x))')
     &        z, ( aimag(EValT(j,iT)),j = 1,M ), 0.0
      end do
      close(unit=70)
c     -------------------------------------------------------------------------
C                                Wavenumbers
c     -------------------------------------------------------------------------
      filename =
     & TRIM(path)//'WavenumberRealLinModel'//TRIM(charM)//'.dat'
      open(unit=70,file=trim(filename),status='unknown')
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        z = zmin + (i-1)*dz
        write(70,'(1pe11.4,3x,7(1pe13.4,2x))')
     &        z, ( dble(KRCT(j,iT)),j = 1,M ), 0.0
      end do
      close(unit=70)
C
      filename =
     &  TRIM(path)//'WavenumberImagAscDescGWLinModel'//
     &  TRIM(charM)//'.dat'
      open(unit=70,file=trim(filename),status='unknown')
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        z = zmin + (i-1)*dz
        write(70,'(1pe11.4,3x,2(1pe13.4,2x))')
     &        z, aimag(KRCT(1,iT)), aimag(KRCT(4,iT))
      end do
      close(unit=70)
C
      filename =
     &   TRIM(path)//'WavenumberImagLinModel'//TRIM(charM)//'.dat'
      open(unit=70,file=trim(filename),status='unknown')
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        z = zmin + (i-1)*dz
        write(70,'(1pe11.4,3x,7(1pe13.4,2x))')
     &        z, ( aimag(KRCT(j,iT)),j = 1,M ), 0.0
      end do
      close(unit=70)

      write(6,'(a)')
     &   'PRINT MESSAGE: Eigenvalues and vertical wavenumbers '//
     &   'have been written to the directory WaveNumbersSF'
      write(6,'(15x,a)') 'in the files: EigenValRealLinModelX.dat, '//
     &   'EigenValRealAscDescGWLinModelX.dat, etc.'
      end
C     *************************************************************************
      subroutine PrintWaveParametersModelGMMAComplete
     i         ( TypeLinModel, nzWT, izT, zmin, dz, OmegaIm,
     i           UThat, UThat0, UThatPlus, UThatMinus,
     i           URhat, URhat0, URhatPlus, URhatMinus,
     i           TWhat, TWhat0, TWhatPlus, TWhatMinus,
     i           pWhat, roWhat )
      implicit none
      integer TypeLinModel, nzWT
      integer izT(nzWT+1)
      real*8  zmin, dz, OmegaIm

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

c     --- local ---
      integer iT, i
      real*8  z, MaxP
      logical PrintPressDens
      character*1  charM
      character*72 path, path0, filename
c
      PrintPressDens = .false.

      write(charM,'(I1)') TypeLinModel
c
      path  = 'WaveResultsSF/'
      path0 = 'SolMetGMMALinModel'//TRIM(charM)//'.dat'
!     -------------------------------------------------------------------------
!                             Horizontal velocity
!     -------------------------------------------------------------------------
      filename = TRIM(path)//'UVelocity'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)') '# Horizontal velocity versus altitude'
      MaxP = -1.d+20
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(UThat(iT)), z
        MaxP = max( MaxP, abs(dble(UThat(iT))) )
      end do
      write(6,*)
      write(6,'(a,1pe11.4)')
     &'Single-frequency wave results for OmegaIm = ', OmegaIm
      write(6,'(4x,a,1pe11.4)')'- Maximum horizontal velocity [m/s] = ',
     &      MaxP
      close(13)

      filename = TRIM(path)//'UVelocityAscPlusDesc'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)')
     & '# Sum of ascending and descending horizontal '//
     & 'velocities versus altitude'
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(UThat0(iT)), z
      end do
      close(13)

      filename = TRIM(path)//'UVelocityAsc'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)') '# Ascending horizontal velocity versus altitude'
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(UThatPlus(iT)), z
      end do
      close(13)

      filename = TRIM(path)//'UVelocityDesc'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)')'# Desecending horizontal velocity versus altitude'
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(UThatMinus(iT)), z
      end do
      close(13)
!     -------------------------------------------------------------------------
!                            Vertical velocity
!     -------------------------------------------------------------------------
      filename = TRIM(path)//'WVelocity'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)') '# Vertical velocity versus altitude'
      MaxP = -1.d+20
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(URhat(iT)), z
        MaxP = max( MaxP, abs(dble(URhat(iT))) )
      end do
      write(6,'(4x,a,1pe11.4)')'- Maximum vertical velocity [m/s]   = ',
     &      MaxP
      close(13)

      filename = TRIM(path)//'WVelocityAscPlusDesc'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)')
     & '# Sum of ascending and descending vertical '//
     & 'velocities versus altitude'
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(URhat0(iT)), z
      end do
      close(13)

      filename = TRIM(path)//'WVelocityAsc'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)') '# Ascending vertical velocity versus altitude'
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(URhatPlus(iT)), z
      end do
      close(13)

      filename = TRIM(path)//'WVelocityDesc'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)') '# Desecending vertical velocity versus altitude'
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(URhatMinus(iT)), z
      end do
      close(13)
!     -------------------------------------------------------------------------
!                             Temperature
!     -------------------------------------------------------------------------
      filename = TRIM(path)//'Temperature'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)') '# Temperature versus altitude'
      MaxP = -1.d+20
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(TWhat(iT)), z
        MaxP = max( MaxP, abs(dble(TWhat(iT))) )
      end do
      write(6,'(4x,a,1pe11.4)')'- Maximum temperature [K]           = ',
     &      MaxP
      close(13)

      filename = TRIM(path)//'TemperatureAscPlusDesc'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)')
     & '# Sum of ascending and descending wave temperatures '//
     & 'versus altitude'
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(TWhat0(iT)), z
      end do
      close(13)

      filename = TRIM(path)//'TemperatureAsc'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)') '# Ascending wave temperature versus altitude'
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(TWhatPlus(iT)), z
      end do
      close(13)

      filename = TRIM(path)//'TemperatureDesc'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)') '# Desecending wave temperature versus altitude'
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(TWhatMinus(iT)), z
      end do
      close(13)
!     -----------------------------------------------------------------------
!                                  Pressure
!     -----------------------------------------------------------------------
      if ( PrintPressDens ) then
        filename = TRIM(path)//'Pressure'//TRIM(path0)
        open(unit=13,file=trim(filename),status='unknown')
        write(13,'(a)') '# Pressure versus altitude'
        do iT = 2, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz
          write(13,'(4(2x,1pe18.8))') dble(pWhat(iT)), z
        end do
        close(13)
      end if
!     -----------------------------------------------------------------------
!                                Mass density
!     -----------------------------------------------------------------------
      if ( PrintPressDens ) then
        filename = TRIM(path)//'MassDensity'//TRIM(path0)
        open(unit=13,file=trim(filename),status='unknown')
        write(13,'(a)') '# Mass density versus altitude'
        do iT = 2, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz
          write(13,'(4(2x,1pe18.8))') dble(roWhat(iT)), z
        end do
        close(13)
      end if
c
      write(6,*)
      write(6,'(a)')
     &   'PRINT MESSAGE: Wave parameters '//
     &   '(velocities, temperature, pressure, and mass density) '//
     &   'have been written'
      write(6,'(15x,a)')
     &   'to the directory WaveResultsSF in the files: '//
     &   'UVelocitySolMetGMMALinModelX.dat, etc.'
      end
c     *************************************************************************
      subroutine PrintWaveParametersModelGMMNandSMMA
     i         ( TypeLinModel, TypeSolMet, nzWT, izT, zmin, dz, OmegaIm,
     i           UThat, URhat, TWhat, pWhat, roWhat )
      implicit none
      integer TypeLinModel, TypeSolMet, nzWT
      integer izT(nzWT+1)
      real*8  zmin, dz, OmegaIm
C
      complex*16 UThat(nzWT)
      complex*16 URhat(nzWT)
      complex*16 TWhat(nzWT)
      complex*16 pWhat(nzWT)
      complex*16 roWhat(nzWT)

c     --- local ---
      integer iT, i
      real*8  z, MaxP
      logical PrintPressDens
      character*1  charM
      character*72 path, path0, filename
c
      PrintPressDens = .false.
c
      write(charM,'(I1)') TypeLinModel
c
      path  = 'WaveResultsSF/'
      if (TypeSolMet .eq. 2 ) then
        path0 = 'SolMetGMMNLinModel'//TRIM(charM)//'.dat'
      else
        path0 = 'SolMetSMMALinModel'//TRIM(charM)//'.dat'
      end if
c     -------------------------------------------------------------------------
C                             Horizontal velocity
c     -------------------------------------------------------------------------
      filename = TRIM(path)//'UVelocity'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)') '# Horizontal velocity versus altitude'
      MaxP = -1.d+20
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(UThat(iT)), z
        MaxP = max( MaxP, abs(dble(UThat(iT))) )
      end do
      write(6,*)
      write(6,'(a,1pe11.4)')
     &'Single-frequency wave results for OmegaIm = ', OmegaIm
      write(6,'(4x,a,1pe11.4)')'- Maximum horizontal velocity [m/s] = ',
     &      MaxP
      close(13)
c     -------------------------------------------------------------------------
C                              Vertical velocity
c     -------------------------------------------------------------------------
      filename = TRIM(path)//'WVelocity'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)') '# Vertical velocity versus altitude'
      MaxP = -1.d+20
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(URhat(iT)), z
        MaxP = max( MaxP, abs(dble(URhat(iT))) )
      end do
      write(6,'(4x,a,1pe11.4)')'- Maximum vertical velocity [m/s]   = ',
     &      MaxP
      close(13)
c     -------------------------------------------------------------------------
C                                Temperature
c     -------------------------------------------------------------------------
      filename = TRIM(path)//'Temperature'//TRIM(path0)
      open(unit=13,file=trim(filename),status='unknown')
      write(13,'(a)') '# Temperature versus altitude'
      MaxP = -1.d+20
      do iT = 2, nzWT
        i = izT(iT)
        z = zmin + (i-1)*dz
        write(13,'(4(2x,1pe18.8))') dble(TWhat(iT)), z
        MaxP = max( MaxP, abs(dble(TWhat(iT))) )
      end do
      write(6,'(4x,a,1pe11.4)')'- Maximum temperature [K]           = ',
     &      MaxP
      close(13)
c     -------------------------------------------------------------------------
C                                 Pressure
c     -------------------------------------------------------------------------
      if ( PrintPressDens ) then
        filename = TRIM(path)//'Pressure'//TRIM(path0)
        open(unit=13,file=trim(filename),status='unknown')
        write(13,'(a)') '# Pressure versus altitude'
        do iT = 2, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz
          write(13,'(4(2x,1pe18.8))') dble(pWhat(iT)), z
        end do
        close(13)
      end if
c     -------------------------------------------------------------------------
C                               Mass density
c     -------------------------------------------------------------------------
      if ( PrintPressDens ) then
        filename = TRIM(path)//'MassDensity'//TRIM(path0)
        open(unit=13,file=trim(filename),status='unknown')
        write(13,'(a)') '# Mass density versus altitude'
        do iT = 2, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz
          write(13,'(4(2x,1pe18.8))') dble(roWhat(iT)), z
        end do
        close(13)
      end if
c
      write(6,*)
      write(6,'(a)')
     &   'PRINT MESSAGE: Wave parameters '//
     &   '(velocities, temperature, pressure, and mass density) '//
     &   'have been written'
      if (TypeSolMet .eq. 2 ) then
        write(6,'(15x,a)')
     &   'to the directory WaveResultsSF in the files '//
     &   'UVelocitySolMetGMMNLinModelX.dat, etc.'
      else
        write(6,'(15x,a)')
     &   'to the directory WaveResultsSF in the files '//
     &   'UVelocitySolMetSMMALinModelX.dat, etc.'
      end if
      end
c     *************************************************************************
      subroutine PrintWaveParsModelWP
     i         ( NFFT, nzWT, izT, zmin, dz,
     i           TimeShift, Omega0, RatioSigmaOmega, TimeK,
     i           UTTime, URTime, TTime  )
      implicit none
      integer NFFT, nzWT
      integer izT(nzWT+1)
      real*8  zmin, dz
      real*8  TimeShift, Omega0, RatioSigmaOmega
      real*8  TimeK( NFFT )
c
      complex*16 UTTime(nzWT, NFFT)
      complex*16 URTime(nzWT, NFFT)
      complex*16 TTime(nzWT, NFFT)
C
C     --- locals ---
      integer i, j, iT, JMax, iTMax
      real*8  z, SigmaOmega, MaxP, Src
      character*72 path, filename
      complex*16 SrcTimeAmplit
c
      path  = 'WaveResultsWP/'
c
c     --- omega standard deviation ---
      SigmaOmega = Omega0 / RatioSigmaOmega
c     -------------------------------------------------------------------------
c                            Horizontal velocity
c     -------------------------------------------------------------------------
      filename = TRIM(path)//'UVelocityWP.dat'
      open(unit=44,file=trim(filename),status='unknown')
      MaxP = -1.d+20
      do j = 1, NFFT
        do iT = 2, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz !in km
          write(44,'(3(1pe13.4,2x))')
     &    TimeK(j)/3600.0, z, dble(UTTime(iT, j))
c
          if ( dble(UTTime(iT, j)) .gt. MaxP ) then
            MaxP  = dble(UTTime(iT, j))
            JMax  = j
            iTMax = iT
          end if
        end do
        write(44,*)
      end do
      close(44)
c
c     --- gnuplot ---
      CALL SYSTEM('gnuplot -p GnuPlot/data_plotU.plt')
c
      filename = TRIM(path)//'UVelocityWPMax.dat'
      open(unit=44,file=trim(filename),status='unknown')
      do j = 1, NFFT
        write(44,'(2(1pe13.4,2x))')
     &  TimeK(j)/3600.0, dble(UTTime(iTMax, j))
      end do
      close(44)
c
      z = zmin + ( izT(iTMax) - 1 )*dz
      write(6,'(4x,3(a,1pe11.4))')
     &'Maximum horizontal velocity [m/s] = ', MaxP,
     &' at time [hr] = ', TimeK(JMax)/3600.d0,
     &' and level [km] = ', z
c     -------------------------------------------------------------------------
c                             Vertical velocity
c     -------------------------------------------------------------------------
      filename = TRIM(path)//'WVelocityWP.dat'
      open(unit=44,file=trim(filename),status='unknown')
      MaxP = -1.d+20
      do j = 1, NFFT
        do iT = 2, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz !in km
          write(44,'(3(1pe13.4,2x))')
     &    TimeK(j)/3600.0, z, dble(URTime(iT, j))
c
          if ( dble(URTime(iT, j)) .gt. MaxP ) then
            MaxP  = dble(URTime(iT, j))
            JMax  = j
            iTMax = iT
          end if
        end do
        write(44,*)
      end do
      close(44)
c
c     --- gnuplot ---
      CALL SYSTEM('gnuplot -p GnuPlot/data_plotW.plt')
c
      filename = TRIM(path)//'WVelocityWPMax.dat'
      open(unit=44,file=trim(filename),status='unknown')
      do j = 1, NFFT
        write(44,'(2(1pe13.4,2x))')
     &  TimeK(j)/3600.0, dble(URTime(iTMax, j))
      end do
      close(44)
c
      z = zmin + ( izT(iTMax) - 1 )*dz
      write(6,'(4x,3(a,1pe11.4))')
     &'Maximum vertical velocity [m/s]   = ', MaxP,
     &' at time [hr] = ', TimeK(JMax)/3600.d0,
     &' and level [km] = ', z
c     -------------------------------------------------------------------------
c                             Temperature
c     -------------------------------------------------------------------------
      filename = TRIM(path)//'TemperatureWP.dat'
      open(unit=44,file=trim(filename),status='unknown')
      MaxP = -1.d+20
      do j = 1, NFFT
        do iT = 2, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz !in km
          write(44,'(3(1pe13.4,2x))')
     &    TimeK(j)/3600.0, z, dble(TTime(iT, j))
c
          if ( dble(TTime(iT, j)) .gt. MaxP ) then
            MaxP  = dble(TTime(iT, j))
            JMax  = j
            iTMax = iT
          end if
        end do
        write(44,*)
      end do
      close(44)
c
c     --- gnuplot ---
      CALL SYSTEM('gnuplot -p GnuPlot/data_plotT.plt')
c
      filename = TRIM(path)//'TemperatureWPMax.dat'
      open(unit=44,file=trim(filename),status='unknown')
      do j = 1, NFFT
        write(44,'(2(1pe13.4,2x))')
     &  TimeK(j)/3600.0, dble(TTime(iTMax, j))
      end do
      close(44)
c
      z = zmin + ( izT(iTMax) - 1 )*dz
      write(6,'(4x,3(a,1pe11.4))')
     &'Maximum temperature [K]           = ', MaxP,
     &' at time [hr] = ', TimeK(JMax)/3600.d0,
     &' and level [km] = ', z
c     -------------------------------------------------------------------------
C                               Source function
c     -------------------------------------------------------------------------
      filename = TRIM(path)//'SrcFctWP.dat'
      open(unit=44,file=trim(filename),status='unknown')
      do j = 1, NFFT
        do iT = 2, nzWT
          i = izT(iT)
          z = zmin + (i-1)*dz !in km
          Src = dble( SrcTimeAmplit( TimeK(j), TimeShift,
     &                Omega0, SigmaOmega ) )
          write(44,'(3(1pe13.4,2x))')
     &    TimeK(j)/3600.0, z, Src
        end do
        write(44,*)
      end do
      close(44)
C
C     --- gnuplot ---
      CALL SYSTEM('gnuplot -p GnuPlot/data_plotSignal.plt')
c
      filename = TRIM(path)//'SrcFctWPMax.dat'
      open(unit=44,file=trim(filename),status='unknown')
      MaxP = -1.d+20
      do j = 1, NFFT
        Src = dble( SrcTimeAmplit( TimeK(j), TimeShift,
     &              Omega0, SigmaOmega ) )
        write(44,'(2(1pe13.4,2x))')
     &  TimeK(j)/3600.0, Src

        if ( Src .gt. MaxP ) then
          MaxP  = Src
          JMax  = j
        end if
      end do
      close(44)
c
      write(6,'(4x,3(a,1pe11.4))')
     &'Maximum source function           = ', MaxP,
     &' at time [hr] = ', TimeK(JMax)/3600.d0
c
      write(6,*)
      write(6,'(A)') 'PRINT MESSAGE: Plots generated with Gnuplot '//
     &               'for velocities and temperature have been written'
      write(6,'(15X,A)')
     &    'to the directory WaveResultsWP in the files: '//
     &    'UVelocityWP.png, etc.'
      end
C     *************************************************************************
      subroutine PrintWAvenumbersAndEigenvaluesWPGlobalCausality
     i         ( TypeLinModel, NModes, nzWT, IdzT, izT, zmin, dz,
     i            KRCT, EvalT, OmegaIm )
      implicit none
      integer    TypeLinModel, NModes, nzWT, IdzT
      integer    izT(nzWT+1)
      real*8     zmin, dz, OmegaIm
      complex*16 KRCT(2*NModes, nzWT)
      complex*16 EValT(2*NModes, nzWT)
c
      integer    i, iT, M
      real*8     z
      character*1  charM
      character*72 filename, path
c
      M = 2*NModes
c
      write(charM,'(I1)') TypeLinModel
      path = 'EigenValWaveNumWP/'
C
      filename = TRIM(path)//'EigenValRealAscDescGWCausalityLinModel'//
     &           TRIM(charM)//'.dat'
      open(unit=70,file=trim(filename),status='unknown')
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        z = zmin + (i-1)*dz
        write(70,'(1pe11.4,3x,2(1pe13.4,2x))')
     &        z, dble(EValT(1,iT)), dble(EValT(4,iT))
      end do
      close(unit=70)
C
      filename=TRIM(path)//'WavenumberImagAscDescGWCausalityLinModel'//
     &         TRIM(charM)//'.dat'
      open(unit=70,file=trim(filename),status='unknown')
      do iT = 1, nzWT
        i = izT(iT) + IdzT / 2
        z = zmin + (i-1)*dz
        write(70,'(1pe11.4,3x,2(1pe13.4,2x))')
     &        z, aimag(KRCT(1,iT)), aimag(KRCT(4,iT))
      end do
      close(unit=70)
      write(6,'(a,1pe11.4,a)')
     &   'PRINT MESSAGE: Eigenvalues and vertical wavenumbers versus '//
     &   'altitude, for OmegaIm = ', OmegaIm, ' satisfying the'
      write(6,'(15x,a)')
     &   'strong causality condition have been '//
     &   'written to the directory EigenValWaveNumWP in the files:'
      write(6,'(15x,a)')
     &   'EigenValRealAscDescGWCausalityLinModelX.dat and '//
     &   'WavenumberImagAscDescGWCausalityLinModelX.dat'
      end
