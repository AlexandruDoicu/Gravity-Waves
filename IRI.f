c     *************************************************************************
c                                IRI test files:
c     1. Latitude and longitudes:
c        A. Jicamarca(Peru)	                     −12   283   300 km
c        B. Arecibo(Puerto Rico)	             +18   293   300 km
c        C. Millstone Hill(USA)	                 +42   288   300 km
c        D. Saint-Santin(France)	             +44     2   300 km
c        E. EISCAT TromsøEISCAT Tromsø (Auroral) +70    19   300 km
c
c
c     2. Date and Time: 11.02.2012; 10:00, i.e.,
c       - year   = 2012
c       - month  = 2
c       - day    = 11
c       - hour   = 10
c       - minute = 0
c       - second = 0
c       All variables are integers
c
c     3. Altitudes: zmin = 80 km, zmax = 500 km, dz = 1.0 km
c     *************************************************************************
      subroutine test_read_iri
      implicit none
      integer   NZW, I
      real*8    ZMIN, ZMAX, DZ, Z
      real      F10P7, FBAR
      real*8    DipDeg
      real*8, allocatable:: NOW_M3( : )
      character*200 FNAME
c
      integer YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
      real    LatRef, LonRef
      integer IYD
      real    SEC, HRL
c
c     --- set your file name
      FNAME = 'IRIData/iri_outputE.txt'
      ZMIN  = 80.d0
      ZMAX  = 500.d0
      NZW   = 801
c
      allocate( NOW_M3(NZW) )
      call READ_IRI_IONFILE
     i   ( FNAME,
     i     NZW, ZMIN, ZMAX,
     o     YEAR, MONTH, DAY, HOUR, MINUTE, SECOND,
     o     LatRef, LonRef,
     o     DipDeg, F10P7, FBAR,
     o     NOW_M3 )
c
      call IRI2HWM07_SETUP
     i   ( YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, LonRef,
     o     IYD, SEC, HRL )
c
      write(6,'(a,i4,1x,i2,1x,i2)') 'Date (Y M D) = ',
     &      YEAR, MONTH, DAY
      write(6,'(a,i2,1x,i2,1x,i2)') 'Time (H M S) = ',
     &      HOUR, MINUTE, SECOND
      write(6,'(a,i5)')    'IYD   (YYDDD)        = ', IYD
      write(6,'(a,f12.3)') 'SEC   (UT seconds)   = ', SEC
      write(6,'(a,f12.6)') 'HRL   (solar LT hrs) = ', HRL
      write(6,'(a,f10.3)') 'LatRef = ', LatRef
      write(6,'(a,f10.3)') 'LonRef = ', LonRef
      write(6,'(a,f10.3)') 'DipDeg = ', DipDeg
      write(6,'(a,f10.3)') 'F10P7  = ', F10P7
      write(6,'(a,f10.3)') 'FBAR   = ', FBAR
C
c     --- print as a check ---
      DZ = (ZMAX - ZMIN) / REAL(NZW - 1)
      do I = 1, NZW
        Z = ZMIN + (I-1)*DZ
         write(6,'(I5, 2X, f8.1,1x,1pe12.4)') I, Z, NOW_M3(I)
         if ( mod(I,20).eq.0 ) read(5,*)
      end do

      deallocate( NOW_M3 )
      end
c     *************************************************************************
      subroutine READ_IRI_IONFILE
     i         ( FNAME,
     i           NZW, ZMIN, ZMAX,
     o           YEAR, MONTH, DAY, HOUR, MINUTE, SECOND,
     o           LatRef, LonRef,
     o           DipDeg, F10P7, FBAR,
     O           NOW_M3 )
c     -------------------------------------------------------------------------
c     Read IRI ASCII output (like iri_output.txt) and compute ion densities.
c
c     Reads:
c       Dip (deg)
c       F10.7 daily
c       F10.7 81-day average
c       Table: altitude, Ne (cm-3), and ion percentages [%]*10
c
c     Computes (1/m^3):
c        Ne_m3, n(O+), n(N+), n(H+), n(He+), n(O2+), n(NO+), n(Cluster)
c
c     IMPORTANT (matches IRI table header "ION PERCENTAGES[%]*10"):
c        fraction(species) = ion_col / 1000.0
c        => n_species = Ne_m3 * ion_col / 1000.0
c
c     -------------------------------------------------------------------------
      implicit none
c
c     --- arguments ---
      integer  NZW
      real*8   ZMIN, ZMAX
      character*(*) FNAME
C
      integer YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
      real    LatRef, LonRef

      real   F10P7, FBAR
      real*8 DipDeg
      real*8 NOW_M3(NZW)
c
c     ---- locals ---
      integer NMAX, N, IU, IOS
      character*200 LINE
      logical INTABLE, FOUNDHDR, ISBLANK
      integer P_O, P_N, P_H, P_HE, P_O2, P_NO, P_CL
      real    AltRef, Z, NECM3, DUM, TN, TI, TE, TEC, TPCT, FAC
      real    DipDegR
      real*8  DZ
      real, allocatable:: NE_M3(:), NN_M3(:), NH_M3(:), NHE_M3(:)
      real, allocatable:: NO2_M3(:), NNO_M3(:), NCL_M3(:)

      real*8, allocatable:: ZKM(:), NO_M3(:)
c
      NMAX = 5000
      allocate( NE_M3(NMAX), NN_M3(NMAX), NH_M3(NMAX), NHE_M3(NMAX) )
      allocate( NO2_M3(NMAX), NNO_M3(NMAX), NCL_M3(NMAX) )
      allocate( ZKM(NMAX), NO_M3(NMAX) )
c
      N = 0
      DipDegR = -999.0
      F10P7   = -999.0
      FBAR    = -999.0
      INTABLE = .false.
      FOUNDHDR= .false.
c
c     --- open file ---
      IU = 10
      open(unit=IU, file=FNAME, status='old', iostat=IOS)
      if (IOS .ne. 0) then
         write(6,'(a)') 'Error when opening the IRI-data file'
         stop
      end if
c
c     --- find first NON-empty line (header line); skip leading blanks ---
 5    continue
      read(IU,'(A)',iostat=IOS) LINE
      if (IOS .ne. 0) then
         write(6,'(a)')
     &  'Error when reading the first line in the IRI-data file'
         stop
      end if
      if (ISBLANK(LINE)) goto 5
c
c     --- parse that first non-empty line ---
      call PARSE_IRI_FIRSTLINE
     &   ( LINE,
     &     YEAR, MONTH, DAY, HOUR, MINUTE, SECOND,
     &     LatRef, LonRef, AltRef )
      FOUNDHDR = .true.
c
c     --- scan file line by line ---
 10   continue
      read(IU,'(A)',iostat=IOS) LINE
      if (IOS .ne. 0) goto 90
c
c     --- pick scalar values ---
      if (index(LINE,'Dip (Magnetic Inclination)/degree') .gt. 0) then
         call READ_LAST_REAL(LINE, DipDegR)
      else if (index(LINE,'Solar radio flux F10.7 (daily)') .gt. 0) then
         call READ_LAST_REAL(LINE, F10P7)
      else if (index(LINE,
     & 'Solar radio flux F10.7 (81-day average)') .gt. 0) then
         call READ_LAST_REAL(LINE, FBAR)
      end if
c
c     --- detect start of table (header line with "km  Ne/cm-3") ---
      if (.not. INTABLE) then
         if (index(LINE,'km  Ne/cm-3') .gt. 0) then
            INTABLE = .true.
c           next line(s) may be dashed; keep reading
         end if
         goto 10
      end if
c
c     --- in table: skip empty/separator lines ---
c     if (LINE .eq. ' ') goto 10
      if (ISBLANK(LINE)) goto 10
      if (index(LINE,'-') .eq. 1) goto 10
c     -------------------------------------------------------------------------
c                  try to read one data row (list-directed)
c     Expected columns (from IRI):
c       z  Ne/cm3  Ne/NmF2  Tn Ti Te  O+ N+ H+ He+ O2+ NO+ Clust  TEC  t/%
c     -------------------------------------------------------------------------
      IOS = 0
      read(LINE,*,iostat=IOS) Z, NECM3, DUM, TN, TI, TE,
     &                        P_O, P_N, P_H, P_HE, P_O2, P_NO, P_CL,
     &                        TEC, TPCT
      if (IOS .ne. 0) then
c        if we can no longer parse numeric rows, table ended
         goto 90
      end if

      if (N .ge. NMAX) then
         write(6,'(a)')
     &  'Number of data in the IRI file is larger than assumed'
         stop
      end if

      N = N + 1
      ZKM(N) = dble(Z)
c
c     --- Ne: cm^-3 -> m^-3 ---
      NE_M3(N) = NECM3 * 1.0E6
c
c     --- IRI ion columns are [%]*10, so fraction = col/1000 ---
      FAC = NE_M3(N) / 1000.0
c
c     --- treat negative percentages (e.g. -1) as 0 ---
      if (P_O  .lt. 0) P_O  = 0
      if (P_N  .lt. 0) P_N  = 0
      if (P_H  .lt. 0) P_H  = 0
      if (P_HE .lt. 0) P_HE = 0
      if (P_O2 .lt. 0) P_O2 = 0
      if (P_NO .lt. 0) P_NO = 0
      if (P_CL .lt. 0) P_CL = 0

      NO_M3(N)  = dble( FAC * P_O )
      NN_M3(N)  = FAC * real(P_N)
      NH_M3(N)  = FAC * real(P_H)
      NHE_M3(N) = FAC * real(P_HE)
      NO2_M3(N) = FAC * real(P_O2)
      NNO_M3(N) = FAC * real(P_NO)
      NCL_M3(N) = FAC * real(P_CL)

      goto 10

 90   continue
      close(IU)
C
C     --- interpolate on the output grid ---
      DZ = (ZMAX - ZMIN) / REAL(NZW - 1)
      call LININT_UNIFORM
     i   ( N, ZKM, NO_M3,
     &     NZW, ZMIN, DZ,
     &     NOW_M3 )
c
c     --- pass to double ---
      DipDeg = dble(DipDegR)
c
      deallocate( NE_M3, NN_M3, NH_M3, NHE_M3 )
      deallocate( NO2_M3, NNO_M3, NCL_M3 )
      deallocate( ZKM, NO_M3 )
      return
      end
c     *************************************************************************
      logical function ISBLANK(LINE)
c     -------------------------------------------------------------------------
c             Return .true. if LINE is blank or only spaces/tabs
c      ------------------------------------------------------------------------
      implicit none
      character*(*) LINE
      integer I, N
      N = len(LINE)
      ISBLANK = .true.
      do I = 1, N
         if (LINE(I:I) .ne. ' ' .and. LINE(I:I) .ne. char(9)) then
            ISBLANK = .false.
            return
         end if
      end do
      return
      end
c     *************************************************************************
      subroutine READ_LAST_REAL(LINE, X)
      implicit none
      character*(*) LINE
      real X
      integer L, I, J, IOS
      character*64 TOK

      X = -999.0
      L = len(LINE)
c
c     --- scan from end to find start of last token ----
      I = L
 10   if (I .le. 1) goto 90
      if (LINE(I:I) .ne. ' ') goto 20
      I = I - 1
      goto 10

 20   continue
      J = I
 30   if (J .le. 1) goto 40
      if (LINE(J:J) .eq. ' ') goto 40
      J = J - 1
      goto 30

 40   continue
      TOK = ' '
      TOK = LINE(J+1:I)

      IOS = 0
      read(TOK,*,iostat=IOS) X
      if (IOS .ne. 0) X = -999.0

 90   continue
      return
      end
c     *************************************************************************
      subroutine LININT_UNIFORM
     i         ( N, X, Y, NZW, ZMIN, DZ,
     o           YW )
c     -------------------------------------------------------------------------
c     Linear interpolation from irregular grid (x(i),y(i)) onto uniform grid.
c
c     Inputs:
c       N      number of input points
c       X(N)   input grid (must be increasing)
c       Y(N)   input values
c       NZW    number of output points
c       ZMIN   output grid start
c       DZ     output grid spacing
c
c     Outputs:
c       YW(NZW)     interpolated values
c
c      Behavior:
c      - uses piecewise linear interpolation
c      - endpoints match input endpoints by construction (ZMIN=X(1), ZMAX=X(N))
c     -------------------------------------------------------------------------
      implicit none
      integer N, NZW
      real*8  X(N), Y(N), ZMIN, DZ
      real*8  YW(NZW)

      integer I, K
      real*8  Z, T

      K = 1

      do I = 1, NZW
         Z = ZMIN + REAL(I-1)*DZ
c
c        ---- advance bracket index K so that X(K) <= Z <= X(K+1) ---
 10      if ( K .lt. N - 1 ) then
            if ( Z .gt. X(K+1) ) then
               K = K + 1
               goto 10
            end if
         end if
c
c        --- handle edges (should be exact by construction, but keep safe) ---
         if (Z .le. X(1)) then
            YW(I) = Y(1)
         else if (Z .ge. X(N)) then
            YW(I) = Y(N)
c
c        --- linear interpolation in interval [X(K),X(K+1)] ----
         else
            T = (Z - X(K)) / (X(K+1) - X(K))
            YW(I) = (1.d0 - T)*Y(K) + T*Y(K+1)
         end if
      end do
      return
      end
c     *************************************************************************
      subroutine PARSE_IRI_FIRSTLINE
     &         ( LINE,
     &           YEAR, MONTH, DAY, HOUR, MINUTE, SECOND,
     &           LatRef, LonRef, AltRef )
c     -------------------------------------------------------------------------
c     Parse the first non-empty IRI line:
c          ...):YYYY/ <MMDD or -DDD>/<HH.H>UT  geog Lat/Long/Alt= lat/lon/alt
c     -------------------------------------------------------------------------
      implicit none
      character*(*) LINE
      integer YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
      real    LatRef, LonRef, AltRef

      integer I, J, K, L
      integer DATECODE, DOY
      real    UTH, X
      character*200 S
      character*60  TOK

      S = LINE

      I = index(S, '):')
      if (I .le. 0) then
         write(6,'(a)') 'PARSE_IRI_FIRSTLINE: cannot find "):"'
         write(6,'(a)') LINE
         stop
      end if
      J = I + 2

      K = index(S(J:), 'UT')
      if (K .le. 0) then
         write(6,'(a)') 'PARSE_IRI_FIRSTLINE: cannot find "UT"'
         write(6,'(a)') LINE
         stop
      end if

      TOK = ' '
      TOK = S(J:J+K-2)
      call REPLCHAR(TOK, '/', ' ')

      read(TOK,*,err=900) YEAR, DATECODE, UTH

      if (DATECODE .lt. 0) then
         DOY = -DATECODE
         call DOY_TO_MD(YEAR, DOY, MONTH, DAY)
      else
         MONTH = DATECODE / 100
         DAY   = mod(DATECODE,100)
      end if

      if (UTH .lt. 0.0) UTH = 0.0
      HOUR = int(UTH)
      X = (UTH - real(HOUR))*60.0
      MINUTE = int(X)
      SECOND = int((X - real(MINUTE))*60.0 + 0.5)

      if (SECOND .ge. 60) then
         SECOND = SECOND - 60
         MINUTE = MINUTE + 1
      end if
      if (MINUTE .ge. 60) then
         MINUTE = MINUTE - 60
         HOUR   = HOUR + 1
      end if
      if (HOUR .ge. 24) then
         HOUR = HOUR - 24
      end if

      L = index(S, 'geog Lat/Long/Alt=')
      if (L .le. 0) then
         write(6,'(a)')
     &  'PARSE_IRI_FIRSTLINE: cannot find geog Lat/Long/Alt='
         write(6,'(a)') LINE
         stop
      end if
      L = L + len('geog Lat/Long/Alt=')

      TOK = ' '
      TOK = S(L:)
      call REPLCHAR(TOK, '/', ' ')
      read(TOK,*,err=901) LatRef, LonRef, AltRef

      return

 900  continue
      write(6,'(a)')
     &'PARSE_IRI_FIRSTLINE: error parsing YEAR/DATECODE/UTH'
      write(6,'(a)') LINE
      stop
 901  continue
      write(6,'(a)') 'PARSE_IRI_FIRSTLINE: error parsing Lat/Lon/Alt'
      write(6,'(a)') LINE
      stop
      end
c     *************************************************************************
      subroutine REPLCHAR(STR, COLD, CNEW)
c     -------------------------------------------------------------------------
c             Replace all occurrences of COLD with CNEW in STR (in place)
c     -------------------------------------------------------------------------
      implicit none
      character*(*) STR
      character*1 COLD, CNEW
      integer I, N
      N = len(STR)
      do I = 1, N
         if (STR(I:I) .eq. COLD) STR(I:I) = CNEW
      end do
      return
      end
c     *************************************************************************
      subroutine DOY_TO_MD(YEAR, DOY, MONTH, DAY)
c     -------------------------------------------------------------------------
c                              DOY -> (MONTH,DAY)
c     -------------------------------------------------------------------------
      implicit none
      integer YEAR, DOY, MONTH, DAY
      integer NDAYS(12), M, DLEFT
      logical LEAP
      data NDAYS /31,28,31,30,31,30,31,31,30,31,30,31/

      call ISLEAPYEAR(YEAR, LEAP)

      DLEFT = DOY
      do M = 1, 12
         if (M .eq. 2 .and. LEAP) then
            if (DLEFT .le. 29) then
               MONTH = 2
               DAY   = DLEFT
               return
            else
               DLEFT = DLEFT - 29
            end if
         else
            if (DLEFT .le. NDAYS(M)) then
               MONTH = M
               DAY   = DLEFT
               return
            else
               DLEFT = DLEFT - NDAYS(M)
            end if
         end if
      end do

      MONTH = 12
      DAY   = 31
      return
      end
c     *************************************************************************
      subroutine IRI2HWM07_SETUP
     i         ( YEAR, MONTH, DAY, HOUR, MINUTE, SECOND, LonRef,
     o           IYD, SECUT, HRL)
c     -------------------------------------------------------------------------
c     Purpose:
c     - Take IRI-style UTC calendar inputs (YYYY,MM,DD,HH:MM:SS, lat, lon)
c     - Produce consistent HWM07 inputs:
c         IYD   = YYDDD  (year/day)
c         SECUT = UT seconds since 00:00
c         HRL   = solar local time [0,24) consistent with SECUT and LonRef
c
c      Notes:
c      - HRL here is the longitude-based solar local time:
c            HRL = (UT_hours + Lon/15) mod 24
c      - Longitude is assumed degrees EAST-positive. Any range is accepted.
c      - Lat is passed through (not used by these routines).
c     -------------------------------------------------------------------------
      implicit none
c
c     --- inputs (IRI-style) ---
      integer YEAR, MONTH, DAY, HOUR, MINUTE, SECOND
      real    LonRef
c
c     --- outputs (HWM07-style) ---
      integer IYD
      real    SECUT, HRL
c
c     --- locals ---
      integer DOY, YY
c
c     --- day-of-year ---
      call DAYOFYEAR(YEAR, MONTH, DAY, DOY)
c
c     --- IYD = YYDDD ---
      YY  = MOD(YEAR,100)
      IYD = YY*1000 + DOY
c
c     --- UT seconds since 00:00 ---
      SECUT = 3600.0*REAL(HOUR) + 60.0*REAL(MINUTE)
     &      + REAL(SECOND)
c
c     --- compute HRL consistent with SECUT and LonRef ---
      call UTSEC_TO_HRL(SECUT, LonRef, HRL)
      return
      end
c     ************************************************************************
      subroutine UTSEC_TO_HRL(SECUT, LONDEG, HRL)
c     ------------------------------------------------------------------------
c     Compute solar local time HRL consistent with UT seconds and longitude.
c          HRL = (UT_hours + Lon/15) mod 24, Lon in degrees east.
c     ------------------------------------------------------------------------
      implicit none
      real SECUT, LONDEG, HRL
      real UTH, LONN, TMP
c
c     --- normalize longitude to [0,360) ---
      LONN = MOD(LONDEG,360.0)
      if (LONN .lt. 0.0) LONN = LONN + 360.0
c
c     --- UT in hours ---
      UTH = SECUT / 3600.0
c
c     --- solar local time in hours ---
      TMP = UTH + LONN/15.0
      HRL = MOD(TMP,24.0)
      if (HRL .lt. 0.0) HRL = HRL + 24.0
      return
      end
c     *************************************************************************
      subroutine DAYOFYEAR(YEAR, MONTH, DAY, DOY)
c     -------------------------------------------------------------------------
c              DAYOFYEAR: day-of-year from calendar date (Gregorian).
c     -------------------------------------------------------------------------
      implicit none
      integer YEAR, MONTH, DAY, DOY
      integer M
      integer NDAYS(12)
      logical LEAP

      data NDAYS /31,28,31,30,31,30,31,31,30,31,30,31/

      call ISLEAPYEAR(YEAR, LEAP)

      DOY = DAY
      do M = 1, MONTH - 1
         DOY = DOY + NDAYS(M)
      end do
      if (LEAP .and. MONTH .gt. 2) DOY = DOY + 1
      return
      end
c     *************************************************************************
      subroutine ISLEAPYEAR(YEAR, LEAP)
c     -------------------------------------------------------------------------
c                      ISLEAPYEAR: Gregorian leap-year test.
c     ------------------------------------------------------------------------
      implicit none
      integer YEAR
      logical LEAP

      LEAP = .false.
      if (MOD(YEAR,4) .ne. 0) then
         LEAP = .false.
      else if (MOD(YEAR,100) .ne. 0) then
         LEAP = .true.
      else if (MOD(YEAR,400) .ne. 0) then
         LEAP = .false.
      else
         LEAP = .true.
      end if
      return
      end
