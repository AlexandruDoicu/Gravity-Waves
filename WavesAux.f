C     *************************************************************************
C
c
C                           Auxiliary routines
C           Lapack: ZGEEV0, ZGBTRF, ZGBTRS, ZGETRF, ZGETRI, ZGEMM
C
C
c     *************************************************************************
      real*8 function Cv_vib_molar(T, THETA, Rgas)
c     -------------------------------------------------------------------------
C     Vibrational molar heat capacity (harmonic oscillator, one mode)
C         Cv_vib = R * (x^2 * exp(x)) / (exp(x) - 1)^2, with x = THETA/T
c     -------------------------------------------------------------------------
      implicit none
      real*8 T, THETA, Rgas
      real*8 x, ex, denom
C
      x = THETA / T
C
C     --- guard extreme values ---
      if (x .gt. 700.d0) then
         Cv_vib_molar = 0.d0
         return
      end if

      ex = exp(x)
      denom = ex - 1.d0
      Cv_vib_molar = Rgas * (x*x * ex) / (denom*denom)
      end
C     *************************************************************************
      subroutine csplineD
     i         ( doAbsoluteVal, n, y, lam,
     o           x )
c     ------------------------------------------------------------------------
c                          Spline interpolation
c     ------------------------------------------------------------------------
      implicit none
!
!     --- input arguments ---
      integer n
      real*8  lam
      real*8  y(n) !data to be smoothed
      logical doAbsoluteVal
!
!     --- Output arguments
      real*8 x(n) !smoothed data
!
!     --- Local variables ---
      integer    r, JJ
      real*8     T, score
!
      real*8     a0, a1, d
      complex*16 a2, x1, x2, alpha, beta, atmp
      real*8     fac1, fac2
      real*8     elim, flim, glim, hlim, qlim
      real*8     g1, g2, h
      real*8     j1, j2, j3, j4
      real*8     lamc, lamr, mu, q, sq
      real*8     Tcu
      real*8     tmp1, tmp2, tmp3, tmp4
      real*8     tr, tr1, tr2, tr3
      real*8     v
      integer    ir, nc
      integer    r2, rn, rsq
      integer    i, j, k
      integer    NN, m

      real*8, allocatable :: c(:), e(:), f(:), w(:), z(:)
!
      r  = 1
      JJ = 6
      T  = 0.10d0
!
      nc = CEILING(dble(n)/2.d0)
      Tcu = T**3
      rn = r*n
!
!     --- Initialize output vector 'x' with zero for all 'm' members
      m = r*n+r-1
      x = 0.0d0
!
      allocate(c(n))
      c = 0.0d0
      allocate(w(n-2))
      w = 0.0d0
      allocate(z(n))
      z = 0.0d0

      do j=1,n-2
        w(j) = y(j)-2*y(j+1)+y(j+2)
      end do

      lamr = lam*Tcu
      a0 = 6.d0 + lamr*2.d0/3.d0
      a1 = 4.d0 - lamr/6.d0
      atmp = a1**2 - 4.d0*(a0-2.d0)
      a2 = SQRT(atmp)
      x1 =  (a1 + a2) / 2.d0
      x2 = -(a2 - a1) / 2.d0

      if (lamr .gt. 24.0) then
        alpha = 0.5d0*(x1 + SQRT(x1**2-4.d0))
        beta  = 0.5d0*(x2 + SQRT(x2**2-4.d0))
      else if (lamr .lt. 24.d0) then
        alpha = 0.5d0*(x1 - SQRT(x1**2-4.d0))
        beta  = 0.5d0*(x2 - SQRT(x2**2-4.d0))
      else
        alpha = 0.5d0*(x1 - SQRT(x1**2-4.d0))
        beta  = 0.5d0*(x2 + SQRT(x2**2-4.d0))
      end if
      if (JJ .gt. LOG10(realpart(alpha*beta))
     &          - (nc-1)*2*LOG10(ABS(alpha))) then

!        if (allocateD(e)) deallocate(e)
        allocate(e(n-2))
        e = 0.0d0

!        if (allocateD(f)) deallocate(f)
        allocate(f(n-2))
        f = 0.0d0

        d = a0
        f(1) = 1.d0 / d
        c(2) = f(1)*w(1)
        mu = a1
        e(1) = mu*f(1)
        d = a0-mu*e(1)
        f(2) = 1.d0 / d
        c(3) = f(2)*(w(2)+mu*c(2))
        mu = a1 - e(1)
        e(2) = mu*f(2)

        do j=3,n-2
          d = a0-mu*e(j-1)-f(j-2)
          f(j) = 1.d0 / d
          c(j+1) = f(j)*(w(j)+mu*c(j)-c(j-1))
          mu = a1 - e(j-1)
          e(j) = mu * f(j)
        end do
        c(n-2) = c(n-2)+e(n-3)*c(n-1)

        do j=n-4,1,-1
          c(j+1) = c(j+1) + e(j)*c(j+2) - f(j)*c(j+3)
        end do

        g2 = f(n-2)
        tr1 = g2
        h = e(n-3)*g2
        tr2 = h
        g1 = f(n-3) + e(n-3)*h
        tr1 = tr1 + g1
        tr3=0

        do k=n-4,n-nc,-1
          q = e(k)*h - f(k)*g2
          tr3 = tr3 + q
          h = e(k)*g1 - f(k)*h
          tr2 = tr2 + h
          g2 = g1
          g1 = f(k)*(1-q) + e(k)*h
          tr1 = tr1 + g1
        end do

        q = e(n-nc-1)*h - f(n-nc-1)*g2
        tr3 = tr3 + q
        h = e(n-nc-1)*g1 - f(n-nc-1)*h
        tr2 = tr2 + h

        tr1 =  6.d0*(2.d0*tr1 - dble(MOD(n,2))*g1)
        tr2 = -8.d0*(2.d0*tr2 - (1 + dble(MOD(n,2)) )*h)
        tr3 =  2.d0*(2.d0*tr3 - dble(MOD(n,2))*q)
        tr = (tr1+tr2+tr3)/dble(n)

      else
        flim = dble( alpha*beta )
        elim = dble( alpha + beta )
        glim = flim*(1.d0+flim)
     &       / ((1.d0 - flim)*((1.d0+flim)**2 - elim**2))
        hlim = elim*glim/(1.d0+flim)
        qlim = elim*hlim - flim*glim
        NN = CEILING((LOG10(flim) - JJ) / (2*LOG10(ABS(alpha))))

        if (allocateD(e)) deallocate(e)
        allocate(e(NN))
        e = 0.0d0
        if (allocateD(f)) deallocate(f)
        allocate(f(NN))
        f = 0.0d0
        d = a0
        f(1) = 1.d0 / d
        c(2) = f(1)*w(1)
        mu = a1
        e(1) = mu*f(1)
        d = a0 - mu*e(1)
        f(2) = 1.d0 / d
        c(3) = f(2)*(w(2) + mu*c(2))
        mu = a1 - e(1)
        e(2) = mu*f(2)
        g2 = flim
        tr1 = g2
        h = elim*g2
        tr2 = h
        g1 = flim+elim*h
        tr1 = tr1+g1
        tr3 = 0.d0
        do j=3,NN
          d = a0 - mu*e(j-1) - f(j-2)
          f(j) = 1.d0 / d
          c(j+1) = f(j)*(w(j) + mu*c(j) - c(j-1))
          mu = a1 - e(j-1)
          e(j) = mu*f(j)
          q = elim*h - flim*g2
          tr3 = tr3 + q
          h = elim*g1 - flim*h
          tr2 = tr2+h
          g2 = g1
          g1 = flim*(1-q) + elim*h
          tr1 = tr1 + g1
        end do
        tr1 = tr1+(nc-NN-1)*glim
        tr2 = tr2+(nc-NN)*hlim
        tr3 = tr3+(nc-NN)*qlim
        tr1 =  6.d0*(2.d0*tr1 - dble(MOD(n,2))*glim)
        tr2 = -8.d0*(2.d0*tr2 - (1.0 + dble(MOD(n,2)) )*hlim)
        tr3 =  2.d0*(2.d0*tr3 - dble(MOD(n,2))*qlim)
        tr = (tr1 + tr2 + tr3)/dble(n)
        mu = a1 - elim
        do j=NN+1,n-2
          c(j+1) = flim*(w(j) + mu*c(j)-c(j-1))
        end do
        c(n-2) = c(n-2) + elim*c(n-1)

        do j=n-3,NN+2,-1
          c(j) = c(j) + elim*c(j+1) - flim*c(j+2)
        end do

        do j=NN,1,-1
          c(j+1) = c(j+1) + e(j)*c(j+2) - f(j)*c(j+3)
        end do
      end if

      z(1) = c(2)
      z(2) = c(3)-2*c(2)
      do j=3,n-2
        z(j) = c(j-1)-2*c(j)+c(j+1)
      end do

      z(n-1) = c(n-2)-2*c(n-1)
      z(n) = c(n-1)

      sq = sum (z * z)
      score = sq / tr**2

!     x(r:rn:r) = y-z
      do i = 1, n
        x(i) = y(i) - z(i)
      end do
!
      if (r .lt. 8) then
        fac1 = x(2*r) - x(r)    - lamr*c(2)/6.d0
        fac2 = x(rn)  - x(rn-r) + lamr*c(n-1)/6.d0
        do j=1,r-1
          j1 = dble(j)/dble(r)
          j2 = 1.d0 - j1
          v = lamr*j1*j2 / 6.d0
          j3 = v*(1.d0 + j1)
          j4 = v*(2.d0 - j1)
          do i=1,n-1
            ir = i*r
            x(ir+j) = j2*x(ir) + j1*x(ir+r) - j3*c(i+1) - j4*c(i)
          end do
          x(j)    = x(r)  - j2*fac1
          x(rn+j) = x(rn) + j1*fac2
        end do

      else
        lamc = lamr/(6*r**3)
        r2 = 2*r
        rsq = r**2
        do i=1,n-1
          ir = i*r
          tmp1 = x(ir)   / r
          tmp2 = x(ir+r) / r
          tmp3 = lamc*c(i+1)
          tmp4 = lamc*c(i)
          do j=1,r-1
            x(ir+j) = dble(r-j)*tmp1 + dble(j)*tmp2
     &              - dble(j)*dble(rsq-j*j)*tmp3
     &              - dble(j)*dble(r-j)*dble(r2-j)*tmp4
          end do
        end do
        tmp1 = x(r) - x(r+1)
        tmp2 = x(rn) - x(rn-1)
        do j = 1,r-1
          x(j)    = x(r) + (r-j)*tmp1
          x(rn+j) = x(rn) + j*tmp2
        end do
      end if
      if ( doAbsoluteVal ) then
        do i = 1, n
           if ( x(i) .lt. 0.d0 ) x(i) = abs(x(i))
        end do
      end if
!
      if (allocated(c)) deallocate(c)
      if (allocated(e)) deallocate(e)
      if (allocated(f)) deallocate(f)
      if (allocated(w)) deallocate(w)
      if (allocated(z)) deallocate(z)
      end
!     *******************************************************************
      subroutine SortEigenvalues( N, NMAX, RA, RB, RC )
      implicit none
      integer    N, NMAX
      real*8     RA( NMAX )
      complex*16 RB( NMAX )
      complex*16 RC( NMAX, NMAX )
C
      integer    L, IR, I, J, k
      real*8     RRA
      complex*16 RRB
      complex*16, allocatable:: RRC(:)
      logical  MORE
C
      allocate( RRC(N) )
C
      L  = N / 2 + 1
      IR = N
      MORE = .true.
      do while ( MORE )
        if ( L .GT. 1 ) then
          L   = L - 1
          RRA = RA(L)
          RRB = RB(L)
          do k = 1, N
            RRC(k) = RC(k,L)
          end do
        else
          RRA = RA(IR)
          RRB = RB(IR)
          do k = 1, N
            RRC(k) = RC(k,IR)
          end do
C
          RA(IR) = RA(1)
          RB(IR) = RB(1)
          do k = 1, N
            RC(k,IR) = RC(k,1)
          end do
C
          IR = IR -1
          if ( IR .EQ. 1 ) then
            RA(1) = RRA
            RB(1) = RRB
            do k = 1, N
              RC(k,1) = RRC(k)
            end do
            deallocate( RRC )
            return
          end if
        end if
        I = L
        J = L + L
        do while ( J .LE. IR )
          if (J .LT. IR) then
            if ( RA(J) .LT. RA(J+1) ) J = J + 1
          end if
          if ( RRA .LT. RA(J) ) then
            RA(I) = RA(J)
            RB(I) = RB(J)
            do k = 1, N
              RC(k,I) = RC(k,J)
            end do
C
            I = J
            J = J + J
          else
            J = IR + 1
          end if
        end do
        RA(I) = RRA
        RB(I) = RRB
        do k = 1, N
          RC(k,I) = RRC(k)
        end do
      end do
C
      deallocate( RRC )
      end
C     *************************************************************************
      subroutine SortEigenValEigenVct( N, NMAX, RA, RB, RD, RC )
      implicit none
      integer    N, NMAX
      real*8     RA( NMAX )
      complex*16 RB( NMAX ), RD( NMAX )
      complex*16 RC( NMAX, NMAX )
C
      integer    L, IR, I, J, k
      real*8     RRA
      complex*16 RRB, RRD
      complex*16, allocatable:: RRC(:)
      logical  MORE
C
      allocate( RRC(N) )
C
      L  = N / 2 + 1
      IR = N
      MORE = .true.
      do while ( MORE )
        if ( L .GT. 1 ) then
          L   = L - 1
          RRA = RA(L)
          RRB = RB(L)
          RRD = RD(L)
          do k = 1, N
            RRC(k) = RC(k,L)
          end do
        else
          RRA = RA(IR)
          RRB = RB(IR)
          RRD = RD(IR)
          do k = 1, N
            RRC(k) = RC(k,IR)
          end do
C
          RA(IR) = RA(1)
          RB(IR) = RB(1)
          RD(IR) = RD(1)
          do k = 1, N
            RC(k,IR) = RC(k,1)
          end do
C
          IR = IR -1
          if ( IR .EQ. 1 ) then
            RA(1) = RRA
            RB(1) = RRB
            RD(1) = RRD
            do k = 1, N
              RC(k,1) = RRC(k)
            end do
            deallocate( RRC )
            return
          end if
        end if
        I = L
        J = L + L
        do while ( J .LE. IR )
          if (J .LT. IR) then
            if ( RA(J) .LT. RA(J+1) ) J = J + 1
          end if
          if ( RRA .LT. RA(J) ) then
            RA(I) = RA(J)
            RB(I) = RB(J)
            RD(I) = RD(J)
            do k = 1, N
              RC(k,I) = RC(k,J)
            end do
C
            I = J
            J = J + J
          else
            J = IR + 1
          end if
        end do
        RA(I) = RRA
        RB(I) = RRB
        RD(I) = RRD
        do k = 1, N
          RC(k,I) = RRC(k)
        end do
      end do
C
      deallocate( RRC )
      end
C     *************************************************************************
      subroutine SortOneArrayPlusFiveArrays
     &         ( N, NMAX, RA, RB, RC, RD, RE, RF  )
      implicit none
      integer N, NMAX
      real    RA( NMAX ), RB( NMAX ), RC( NMAX ), RD( NMAX ),
     &        RE( NMAX ), RF( NMAX )
C
      integer  L, IR, I, J
      real     RRA, RRB, RRC, RRD, RRE, RRF
      logical  MORE
C
      L  = N / 2 + 1
      IR = N
      MORE = .true.
      do while ( MORE )
        if ( L .GT. 1 ) then
          L   = L - 1
          RRA = RA(L)
          RRB = RB(L)
          RRC = RC(L)
          RRD = RD(L)
          RRE = RE(L)
          RRF = RF(L)
        else
          RRA = RA(IR)
          RRB = RB(IR)
          RRC = RC(IR)
          RRD = RD(IR)
          RRE = RE(IR)
          RRF = RF(IR)
C
          RA(IR) = RA(1)
          RB(IR) = RB(1)
          RC(IR) = RC(1)
          RD(IR) = RD(1)
          RE(IR) = RE(1)
          RF(IR) = RF(1)
C
          IR = IR -1
          if ( IR .EQ. 1 ) then
            RA(1) = RRA
            RB(1) = RRB
            RC(1) = RRC
            RD(1) = RRD
            RE(1) = RRE
            RF(1) = RRF
            return
          end if
        end if
        I = L
        J = L + L
        do while ( J .LE. IR )
          if (J .LT. IR) then
            if ( RA(J) .LT. RA(J+1) ) J = J + 1
          end if
          if ( RRA .LT. RA(J) ) then
            RA(I) = RA(J)
            RB(I) = RB(J)
            RC(I) = RC(J)
            RD(I) = RD(J)
            RE(I) = RE(J)
            RF(I) = RF(J)
C
            I = J
            J = J + J
          else
            J = IR + 1
          end if
        end do
        RA(I) = RRA
        RB(I) = RRB
        RC(I) = RRC
        RD(I) = RRD
        RE(I) = RRE
        RF(I) = RRF
      end do
      end
C     *************************************************************************
      subroutine SortOneArrayPlusTwoArrays( N, NMAX, RA, RB, RC  )
      implicit none
      integer   N, NMAX
      real*8    RA( NMAX ), RB( NMAX ), RC( NMAX )
C
      integer    L, IR, I, J
      real*8     RRA, RRB, RRC
      logical    MORE
C
      L  = N / 2 + 1
      IR = N
      MORE = .true.
      do while ( MORE )
        if ( L .GT. 1 ) then
          L   = L - 1
          RRA = RA(L)
          RRB = RB(L)
          RRC = RC(L)
        else
          RRA = RA(IR)
          RRB = RB(IR)
          RRC = RC(IR)
C
          RA(IR) = RA(1)
          RB(IR) = RB(1)
          RC(IR) = RC(1)
C
          IR = IR -1
          if ( IR .EQ. 1 ) then
            RA(1) = RRA
            RB(1) = RRB
            RC(1) = RRC
            return
          end if
        end if
        I = L
        J = L + L
        do while ( J .LE. IR )
          if (J .LT. IR) then
            if ( RA(J) .LT. RA(J+1) ) J = J + 1
          end if
          if ( RRA .LT. RA(J) ) then
            RA(I) = RA(J)
            RB(I) = RB(J)
            RC(I) = RC(J)
C
            I = J
            J = J + J
          else
            J = IR + 1
          end if
        end do
        RA(I) = RRA
        RB(I) = RRB
        RC(I) = RRC
      end do
      end
C     ****************************************************************
      subroutine SortOneArrayPlusTwoMatrices
     &         ( N, M, NMAX, RA, RB, RC )
      implicit none
      integer  N, M, NMAX
      real*8   RA(NMAX), RB(M,NMAX), RC(M,NMAX)
C
      integer  L, IR, I, J, K
      real*8   RRA
      logical  MORE
      real*8, allocatable:: RRB( : ), RRC( : )
C
      allocate( RRB(M), RRC(M) )
      L  = N / 2 + 1
      IR = N
      MORE = .TRUE.
      do while ( MORE )
        if ( L .GT. 1 ) then
          L   = L - 1
          RRA = RA(L)
          do K = 1, M
            RRB(K) = RB(K,L)
            RRC(K) = RC(K,L)
          end do
        else
          RRA = RA(IR)
          do K = 1, M
            RRB(K) = RB(K,IR)
            RRC(K) = RC(K,IR)
          end do

          RA(IR) = RA(1)
          do K = 1, M
            RB(K,IR) = RB(K,1)
            RC(K,IR) = RC(K,1)
          end do

          IR = IR -1
          if ( IR .EQ. 1 ) then
            RA(1) = RRA
            do K = 1, M
              RB(K,1) = RRB(K)
              RC(K,1) = RRC(K)
            end do
            deallocate( RRB, RRC )
           return
         end if
       end if
       I = L
       J = L + L
       do while ( J .LE. IR )
         if (J .LT. IR) then
           if ( RA(J) .LT. RA(J+1) ) J = J + 1
         end if
         if ( RRA .LT. RA(J) ) then
           RA(I) = RA(J)
           do K = 1, M
             RB(K,I) = RB(K,J)
             RC(K,I) = RC(K,J)
           end do

           I = J
           J = J + J
         else
           J = IR + 1
         end if
       end do
       RA(I) = RRA
       do K = 1, M
         RB(K,I) = RRB(K)
         RC(K,I) = RRC(K)
       end do
      end do
C
      deallocate( RRB, RRC )
      end
