!    _______              __
!   / ___/ /  ___  __ _  / /  ___
!  / /__/ _ \/ _ \/  V \/ _ \/ _ \
!  \___/_//_/\___/_/_/_/_.__/\___/
!
!  CHOMBO Copyright (c) 2000-2004, The Regents of the University of
!  California, through Lawrence Berkeley National Laboratory (subject to
!  receipt of any required approvals from U.S. Dept. of Energy).  All
!  rights reserved.
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions are met:
!
!  (1) Redistributions of source code must retain the above copyright
!  notice, this list of conditions and the following disclaimer.
!  (2) Redistributions in binary form must reproduce the above copyright
!  notice, this list of conditions and the following disclaimer in the
!  documentation and/or other materials provided with the distribution.
!  (3) Neither the name of Lawrence Berkeley National Laboratory, U.S.
!  Dept. of Energy nor the names of its contributors may be used to endorse
!  or promote products derived from this software without specific prior
!  written permission.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
!  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
!  PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!  OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!  You are under no obligation whatsoever to provide any bug fixes,
!  patches, or upgrades to the features, functionality or performance of
!  the source code ("Enhancements") to anyone; however, if you choose to
!  make your Enhancements available either publicly, or directly to
!  Lawrence Berkeley National Laboratory, without imposing a separate
!  written license agreement for such Enhancements, then you hereby grant
!  the following license: a non-exclusive, royalty-free perpetual license
!  to install, use, modify, prepare derivative works, incorporate into
!  other computer software, distribute, and sublicense such Enhancements or
!  derivative works thereof, in binary and source code form.
!
!  TRADEMARKS. Product and company names mentioned herein may be the
!  trademarks of their respective owners.  Any rights not expressly granted
!  herein are reserved.

      subroutine POISSON3D( f ,N1,N2,N3 ,H ,Op )

! Node-centered Poisson solver with homogenous Dirichlet boundary conditions.
!
! inputs:
! N1,N2,N3 : number of cells in each dimension; should be highly
!     composite, must be even
! f(0:N1,0:N2,0:N3): rhs on input, solution to Poisson's equation on
!                 output. Rhs is defined on interior nodes [1:N*-1]^3.
! H : mesh spacing.
! Op: operator type: 0 for spectral, 1 for 19-point Mehrstellen (4th order),
!         2 for 7-point (2nd order), 3 for 27-point Mehrstellen (6th order))
!
! The algorithm uses SFFT3D(), a 3d sine transform that uses
! the real-to-complex 1D FFT's in the FFTW library.
!
! The structure of the program is as follows:
! (1) Compute the 3d sine transform of the inputs - SFFT3D()
!     updates f in place.
! (2) Multiply each Fourier component by the symbol of the operator
!     evaluated at that component. These are easily derived using
!     the addition formulas for sin(A+B). In the case of the spectral
!     approximation, multiply the Fourier coefficients by the symbol
!     of the Laplacian. We also scale the coefficients to account for
!     the fact that the sine transform applied twice is (n/2)^3 times
!     the identity operator.
! (3) Apply the sine transform again, which is its own inverse.

      implicit none
!Args
      integer N1,N2,N3 ,Op
      real*8  f(0:N1,0:N2,0:N3)
      real*8  H(3)
!Locals
      integer nmax
      ! Need the max dimension to declare some local vars so
      ! use automatic arrays in a separate subroutine because
      ! g77 does not support allocatable.
      nmax = MAX( N1,N2,N3 )
      call POISSON3D_1( f ,N1,N2,N3 ,H,Op ,nmax )
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine POISSON3D_1( f ,N1,N2,N3 ,H,Op ,Nmax )
      implicit none
!Constants
      real*8     One   ,Two   ,Three   ,Four   ,TwoThirds
      parameter( One=1 ,Two=2 ,Three=3 ,Four=4 ,TwoThirds=Two/Three )
!Args
      integer N1,N2,N3 ,Op ,Nmax
      real*8  f(0:N1,0:N2,0:N3)
      real*8  H(3)
!Locals
      integer d,i,j,k,ll ,n(3) ,nh(3)
      logical isotropicgrid
      real*8 coscoef(Nmax-1,3) ,Pi ,ampl
      real*8 hsqrd(3) ,over18hsqrd(3) ,pioverll(3)
     &      ,overhsqrd ,over30hsqrd ,nh3
!XXX      real*8 coefk ,coefjplusk ,coefjtimesk
      real*8 coefjplusk ,coefjtimesk
      integer itest
      save    itest
      data    itest/10/

      ! The number of grid cells must be even.
      if ((N1/2)*2.NE.N1 .OR. (N2/2)*2.NE.N2 .OR. (N3/2)*2.NE.N3 ) then
         write(6,*) 'error: POISSON3D: N1,N2,N3 must be even, not '
     &      ,N1,N2,N3
         STOP 'POIS3a'
      endif

      ! Compute some setup information.
      Pi = ATAN( One ) * 4
      n (1) = N1
      n (2) = N2
      n (3) = N3
      nh(1) = N1/2
      nh(2) = N2/2
      nh(3) = N3/2
      nh3 = nh(1) * nh(2) * nh(3)
      do d = 1,3
        hsqrd(d) = H(d) * H(d)
        do ll = 1,n(d)-1
            coscoef(ll,d) = COS( ll * pi / n(d) )
        enddo
      enddo
      IsotropicGrid = ABS( H(1) / H(2) - One ) .LT. 1d-6 .AND.
     &                ABS( H(1) / H(3) - One ) .LT. 1d-6

      ! Apply the sine transform.
      call SFFT3D( f ,N1,N2,N3 )

!TEST
      if( itest .EQ. 0 )then
        print *,'dbg: 3dPOISSON: 1st SFFT3d, writing f to fft1.dat'
        call WRITE3D( 'fft1.dat' ,N1,N2,N3,0 ,H ,f )
        itest = itest + 1
      endif
!TEST

      ! Multiply each Fourier coefficient by the symbol of the operator
      ! and divide by the factor of Volume/8 that comes from performing
      ! the sine transform twice

      if( op .EQ. 0 )then  ! spectral discretization
        pioverll(1) = ( Pi * Pi ) / ( n(1) * n(1) * hsqrd(1) )
        pioverll(2) = ( Pi * Pi ) / ( n(2) * n(2) * hsqrd(2) )
        pioverll(3) = ( Pi * Pi ) / ( n(3) * n(3) * hsqrd(3) )
        do k = 1,N3-1
          do j = 1,N2-1
            do i = 1,N1-1
              ampl = -( i*i*pioverll(1) + j*j*pioverll(2)
     &                + k*k*pioverll(3) )
              f(i,j,k) = f(i,j,k) / ( ampl * nh3 )
            enddo
          enddo
        enddo

      elseif( Op .EQ. 1 .AND. .NOT. IsotropicGrid )then  !Mehrstellen 19pt, non-uniform grid
        ! Let:
        !  u   = u(i,j,k)
        !  ui  = u(i-1,j,k) + u(i+1,j,k)
        !  uij = u(i-1,j-1,k) + u(i+1,j-1,k) + u(i-1,j+1,k) + u(i+1,j+1,k)
        ! The stencil for anisotropic grids is:
        !  uxx + uyy + uzz =
        !  ( 8ui - 4uj + 2uk + 4uij - 2ujk +  uik - 24u )/( 18 H(1)^2 ) +
        !  ( 2ui + 8uj - 4uk +  uij + 4ujk - 2uik - 24u )/( 18 H(2)^2 ) +
        !  (-4ui + 2uj + 8uk - 2uij +  ujk + 4uik - 24u )/( 18 H(3)^2 )

        ! The discrete sine transform (DST) manipulation of this expression is long
        ! and boring.  It just repeats the technique used with the 7pt stencil for
        ! all the pairs of terms, with the added complexity that the off-axis terms
        ! have 2 cosine factors (one for each dimension that is off-axis).  To get
        ! the proper cancellation you have to apply the DST to all 4 stencil points
        ! in the off-axis plane.  So for the first line of 1/18h^2 terms we get:
        !    8*2cos(pi*hx*l) - 4*2cos(pi*hy*m) + 2*2cos(pi*hz*n)
        !  + 4*4cos(pi*hx*l)cos(pi*hy*m) - 2*4cos(pi*hy*m)cos(pi*hz*n)
        !  +   4cos(pi*hx*l)cos(pi*hz*n)
        ! and similarly for the others.
        over18hsqrd(1) = One / ( H(1)**2 * 18 )
        over18hsqrd(2) = One / ( H(2)**2 * 18 )
        over18hsqrd(3) = One / ( H(3)**2 * 18 )
        do k = 1,N3-1
          do j = 1,N2-1
            do i = 1,N1-1
              ampl = ( 16*coscoef(i,1) - 8*coscoef(j,2) + 4*coscoef(k,3)
     &               + 16 * coscoef(i,1) * coscoef(j,2)
     &               -  8 * coscoef(j,2) * coscoef(k,3)
     &               +  4 * coscoef(i,1) * coscoef(k,3) - 24
     &               ) * over18hsqrd(1)
     &             + ( 4*coscoef(i,1) + 16*coscoef(j,2) - 8*coscoef(k,3)
     &               +  4 * coscoef(i,1) * coscoef(j,2)
     &               + 16 * coscoef(j,2) * coscoef(k,3)
     &               -  8 * coscoef(i,1) * coscoef(k,3) - 24
     &               ) * over18hsqrd(2)
     &             + (-8*coscoef(i,1) + 4*coscoef(j,2) + 16*coscoef(k,3)
     &               -  8 * coscoef(i,1) * coscoef(j,2)
     &               +  4 * coscoef(j,2) * coscoef(k,3)
     &               + 16 * coscoef(i,1) * coscoef(k,3) - 24
     &               ) * over18hsqrd(3)
              f(i,j,k) = f(i,j,k) / ( ampl * nh3 )
            enddo
          enddo
        enddo

      elseif( Op .EQ. 1 .AND. IsotropicGrid )then  ! Mehrstellen operator, 19 point, O(h^4)
        ! Let H(1)==H(2)==H(3) and simplify the formula for anisotropic grid above:
        ! uxx + uyy + uzz =
        ! (8+2-4)ui +(-4+8+2)uj +(2-4+8)uk +(+4+1-2)uij +(-2+4+1)ujk +(1-2+4)uik -72u / 18h^2
        !==>  6(ui+uj+uk) + 3(uij+ujk+uik) - 72u / 18h^2
        !==>  [ 1/3(ui+uj+uk) + 1/6(uij+ujk+uik) - 4 ] / h^2
        ! Applying DST3() to the above multiplies the ui,uj,uk terms by 2 and the
        ! uij,ujk,uik terms by 4, so we get:
        ! [ 2/3 ( cos(i) + cos(j) + cos(k) )
        ! + 2/3 ( cos(i)cos(j) + cos(j)cos(k) + cos(i)cos(k) ) -4 ] / h^2
        overhsqrd = One / H(1)**2  !all H are equal
        do k = 1,N3-1
          do j = 1,N2-1
            coefjtimesk = coscoef(j,2)*coscoef(k,3)
            coefjplusk = coscoef(j,2)+coscoef(k,3)
            do i = 1,N1-1
              ampl = ( TwoThirds
     &               * ( coscoef(i,1)*coscoef(j,2)
     &                 + coefjtimesk
     &                 + coscoef(i,1)*coscoef(k,3)
     &                 + coscoef(i,1)
     &                 + coefjplusk
     &                 ) - Four
     &               ) * overhsqrd
              f(i,j,k) = f(i,j,k) / ( ampl * nh3 )
            enddo
          enddo
        enddo

      elseif( op .EQ. 2 )then   ! standard 7-point Laplacian operator
        ! the anisotropic stencil for uxx + uyy + uzz is:
        ! ( u(i-1,j,k) + u(i+1,j,k) - 2u(i,j,k) )/h(1)^2 +
        ! ( u(i,j-1,k) + u(i,j+1,k) - 2u(i,j,k) )/h(2)^2 +
        ! ( u(i,j,k-1) + u(i,j,k+1) - 2u(i,j,k) )/h(3)^2
        ! After applying the sin transform (DST) to the right hand side (f)
        ! we have, approximately, the DST of the stencil.
        ! Applying the 3d DST to the stencil gives:
        !    DST( u(i-1) + u(i+1) - 2u )/h(1)^2 + <for j,k dirs>
        ! ==>DST( u(i*hx-hx) + u(i*hx+hx) - 2u )/h(1)^2 + ...
        ! applying DST and using sin(a+b) identity to isolate u(i*hx) and
        ! let ci,cj,ck=appropriate constants, l,m,n=indices in transform
        ! = sum_{lmn}(sin(ci*i*l) cos(ci*l) + cos(-ci*i*l) sin(-ci*l)
        !           + sin(ci*i*l) cos(ci*l) + cos(+ci*i*l) sin(+ci*l)
        !           - 2 sin(ci*i*l) )/h1^2 + ...
        ! the 2nd term cancels because sin(-x)=-sin(x), giving
        ! = sum_{lmn} 2 sin(ci*i*l) cos(ci*l) - 2 sin(ci*i*l)/h1^2
        ! = sum_{lmn} 2 sin(ci*i*l) ( cos(ci*l) - 1 )/h1^2
        ! which is equivalent to DST( u*2*( cos(ci*l) - 1 )/h1^2 )
        ! so we can divide f==DST(rhs) by 2(cos(ci*l)-1)/h1^2 and apply the DST again,
        ! which inverts the transform, and scale by the number of points,
        ! to get the solution to the Poisson equation: u.
!XXX -- unoptimized version
        do k = 1,N3-1
          do j = 1,N2-1
            do i = 1,N1-1
              ampl = ( Two * coscoef(i,1) - Two ) / H(1)**2
     &             + ( Two * coscoef(j,2) - Two ) / H(2)**2
     &             + ( Two * coscoef(k,3) - Two ) / H(3)**2
              f(i,j,k) = f(i,j,k) / ( ampl * nh3 )
            enddo
          enddo
        enddo
!XXX -- optimized version
!XXX        do k = 1,N3-1
!XXX          coefk = Two * ( coscoef(k,3) - One )/hsqrd(3)
!XXX          do j = 1,N2-1
!XXX            coefjplusk = Two * ( coscoef(j,2) - One )/hsqrd(2) + coefk
!XXX            do i = 1,N1-1
!XXX              ampl = Two * ( coscoef(i,1) - One )/hsqrd(1) + coefjplusk
!XXX              f(i,j,k) = f(i,j,k) / ( ampl * nh3 )
!XXX            enddo
!XXX          enddo
!XXX        enddo

      elseif( Op .EQ. 3 .AND. IsotropicGrid )then ! Mehrstellen operator, 27 point, O(h^6)
        over30hsqrd = One / ( H(1) * H(1) * 30 )
!XXX -- unoptimized version from Phil Colella
!XXX                        ampl = ( coscoef(i)*coscoef(j)*coscoef(k)*8
!XXX     &                         + 12*( coscoef(i)*coscoef(j)
!XXX     &                              + coscoef(j)*coscoef(k)
!XXX     &                              + coscoef(i)*coscoef(k) )
!XXX     &                         + 28*( coscoef(i)+coscoef(j)+coscoef(k) )
!XXX     &                         - 128.d0 ) * Over30H2
        do k = 1,N3-1
          do j = 1,N2-1
            coefjtimesk = coscoef(j,2) * coscoef(k,3)
            coefjplusk  = coscoef(j,2) + coscoef(k,3)
            do i = 1,N1-1
              ampl = (  8 * coscoef(i,1) * coefjtimesk
     &               + 12 * ( coscoef(i,1) * coscoef(j,2)
     &                      + coefjtimesk
     &                      + coscoef(i,1) * coscoef(k,3) )
     &               + 28 * ( coscoef(i,1) + coefjplusk )
     &               - 128 ) * over30hsqrd
              f(i,j,k) = f(i,j,k) / ( ampl * nh3 )
            enddo
          enddo
        enddo

      elseif( Op .EQ. 3 .AND. .NOT. IsotropicGrid )then ! Mehrstellen operator, 27 point, O(h^6)
        write(6,*) 'error: POISSON3D(): 6th order Mehr operator not '
     &            ,'implemented for non-isotropic grids.'
        stop 'POIS3c'
      else
        write(6,*) 'error: POISSON3D(): unknown Op argument: ',Op
        stop 'POIS3b'
      endif

!TEST
      if( itest .EQ. 1 )then
        print *,'dbg: 3dPOISSON: b4 2nd SFFT3d, writing f to fft2.dat'
        call WRITE3D( 'fft2.dat' ,N1,N2,N3,0 ,H ,f )
        itest = itest + 2
      endif
!TEST
      ! Compute the inverse sine transform.
      call SFFT3D( f ,N1,N2,N3 )

!TEST
      if( itest .EQ. 3 )then
        print *,'dbg: 3dPOISSON: 2nd SFFT3d, writing f to fft3.dat'
        call WRITE3D( 'fft3.dat' ,N1,N2,N3,0 ,H ,f )
        itest = itest + 1
      endif
!TEST
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine WRITE3D( Fn ,N1,N2,N3 ,S ,H ,D3d )
      implicit none
      character Fn*(*)
      integer N1,N2,N3 ,S
      real*8  H(3) ,D3d(-S:N1+S,-S:N2+S,-S:N3+S)
      integer i,j,k
      open(1,file=Fn)
      do k = -S,N3+S
        do j = -S,N2+S
          do i = -S,N1+S
            write(1,1) i,j,k,H(1)*i,H(2)*j,H(3)*k,D3d(i,j,k)
          enddo
        enddo
      enddo
      close(1)
      return
    1 format(3(1x,I4),3(1x,F9.6),1P,1X,E18.11)
      end

