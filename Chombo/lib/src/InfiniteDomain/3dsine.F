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

! Purpose:
!   Compute a 3D sine transform on a rectangular grid using a series
!   of 1D sine transforms with FFTW.
!
! Notes:
!   The grid must have an even number of cells in each dimension; they should be
!     highly composite.
!   It only works for double precision.
!   The cpp macro variable FFTW_PLAN_TYPE must be defined to an integer declaration
!     of the same length as the default C pointer type.
!
! Modification History:
!   Originally written by Peter McCourquodale for cubic grids
!   FFTW version 3 written by Dan Martin.
!   Rectangular grid and 64bit support by David Serafini, Jun2004.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#include "REAL.H"
!#include "SPACE.H"
#include "CONSTANTS.H"

      subroutine SFFT3D( f ,N1,N2,N3 )
      implicit none
      integer N1,N2,N3
      REAL_T f(0:N1,0:N2,0:N3)
      integer nmax
      ! Need the max dimension to declare some local vars so
      ! use automatic arrays in a separate subroutine because
      ! g77 does not support allocatable.
      nmax = MAX( N1,N2,N3 )
#ifdef FFTW_3
      call SFFT3D_1( f ,N1,N2,N3 ,nmax )
#else
      ! Need FFTW!
      call MAYDAY_ERROR()
#endif
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine SFFT3D_1( f ,N1,N2,N3 ,NMax )
      implicit none
#ifdef FFTW_2
#include "fftw_f77.i"
#elif defined( FFTW_3 )
#include "fftw3.f"
#endif

      integer N1,N2,N3 ,NMax
      REAL_T f(0:N1,0:N2,0:N3)

      REAL_T in1d(0:NMax) ,out1d(0:NMax)

      integer d,i,j,k ,n(3)
#ifdef FFTW_2
      integer ll
      REAL_T pi_over_n, sincoef(0:NMax-1,3)
#endif

#ifdef FFTW_3
      ! this must be the same length as a pointer
      FFTW_PLAN_TYPE plan(3)
#endif

      ! N must be an even number.
      if ((N1/2)*2.NE.N1 .OR. (N2/2)*2.NE.N2 .OR. (N3/2)*2.NE.N3 ) then
         write(6,*) 'error: SFFT3D: N1,N2,N3 must be even, not '
     &      ,N1,N2,N3
         stop 'SFFT3D_1a'
      endif

      ! set up constants
      n(1) = N1
      n(2) = N2
      n(3) = N3

#ifdef FFTW_2
      do d = 1,3
        call RFFTW_F77_CREATE_PLAN( plan(d) ,n(d) ,FFTW_REAL_TO_COMPLEX
     &                             ,FFTW_ESTIMATE )
        pi_over_n = Pi/n(d)
        do ll = 0,n(d)-1
          sincoef(ll,d) = SIN( ll * pi_over_n )
        enddo
      enddo

      ! sine transform in first direction
      d = 1
      do k = 0,N3-1
        do j = 0,N2-1
          in1d(0) = 0
          do i = 1,N1-1
            in1d(i) = (f(i,j,k) + f(N1-i,j,k)) * sincoef(i,d)
     &              + (f(i,j,k) - f(N1-i,j,k)) * 0.5
          enddo

          call RFFTW_F77_ONE( plan(d) ,in1d ,out1d )

          f(1,j,k) = out1d(0) * 0.5
          f(0,j,k) = 0

          do i = 1,N1/2-1
            f(2*i  ,j,k) = -out1d(N1-i)
            f(2*i+1,j,k) =  out1d(i)    + f(2*i-1,j,k)
          enddo
        enddo
      enddo

      ! sine transform in the second direction
      d = 2
      do k = 0,N3-1
        do i = 0,N1-1
          in1d(0) = 0
          do j = 1,N2-1
            in1d(j) = ( f(i,j,k) + f(i,N2-j,k) ) * sincoef(j,d)
     &              + ( f(i,j,k) - f(i,N2-j,k) ) * 0.5
          enddo

          call RFFTW_F77_ONE( plan(d) ,in1d ,out1d )

          f(i,1,k) = out1d(0) * 0.5
          f(i,0,k) = 0

          do j = 1,N2/2-1
            f(i,2*j  ,k) = -out1d(N2-j)
            f(i,2*j+1,k) =  out1d(j)    + f(i,2*j-1,k)
          enddo
        enddo
      enddo

      ! sine transform in the third direction
      d = 3
      do j = 0,N2-1
        do i = 0,N1-1
          in1d(0) = 0
          do k = 1,N3-1
            in1d(k) = ( f(i,j,k) + f(i,j,N3-k) ) * sincoef(k,d)
     &              + ( f(i,j,k) - f(i,j,N3-k) ) * 0.5
          enddo

          call RFFTW_F77_ONE( plan(d) ,in1d ,out1d )

          f(i,j,1) = out1d(0) * 0.5
          f(i,j,0) = 0

          do k = 1,N3/2-1
            f(i,j,2*k)   = -out1d(N3-k)
            f(i,j,2*k+1) =  out1d(k)    + f(i,j,2*k-1)
          enddo
        enddo
      enddo

      ![NOTE: should cache the plan for the next time]
      do d = 1,3
        call RFFTW_F77_DESTROY_PLAN( plan(d) )
      enddo

#elif defined (FFTW_3)

! Things are somewhat different with FFTW Version 3.
! We can use the 1-d sine transform, rather than worrying about complex
! transforms.

      do d = 1 ,3
#ifdef CH_USE_DOUBLE
        call DFFTW_PLAN_R2R_1D
#else
        call SFFTW_PLAN_R2R_1D
#endif
     &                        ( plan(d) ,n(d)-1, in1d, out1d
     &                         ,FFTW_RODFT00, FFTW_ESTIMATE )
      enddo

![NOTE: something is wrong with the initializations of in1d()
!       inside the loops, since the first direction initializes
!       the N-1 element, but the other directions initialize the
!       zero element.  Dan Martin doesn't remember why he did this.
!       I'm concerned it's a bug waiting to happen.  <dbs>]

!     Sine transform in the first direction.
      d = 1
      do k=0,N3
        do j = 0,N2
          in1d(N1-1) = 0.d0
          out1d(N1-1) = 0.d0
          do i = 1,N1-1
            in1d(i-1) = f(i,j,k)
            out1d(i-1) = 0.0d0
          enddo
#ifdef CH_USE_DOUBLE
          call DFFTW_EXECUTE( plan(d) )
#else
          call SFFTW_EXECUTE( plan(d) )
#endif
          do i = 0,N1-1
            f(i+1,j,k) = out1d(i)
          enddo
        enddo
      enddo

!     Sine transform in the second direction.
      d = 2
      do k=0,N3
        do i = 0,N1
![NOTE: this is redundant, since the first iteration of the j loop sets in1d(0)]
!            in1d(0) = 0.d0
![NOTE: the non-cubic case uses different numbers of elements in in1d and out1d
!       in different directions so re-initialize each time.]
          in1d(N2-1) = 0
          out1d(N2-1) = 0
          do j = 1,N2-1
            in1d(j-1) = f(i,j,k)
          enddo
#ifdef CH_USE_DOUBLE
          call DFFTW_EXECUTE( plan(d) )
#else
          call SFFTW_EXECUTE( plan(d) )
#endif
          do j = 0,N2-1
            f(i,j+1,k) = out1d(j)
          enddo
        enddo
      enddo

!     Sine transform in the z direction.
      d = 3
      do i = 0,N1
        do j=0,N2
![NOTE: see comment above.]
!            in1d(0) = 0.d0
          in1d(N3-1) = 0
          out1d(N3-1) = 0
          do k = 1,N3-1
            in1d(k-1) = f(i,j,k)
          enddo
#ifdef CH_USE_DOUBLE
          call DFFTW_EXECUTE( plan(d) )
#else
          call SFFTW_EXECUTE( plan(d) )
#endif
          do k = 0,N3-1
            f(i,j,k+1) = out1d(k)
          enddo
        enddo
      enddo

!     Finally, rescale so that scaling matches original transform.
![NOTE: I assume by "original" he means the FFTW Version2 implementation above. <dbs>]
      do k=0, N3
        do j=0,N2
          do i=0,N1
            f(i,j,k) = f(i,j,k)/8.0d0
          enddo
        enddo
      enddo

      do d=1,3
#ifdef CH_USE_DOUBLE
        call DFFTW_DESTROY_PLAN( plan(d) )
#else
        call SFFTW_DESTROY_PLAN( plan(d) )
#endif
      enddo
#else
      ! NEED FFTW                 !!
      call MAYDAY_ERROR()

#endif

      return
      end
