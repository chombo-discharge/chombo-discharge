C      _______              __
C     / ___/ /  ___  __ _  / /  ___
C    / /__/ _ \/ _ \/  V \/ _ \/ _ \
C    \___/_//_/\___/_/_/_/_.__/\___/
C    Please refer to Copyright.txt, in Chombo's root directory.

      integer function lsworksize(neqns, nvars)
!     Returns size of work array required in solvels.
!
!     We solve a least-squares problem derived from equating the value
!     at each point in the neighborhood to a linear combination of
!     moments of displacement from the base point of the neighborhood;
!     we solve for the coefficients in the linear combination, which are
!     the same for every point.
!     There should be more points (equations) than coefficients:
!     that is, the system of equations is overdetermined.
!
!     => inputs,  <= outputs
!     neqns => number of equations; is number of points in each neighborhood
!     nvars => number of variables; this is number of coefficients
!
!     petermc, 21 Dec 2007
      integer neqns, nvars
      integer nrhs, lwork, ier
      parameter (nrhs=1)
#ifdef CH_USE_DOUBLE
      double precision mat(neqns, nvars), vec(neqns, nrhs), work(1)
#else
      real mat(neqns, nvars), vec(neqns, nrhs), work(1)
#endif

      lwork = -1
#ifdef CH_USE_DOUBLE
      call DGELS('N', neqns, nvars, nrhs, mat, neqns, vec, neqns,
     &     work, lwork, ier)
#else
      call SGELS('N', neqns, nvars, nrhs, mat, neqns, vec, neqns,
     &     work, lwork, ier)
#endif

      if (ier .ne. 0) then
#ifdef CH_USE_DOUBLE
         print *, 'Error in lsworksize:  DGELS returned ', ier
#else
         print *, 'Error in lsworksize:  SGELS returned ', ier
#endif
         call MAYDAY_ERROR()
      endif

      lsworksize = work(1)

      return
      end



      subroutine solvels(neqns, nvars, mat, vec, lwork, work)
!     Calls LAPACK routine to solve least-squares problem.
!
!     We solve a least-squares problem derived from equating the value
!     at each point in the neighborhood to a linear combination of
!     moments of displacement from the base point of the neighborhood;
!     we solve for the coefficients in the linear combination, which are
!     the same for every point.
!     There should be more points (equations) than coefficients:
!     that is, the system of equations is overdetermined.
!
!     => inputs,  <= outputs
!     neqns => number of equations; is number of points in each neighborhood
!     nvars => number of variables; this is number of coefficients
!     mat => matrix of coefficients
!     vec <=> input right-hand side, output solution
!     lwork => length of work vector, lsworksize(neqns, nvars)
!     work <=> work vector, first element returned as condition number
!
!     petermc, 21 Dec 2007
      integer nvars, neqns, lwork, ier, nrhs
!      integer ieqn, ivar
      integer iwork(nvars)
      parameter (nrhs=1)
#ifdef CH_USE_DOUBLE
      double precision mat(neqns, nvars), vec(neqns, nrhs), work(lwork)
      double precision rcond, cwork(3*nvars)
#else
      real mat(neqns, nvars), vec(neqns, nrhs), work(lwork)
      real rcond, cwork(3*nvars)
#endif

!     rhs is vec(1:neqns, 1)
!     matrix is mat(1:neqns, 1:nvars)
#ifdef CH_USE_DOUBLE
      call DGELS('N', neqns, nvars, nrhs, mat, neqns, vec, neqns,
     &     work, lwork, ier)
#else
      call SGELS('N', neqns, nvars, nrhs, mat, neqns, vec, neqns,
     &     work, lwork, ier)
#endif
!     solution is now in vec(1:nvars, 1)

!     Recall neqns > nvars in mat(neqns, nvars).
!     On return from {D|S}GELS, mat is overwritten by details of its
!     QR factorization as returned by {D|S}GEQRF.
!     That is:  the elements on and above the diagonal of mat
!     contain the nvars-by-nvars upper triangular matrix R;
!     the elements below the diagonal, with the array TAU, represent
!     the orthogonal matrix Q as a product of nvars elementary reflectors.

      if (ier .ne. 0) then
#ifdef CH_USE_DOUBLE
         print *, 'Error:  DGELS returned ', ier
#else
         print *, 'Error:  SGELS returned ', ier
#endif
         if (ier .lt. 0) then
            print *, 'illegal value in argument number ', -ier
         else
            print *, 'matrix does not have full rank'
         endif
         call MAYDAY_ERROR()
      endif

!     Get condition number.
!     First arg is 'I' for infinity norm or '1' for max norm.
!     Second arg is 'U' for upper triangular or 'L' for lower triangular.
!     Third arg is 'N' to use diagonal entries or 'U' to replace with 1.
#ifdef CH_USE_DOUBLE
      call DTRCON('I', 'U', 'N', nvars, mat, neqns, rcond,
     &     cwork, iwork, ier)
#else
      call STRCON('I', 'U', 'N', nvars, mat, neqns, rcond,
     &     cwork, iwork, ier)
#endif

      work(1) = 1.0 / rcond

      if (ier .lt. 0) then
#ifdef CH_USE_DOUBLE
         print *, 'Error:  DTRCON returned ', ier
#else
         print *, 'Error:  STRCON returned ', ier
#endif
         print *, 'illegal value in argument number ', -ier
         call MAYDAY_ERROR()
      endif

      return
      end
