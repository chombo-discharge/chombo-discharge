
#include "CONSTANTS.H"
#include "REAL.H"

c ----------------------------------------------------------------------
c  Transform an index
c
c  idxout   <=  Transformed index
c  idxin     => Original index
c  xctm      => Partially expanded compact transformation matrix
c               (even -> new index, odd -> sign)
c ----------------------------------------------------------------------

      subroutine ctTransformStoP(idxout, idxin, xctm)
      implicit none
      integer idxout(0:CH_SPACEDIM-1)
      integer idxin(0:CH_SPACEDIM-1)
      integer xctm(0:2*CH_SPACEDIM-1)
      integer i
      do i = 0, CH_SPACEDIM-1
         idxout(xctm(2*i)) = xctm(2*i+1)*idxin(i)
      enddo
      return
      end

c ----------------------------------------------------------------------
c  Get the transformed coordinate direction
c           <=  Transformed coordinate direction
c  i         => Original coordinate direction (0:SpaceDim-1)
c  xctm      => Partially expanded compact transformation matrix
c               (even -> new index, odd -> sign)
c ----------------------------------------------------------------------

      function ctIndexAt(i, xctm)
      implicit none
      integer ctIndexAt
      integer i
      integer xctm(0:2*CH_SPACEDIM-1)
      ctIndexAt = xctm(2*i)
      return
      end

c ----------------------------------------------------------------------
c  Get the sign change for a coordinate direction
c           <=  Possible sign change due to transformation (-1 or 1)
c  i         => Original coordinate direction (0:SpaceDim-1)
c  xctm      => Partially expanded compact transformation matrix
c               (even -> new index, odd -> sign)
c ----------------------------------------------------------------------

      function ctSignAt(i, xctm)
      implicit none
      integer ctSignAt
      integer i
      integer xctm(0:2*CH_SPACEDIM-1)
      ctSignAt = xctm(2*i+1)
      return
      end
      
c ----------------------------------------------------------------------
c  Convert a symmetric matrix stored only on LD to full LDU
c
c  N         => Order of M
c  A        <=  Symmetric matrix stored on LDU
c  A         => Symmetric matrix stored on LD
c  LDA       => Leading dimension of A
c ----------------------------------------------------------------------

      subroutine sy2ge(N, A, LDA)
      implicit none
      real_t A(LDA,*)
      integer N
      integer LDA
      integer i, j
      do j = 2, n
         do i = 1, j-1
            A(i, j) = A(j, i)
         enddo
      enddo
      return
      end

c ----------------------------------------------------------------------
c  Comparison with limit tol^2 as x and y -> 0.
c           <=  T - not equal
c  x         => First float
c  y         => Second float
c  prec      => Number of significant digits to compare
c ----------------------------------------------------------------------

      function FPCompare(x, y, prec)
      implicit none
      logical FPCompare
      integer prec
      real_t x, y, tol
      tol = 10.D0**(-ABS(prec))
      FPCompare = ((ABS(x - y)) .gt. ((MIN(ABS(x), ABS(y)) + tol)*tol))
      return
      end

c ----------------------------------------------------------------------
c  Comparison with limit tol as x and y -> 0.
c           <=  T - not equal
c  x         => First float
c  y         => Second float
c  prec      => Number of significant digits to compare
c ----------------------------------------------------------------------

      function FPCompare1(x, y, prec)
      implicit none
      logical FPCompare1
      integer prec
      real_t x, y, tol
      tol = 10.D0**(-ABS(prec))
      FPCompare1 = ((ABS(x - y)) .gt. (MIN(ABS(x), ABS(y))*tol + tol))
      return
      end

c ----------------------------------------------------------------------
c  Compute a binomial
c           <=  Binomial coefficient
c  n         => Numerator
c  k         => Denominator
c ----------------------------------------------------------------------

      function binomial(n, k)
      implicit none
      integer binomial
      integer n
      integer k
      integer cnum, cden, j
      if (k .eq. 0) then
         cnum = 1
         cden = 1
      else
         cnum = n - k +1
         cden = k
      endif
      do j = 1, k-1
         cnum = cnum*(n - j + 1)
         cden = cden*j
      enddo
      binomial = cnum/cden
      return
      end

c ----------------------------------------------------------------------
c  Find the sequential index of a power (see C++ routine in
c  'FourthOrderMappedFineInterpSup.H' for more information
c           <=  Sequential index
c  m         => Degree of polynomial
c  p         => Power at each spatial index
c ----------------------------------------------------------------------

      function powerIndex(m, p)
      implicit none
      integer powerIndex
      integer m
      integer p(0:CH_SPACEDIM-1)
      integer binomial
      integer idx, mi, idir, remdir
      idx = 0
      mi = m
      do idir = 0, CH_SPACEDIM-1
         remdir = CH_SPACEDIM - idir
         idx = idx + binomial(remdir + mi, remdir)
         mi = mi - p(idir)
         idx = idx - binomial(remdir + mi, remdir)
      enddo
      powerIndex = idx
      return
      end
