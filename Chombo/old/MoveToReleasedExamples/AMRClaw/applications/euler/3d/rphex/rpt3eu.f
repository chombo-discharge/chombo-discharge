      subroutine rpt3  (ixyz,icoor,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,maux,ilr,asdq,
     &                  bmasdq,bpasdq)
c     ==================================================================
c
c     # Riemann solver in the transverse direction for the
c     # Euler equations.
c     #
c     # Uses Roe averages and other quantities which were
c     # computed in rpn3eu and stored in the common block comroe.
c     #
c     #
c     # On input,
c
c     #    ql,qr is the data along some one-dimensional slice, as in rpn3
c     #         This slice is
c     #             in the x-direction if ixyz=1,
c     #             in the y-direction if ixyz=2, or
c     #             in the z-direction if ixyz=3.
c     #    asdq is an array of flux differences (A^*\Dq).
c     #         asdq(i,:) is the flux difference propagating away from
c     #         the interface between cells i-1 and i.
c     #    Note that asdq represents B^*\Dq if ixyz=2 or C^*\Dq if ixyz=3.
c
c     #    ixyz indicates the direction of the original Riemann solve,
c     #         called the x-like direction in the table below:
c
c     #               x-like direction   y-like direction   z-like direction
c     #      ixyz=1:        x                  y                  z
c     #      ixyz=2:        y                  z                  x
c     #      ixyz=3:        z                  x                  y
c
c     #    icoor indicates direction in which the transverse solve should
c     #         be performed.
c     #      icoor=2: split in the y-like direction.
c     #      icoor=3: split in the z-like direction.
c
c     #    For example,
c     #      ixyz=1, icoor=2 means asdq=A^*\Dq, and should be split in y into
c     #                        bmasdq = B^-A^*\Dq,
c     #                        bpasdq = B^+A^*\Dq.
c     #
c     #      ixyz=2, icoor=2 means asdq=B^*\Dq, and should be split in z into
c     #                        bmasdq = C^-B^*\Dq,
c     #                        bpasdq = C^+B^*\Dq.
c
c     #    The parameter imp is generally needed only if aux
c     #    arrays are being used, in order to access the appropriate
c     #    variable coefficients.
c
c
c
      implicit none
      integer ixyz,icoor, maxm, meqn, mwaves,mbc,mx,maux,ilr

      double precision      ql(1-mbc:maxm+mbc, meqn)
      double precision      qr(1-mbc:maxm+mbc, meqn)
      double precision    asdq(1-mbc:maxm+mbc, meqn)
      double precision  bmasdq(1-mbc:maxm+mbc, meqn)
      double precision  bpasdq(1-mbc:maxm+mbc, meqn)
      double precision    aux1(1-mbc:maxm+mbc, maux, 3)
      double precision    aux2(1-mbc:maxm+mbc, maux, 3)
      double precision    aux3(1-mbc:maxm+mbc, maux, 3)
c

      double precision wave_rot(5,3),s_rot(3), asdq_proj(5), uvw(3)
      double precision uvwc(3), rot(9), gamma, gamma1

      double precision pres, enth, area
      integer mcapa, locrot, locarea, i,j, i1, mws,m

      common /param/ gamma, gamma1

      call get_aux_locations(ixyz,icoor,mcapa,locrot,locarea)

c     # Solve Riemann problem in the second coordinate direction

      do i = 2-mbc, mx+mbc
         i1 = i + ilr - 2

c        # compute values needed for Jacobian
         uvwc(1) = ql(i1,2)/ql(i1,1)
         uvwc(2) = ql(i1,3)/ql(i1,1)
         uvwc(3) = ql(i1,4)/ql(i1,1)
         pres = gamma1*(ql(i1,5)  - 0.5d0*(ql(i1,2)**2 +
     &         ql(i1,3)**2 + ql(i1,4)**2)/ql(i1,1))
         enth = (ql(i1,5) + pres)/ql(i1,1)


c        # -------------------------------------------------------
c        # Compute bmasdq
c        # -------------------------------------------------------

c        # Using aux2 works for both y-like and z-like directions
         do j = 1,6
            rot(j) = aux2(i1,locrot+j-1,2)
         enddo
         call compute_binormal(rot)

         area = aux2(i1,locarea,2)

         do m = 1,meqn
            asdq_proj(m) = asdq(i,m)
         enddo
         call rotate3(rot,asdq_proj(2))

         do j = 1,3
            uvw(j) = uvwc(j)
         enddo
         call rotate3(rot,uvw)

         call solve_riemann(uvw,enth,asdq_proj, wave_rot,s_rot)

         do mws = 1,mwaves
            call rotate3_tr(rot,wave_rot(2,mws))
         enddo

         do m=1,meqn
            bmasdq(i,m) = 0.d0
            do mws=1,mwaves
               bmasdq(i,m) = bmasdq(i,m)
     &               + area*min(s_rot(mws), 0.d0) * wave_rot(m,mws)
            enddo
         enddo

c        # -------------------------------------------------------
c        # Compute bpasdq
c        # -------------------------------------------------------

c        # Using aux2 works for both y-like and z-like directions
         if (icoor .eq. 2) then
            do j = 1,6
               rot(j) = aux3(i1,locrot+j-1,2)
            enddo
            area = aux3(i1,locarea,2)
         else
            do j = 1,6
               rot(j) = aux2(i1,locrot+j-1,3)
            enddo
            area = aux2(i1,locarea,3)
         endif
         call compute_binormal(rot)


         do m = 1,meqn
            asdq_proj(m) = asdq(i,m)
         enddo

         call rotate3(rot,asdq_proj(2))

         do j = 1,3
            uvw(j) = uvwc(j)
         enddo
         call rotate3(rot,uvw)

         call solve_riemann(uvw,enth,asdq_proj, wave_rot,s_rot)

         do mws = 1,mwaves
            call rotate3_tr(rot,wave_rot(2,mws))
         enddo

         do m=1,meqn
            bpasdq(i,m) = 0.d0
            do mws=1,mwaves
               bpasdq(i,m) = bpasdq(i,m)
     &               + area*max(s_rot(mws),0.d0) * wave_rot(m,mws)
            enddo
         enddo

      enddo  !! end of i loop

      return
      end
