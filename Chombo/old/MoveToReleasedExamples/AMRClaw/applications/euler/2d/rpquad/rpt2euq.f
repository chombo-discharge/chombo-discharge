c     =====================================================
      subroutine rpt2(ixy,maxm,meqn,mwaves,mbc,mx,
     &                  ql,qr,aux1,aux2,aux3,
     &                  ilr,asdq,bmasdq,bpasdq)
c     =====================================================
      implicit none
c
c     # Riemann solver in the transverse direction for the Euler equations.
c     # Split asdq (= A^* \Delta q, where * = + or -)
c     # into down-going flux difference bmasdq (= B^- A^* \Delta q)
c     #    and up-going flux difference bpasdq (= B^+ A^* \Delta q)
c
c     # Uses Roe averages and other quantities which were
c     # computed in rpn2eu and stored in the common block comroe.
c
      integer ixy, maxm, meqn, mwaves, mbc, mx, ilr

      double precision      ql(1-mbc:maxm+mbc, meqn)
      double precision      qr(1-mbc:maxm+mbc, meqn)
      double precision    asdq(1-mbc:maxm+mbc, meqn)
      double precision  bmasdq(1-mbc:maxm+mbc, meqn)
      double precision  bpasdq(1-mbc:maxm+mbc, meqn)
      double precision  aux1(1-mbc:maxm+mbc,*)
      double precision  aux2(1-mbc:maxm+mbc,*)
      double precision  aux3(1-mbc:maxm+mbc,*)

      double precision gamma, gamma1
      double precision dtcom, dxcom, dycom, tcom, icom, jcom
      integer i, i1, m, mw, locarea, mcapa, locrot
      double precision pres, enth, area

      double precision s_rot(3), wave_rot(4,3), asdq_rot(4), uv(2)
      double precision rot(4), uvc(2)

      common /param/  gamma,gamma1
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom


c     # Get aux locations for mapped grid stuff
      call get_aux_locations(ixy,mcapa,locrot,locarea)

      do i = 2-mbc, mx+mbc
         i1 = i + ilr - 2


c        # Just compute Jacobian from cell centered values.
         uvc(1) = ql(i1,2)/ql(i1,1)
         uvc(2) = ql(i1,3)/ql(i1,1)
         pres = gamma1*(ql(i1,4)  - 0.5d0*(ql(i1,2)**2 +
     &         ql(i1,3)**2)/ql(i1,1))
         enth = (ql(i1,4)+pres) / ql(i1,1)

c        # ------------------------------------------------
c        # Do bottom (minus) edge
c        # ------------------------------------------------

         rot(1) = aux2(i1,locrot)
         rot(2) = aux2(i1,locrot+1)
         call compute_tangent(rot)

         do m = 1,meqn
            asdq_rot(m) = asdq(i,m)
         enddo
         call rotate2(rot,asdq_rot(2))

         uv(1) = uvc(1)
         uv(2) = uvc(2)

         call rotate2(rot,uv)
         call solve_riemann(uv,enth,asdq_rot,wave_rot,s_rot)

c        # rot-T = [rot(1,1)  rot(2,1); rot(1,2) rot(2,2)]
         do mw = 1,mwaves
            call rotate2_tr(rot,wave_rot(2,mw))
         enddo

         area = aux2(i1,locarea)
         do m=1,meqn
            bmasdq(i,m) = 0.d0
            do mw = 1,mwaves
               bmasdq(i,m) = bmasdq(i,m) +
     &               area*min(s_rot(mw),0.d0)*wave_rot(m,mw)
            enddo
         enddo


c        # ----------------------------------------------------
c        # Do upper (plus) edge
c        # ----------------------------------------------------

         rot(1) = aux3(i1,locrot)
         rot(2) = aux3(i1,locrot+1)
         call compute_tangent(rot)

         do m = 1,meqn
            asdq_rot(m) = asdq(i,m)
         enddo

         call rotate2(rot,asdq_rot(2))

         uv(1) = uvc(1)
         uv(2) = uvc(2)
         call rotate2(rot,uv)

         call solve_riemann(uv,enth,asdq_rot,wave_rot,s_rot)

         do mw = 1,mwaves
            call rotate2_tr(rot,wave_rot(2,mw))
         enddo

         area = aux3(i1,locarea)
         do m=1,meqn
            bpasdq(i,m) = 0.d0
            do mw = 1,mwaves
               bpasdq(i,m) = bpasdq(i,m) +
     &               area*max(s_rot(mw),0.d0)*wave_rot(m,mw)
            enddo
         enddo
      enddo

      end
