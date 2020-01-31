
c
c
c     =====================================================
      subroutine rpn2(ixy,maxm,meqn,mwaves,mbc,mx,ql,qr,auxl,auxr,
     &                  wave,s,amdq,apdq)
c     =====================================================
c
c     # Roe-solver for the Euler equations
c     # solve Riemann problems along one slice of data.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # This data is along a slice in the x-direction if ixy=1
c     #                            or the y-direction if ixy=2.
c     # On output, wave contains the waves, s the speeds,
c     # and amdq, apdq the decomposition of the flux difference
c     #   f(qr(i-1)) - f(ql(i))
c     # into leftgoing and rightgoing parts respectively.
c     # With the Roe solver we have
c     #    amdq  =  A^-s \Delta q    and    apdq  =  A^+ \Delta q
c     # where A is the Roe matrix.  An entropy fix can also be incorporated
c     # into the flux differences.
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit none
c
      integer ixy, maxm, meqn, mwaves, mbc, mx

      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision  apdq(1-mbc:maxm+mbc, meqn)
      double precision  amdq(1-mbc:maxm+mbc, meqn)
      double precision auxl(1-mbc:maxm+mbc,*)
      double precision auxr(1-mbc:maxm+mbc,*)
c
c     local arrays -- common block comroe is passed to rpt2eu
c     ------------
      double precision   ql_proj(4), qr_proj(4), uv(2), wave_rot(4,3),
     &      apdq_rot(4), amdq_rot(4), delta(4), rot(4), s_rot(3)

      logical efix

      double precision gamma, gamma1
      double precision dtcom, dxcom, dycom, tcom, icom, jcom
      integer mcapa, locrot, locarea, i, m, mws
      double precision rhsqrtl, rhsqrtr, rhsq2, pl, pr, enth
      double precision rhoim1, pim1, cim1, s0
      double precision rho1, rhou1, rhov1, en1, p1, c1, s1, sfract
      double precision rhoi, pi, ci, s3, rho2, rhou2, rhov2, en2
      double precision p2, c2, s2, df, area, u2v2


      common /param/  gamma,gamma1
      common /comxyt/ dtcom,dxcom,dycom,tcom,icom,jcom

      data efix /.true./    !# use entropy fix for transonic rarefactions
c
c     # compute the Roe-averaged variables needed in the Roe solver.
c     # These are stored in the common block comroe since they are
c     # later used in routine rpt2eu to do the transverse wave splitting.
c

c     # Get aux locations for mapped grid stuff
      call get_aux_locations_n(ixy,mcapa,locrot,locarea)

c     # Compute waves at cell interfaces, by solving Riemann problem
c     # between cells i-1 and i.
      do i = 2-mbc,mx+mbc

         rot(1) = auxl(i,locrot)
         rot(2) = auxl(i,locrot+1)
         call compute_tangent(rot)

c        Project left and right states into normal direction, (ax,ay), and
c        # tangential directions (-ay, ax).
         do m = 1,meqn
            ql_proj(m) = ql(i,m)
            qr_proj(m) = qr(i-1,m)
         enddo

         call rotate2(rot,ql_proj(2))
         call rotate2(rot,qr_proj(2))

c        # Compute Roe averages using rotated data
         rhsqrtl = sqrt(qr_proj(1))
         rhsqrtr = sqrt(ql_proj(1))
         rhsq2 = rhsqrtl + rhsqrtr

         uv(1) = (qr_proj(2)/rhsqrtl + ql_proj(2)/rhsqrtr)/rhsq2
         uv(2) = (qr_proj(3)/rhsqrtl + ql_proj(3)/rhsqrtr)/rhsq2

         pl = gamma1*(qr_proj(4) - 0.5d0*(qr_proj(2)**2 +
     &         qr_proj(3)**2)/qr_proj(1))
         pr = gamma1*(ql_proj(4) - 0.5d0*(ql_proj(2)**2 +
     &         ql_proj(3)**2)/ql_proj(1))

         enth = (((qr_proj(4)+pl)/rhsqrtl
     &         + (ql_proj(4)+pr)/rhsqrtr)) / rhsq2

         do m = 1,meqn
            delta(m) = ql_proj(m) - qr_proj(m)
         enddo

         call solve_riemann(uv,enth,delta,wave_rot,s_rot)

         if (efix) then

c           # check 1-wave:
c           ---------------
            u2v2 = qr(i-1,2)**2 + qr(i-1,3)**2
            rhoim1 = qr_proj(1)
            pim1 = gamma1*(qr_proj(4) - 0.5d0*u2v2*rhoim1)
            cim1 = sqrt(gamma*pim1/rhoim1)
            s0 = qr_proj(2)/rhoim1 - cim1 !# u-c in left state (cell i-1)

c           # check for fully supersonic case:
            if (s0 .ge. 0.d0 .and. s_rot(1).gt.0.d0)  then
c              # everything is right-going
               do m=1,meqn
                  amdq_rot(m) = 0.d0
               enddo
               go to 200
            else
               rho1  = qr_proj(1)  + wave_rot(1,1)
               rhou1 = qr_proj(2)  + wave_rot(2,1)
               rhov1 = qr_proj(3)  + wave_rot(3,1)
               en1   = qr_proj(4)  + wave_rot(4,1)
               p1 = gamma1*(en1 - 0.5d0*(rhou1**2 + rhov1**2)/rho1)
               c1 = sqrt(gamma*p1/rho1)
               s1 = rhou1/rho1 - c1     !# u-c to right of 1-wave
               if (s0 .lt. 0.d0 .and. s1 .gt. 0.d0) then
c                 # transonic rarefaction in the 1-wave
                  sfract = s0 * (s1 - s_rot(1)) / (s1-s0)
               else if (s_rot(1) .lt. 0.d0) then
c                 # 1-wave is leftgoing
                  sfract = s_rot(1)
               else
c                 # 1-wave is rightgoing
                  sfract = 0.d0         !# this shouldn't happen since s0 < 0
               endif
               do m = 1,4
                  amdq_rot(m) = sfract*wave_rot(m,1)
               enddo
            endif
c
c           # check 2-wave:
c           ---------------
c
            if (s_rot(2) .ge. 0.d0) go to 200

c           # This gets added to whatever was computed above.
            do m=1,meqn
               amdq_rot(m) = amdq_rot(m) + s_rot(2)*wave_rot(m,2)
            enddo

c
c           # check 3-wave:
c           ---------------
c

            u2v2 = ql(i,2)**2 + ql(i,3)**2
            rhoi = ql_proj(1)
            pi = gamma1*(ql_proj(4) - 0.5d0*u2v2*rhoi)
            ci = dsqrt(gamma*pi/rhoi)
            s3 = ql_proj(2)/rhoi + ci   !# u+c in right state  (cell i)
c
            rho2  = ql_proj(1) - wave_rot(1,3)
            rhou2 = ql_proj(2) - wave_rot(2,3)
            rhov2 = ql_proj(3) - wave_rot(3,3)
            en2   = ql_proj(4) - wave_rot(4,3)
            p2 = gamma1*(en2 - 0.5d0*(rhou2**2 + rhov2**2)/rho2)
            c2 = dsqrt(gamma*p2/rho2)
            s2 = rhou2/rho2 + c2        !# u+c to left of 3-wave
            if (s2 .lt. 0.d0 .and. s3.gt.0.d0) then
c              # transonic rarefaction in the 3-wave
               sfract = s2 * (s3-s_rot(3)) / (s3-s2)
            else if (s_rot(3) .lt. 0.d0) then
c              # 3-wave is leftgoing
               sfract = s_rot(3)
            else
c              # 3-wave is rightgoing
               go to 200
            endif
c
            do m=1,meqn
               amdq_rot(m) = amdq_rot(m) + sfract*wave_rot(m,3)
            enddo

  200       continue
c           # Compute apdq from amdq : ApDq + AmDq = sum s*waves
            do  m=1,meqn
               df = 0.d0
               do mws=1,mwaves
                  df = df + s_rot(mws)*wave_rot(m,mws)
               enddo
               apdq_rot(m) = df - amdq_rot(m)
            enddo
         else
c           # no efix
            do m=1,meqn
               amdq_rot(m) = 0.d0
               apdq_rot(m) = 0.d0
               do  mws=1,mwaves
                  if (s_rot(mws) .lt. 0.d0) then
                     amdq_rot(m) = amdq_rot(m) +
     &                     s_rot(mws)*wave_rot(m,mws)
                  else
                     apdq_rot(m) = apdq_rot(m) +
     &                     s_rot(mws)*wave_rot(m,mws)
                  endif
               enddo
            enddo

         endif  !! end of efix if

c        # Rotate and scale everything
         area = auxl(i,locarea)
         do mws = 1,mwaves
c           # Now scale the speeds
            s(i,mws) = area*s_rot(mws)

c           # rotate waves
            call rotate2_tr(rot,wave_rot(2,mws))
            do m = 1,meqn
               wave(i,m,mws) = wave_rot(m,mws)
            enddo
         enddo

c        # Rotate amdq, apdq
         call rotate2_tr(rot,apdq_rot(2))
         call rotate2_tr(rot,amdq_rot(2))
         do m = 1,meqn
            apdq(i,m) = area*apdq_rot(m)
            amdq(i,m) = area*amdq_rot(m)
         enddo

      enddo

      return
      end


      subroutine solve_riemann(uv,enth,delta,wave_rot,s_rot)
      implicit none

      double precision uv(2), enth, delta(4), wave_rot(4,3),
     &      s_rot(3)

      double precision u2v2, a2, a, g1a2, euv
      double precision gamma, gamma1, a1, a3, a4
      double precision u,v

      common /param/ gamma, gamma1

      u = uv(1)
      v = uv(2)

      u2v2 = u**2 + v**2
      a2 = gamma1*(enth - .5d0*u2v2)
      if (a2 .lt. 0) then
         write(6,*) 'a2 < 0; a2 = ', a2
         stop
      endif
      a = sqrt(a2)
      g1a2 = gamma1 / a2
      euv = enth - u2v2
c
      a3 = g1a2 *
     &      (euv*delta(1) + u*delta(2) + v*delta(3) - delta(4))
      a2 = delta(3) - v*delta(1)
      a4 = (delta(2) + (a-u)*delta(1) - a*a3)
     &      / (2.d0*a)
      a1 = delta(1) - a3 - a4
c
c     # Compute the waves.
c     # Note that the 2-wave and 3-wave travel at the same speed and
c     # are lumped together in wave(.,.,2).  The 4-wave is then stored in
c     # wave(.,.,3).
c
      wave_rot(1,1)  = a1
      wave_rot(2,1)  = a1*(u-a)
      wave_rot(3,1)  = a1*v
      wave_rot(4,1)  = a1*(enth - u*a)
      s_rot(1)       = u - a

      wave_rot(1,2)  = a3
      wave_rot(2,2)  = a3*u
      wave_rot(3,2)  = a3*v          + a2
      wave_rot(4,2)  = a3*0.5d0*u2v2 + a2*v
      s_rot(2)       = u

      wave_rot(1,3)  = a4
      wave_rot(2,3)  = a4*(u + a)
      wave_rot(3,3)  = a4*v
      wave_rot(4,3)  = a4*(enth + u*a)
      s_rot(3)       = u + a

      end
