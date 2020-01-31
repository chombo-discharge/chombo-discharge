c
c     ==================================================================
      subroutine rpn3  (ixyz,maxm,meqn,mwaves,mbc,mx,ql,qr,
     &			auxl,auxr,maux,wave,s,amdq,apdq)
c     ==================================================================
c
c     # Roe-solver for the Euler equations
c      -----------------------------------------------------------
c
c     # solve Riemann problems along one slice of data.
c     # This data is along a slice in the x-direction if ixyz=1
c     #                               the y-direction if ixyz=2.
c     #                               the z-direction if ixyz=3.
c
c     # On input, ql contains the state vector at the left edge of each cell
c     #           qr contains the state vector at the right edge of each cell
c
c     # On output, wave contains the waves, s the speeds,
c     # and amdq, apdq the left-going and right-going flux differences,
c     # respectively.
c
c     # Note that the i'th Riemann problem has left state qr(i-1,:)
c     #                                    and right state ql(i,:)
c     # From the basic clawpack routines, this routine is called with ql = qr
c
c
      implicit double precision(a-h,o-z)
c
      dimension wave(1-mbc:maxm+mbc, meqn, mwaves)
      dimension    s(1-mbc:maxm+mbc, mwaves)
      dimension   ql(1-mbc:maxm+mbc, meqn)
      dimension   qr(1-mbc:maxm+mbc, meqn)
      dimension amdq(1-mbc:maxm+mbc, meqn)
      dimension apdq(1-mbc:maxm+mbc, meqn)
      dimension auxl(1-mbc:maxm+mbc, maux)
      dimension auxr(1-mbc:maxm+mbc, maux)

      dimension ql_proj(5),qr_proj(5), s_rot(3), wave_rot(5,3)
      dimension amdq_rot(5), apdq_rot(5), rot(9), uvw(3)

c
c     local arrays -- common block comroe is passed to rpt3eu
c
c     ------------
      dimension delta(5)
      logical efix, debug

      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom
      common /param/ gamma, gamma1
c
      data efix /.true./    !# use entropy fix for transonic rarefactions
c
      if (mwaves .ne. 3) then
         write(6,*) '*** Should have mwaves=3 for this Riemann solver'
         stop
      endif

      call get_aux_locations_n(ixyz,mcapa,locrot,locarea)

      do i = 2-mbc, mx+mbc

c        # Get  two stored vectors; third is computed when rotations are done
         do j = 1,6
            rot(j) = auxl(i,locrot+j-1)
         enddo
         call compute_binormal(rot)

         do m = 1,meqn
            ql_proj(m) = ql(i,m)
            qr_proj(m) = qr(i-1,m)
         enddo

         call rotate3(rot,ql_proj(2))
         call rotate3(rot,qr_proj(2))

c        # Use Roe averaged values for normal solves
         rhsqrtl = dsqrt(qr_proj(1))
         rhsqrtr = dsqrt(ql_proj(1))
         rhsq2 = rhsqrtl + rhsqrtr

         uvw(1) = (qr_proj(2)/rhsqrtl + ql_proj(2)/rhsqrtr) / rhsq2
         uvw(2) = (qr_proj(3)/rhsqrtl + ql_proj(3)/rhsqrtr) / rhsq2
         uvw(3) = (qr_proj(4)/rhsqrtl + ql_proj(4)/rhsqrtr) / rhsq2

         pl = gamma1*(qr_proj(5) - 0.5d0*(qr_proj(2)**2 +
     &         qr_proj(3)**2 + qr_proj(4)**2)/qr_proj(1))
         pr = gamma1*(ql_proj(5) - 0.5d0*(ql_proj(2)**2 +
     &         ql_proj(3)**2 + ql_proj(4)**2)/ql_proj(1))

         enth = (((qr_proj(5)+pl)/rhsqrtl
     &         + (ql_proj(5)+pr)/rhsqrtr)) / rhsq2

         do j = 1,5
            delta(j) = ql_proj(j) - qr_proj(j)
         enddo

c        # Solve normal Riemann problem
         call solve_riemann(uvw, enth, delta, wave_rot,s_rot)

         if (efix) then
c
c           # check 1-wave:
c           ---------------
c
            rhoim1 = qr_proj(1)
            pim1 = gamma1*(qr_proj(5) - 0.5d0*(qr_proj(2)**2
     &            + qr_proj(3)**2 + qr_proj(4)**2) / rhoim1)
            cim1 = dsqrt(gamma*pim1/rhoim1)
            s0 = qr_proj(2)/rhoim1 - cim1 !# u-c in left state (cell i-1)
c
c
c           # check for fully supersonic case:
            if (s0.ge.0.d0 .and. s_rot(1).gt.0.d0)then
c              # everything is right-going
               do 60 m=1,meqn
                  amdq_rot(m) = 0.d0
   60          continue
               go to 200
            else
c
               rho1 = qr_proj(1) + wave_rot(1,1)
               rhou1 = qr_proj(2) + wave_rot(2,1)
               rhov1 = qr_proj(3) + wave_rot(3,1)
               rhow1 = qr_proj(4) + wave_rot(4,1)
               en1 = qr_proj(5) + wave_rot(5,1)
               p1 = gamma1*(en1 - 0.5d0*(rhou1**2 + rhov1**2 +
     &               rhow1**2)/rho1)
               c1 = dsqrt(gamma*p1/rho1)
               s1 = rhou1/rho1 - c1     !# u-c to right of 1-wave_rot
               if (s0.lt.0.d0 .and. s1.gt.0.d0) then
c                 # transonic rarefaction in the 1-wave
                  sfract = s0 * (s1-s_rot(1)) / (s1-s0)
               else if (s_rot(1) .lt. 0.d0) then
c                 # 1-wave is leftgoing
                  sfract = s_rot(1)
               else
c                 # 1-wave is rightgoing
                  sfract = 0.d0         !# this shouldn't happen since s0 < 0
               endif
               do m=1,meqn
                  amdq_rot(m) = sfract*wave_rot(m,1)
               enddo
            endif
c
c           # check 2-wave:
c           ---------------
c
            if (s_rot(2) .ge. 0.d0) go to 200 !# 2-,3- and 4- wave_rots are rightgoing

            do m=1,meqn
               amdq_rot(m) = amdq_rot(m) + s_rot(2)*wave_rot(m,2)
            enddo
c
c           # check 3-wave:
c           ---------------
c
            rhoi = ql_proj(1)
            pi = gamma1*(ql_proj(5) - 0.5d0*(ql_proj(2)**2
     &            + ql_proj(3)**2 + ql_proj(4)**2) / rhoi)
            ci = dsqrt(gamma*pi/rhoi)
            s3 = ql_proj(2)/rhoi + ci     !# u+c in right state  (cell i)
c
            rho2 = ql_proj(1) - wave_rot(1,3)
            rhou2 = ql_proj(2) - wave_rot(2,3)
            rhov2 = ql_proj(3) - wave_rot(3,3)
            rhow2 = ql_proj(4) - wave_rot(4,3)
            en2 = ql_proj(5) - wave_rot(5,3)
            p2 = gamma1*(en2 - 0.5d0*(rhou2**2 + rhov2**2 +
     &            rhow2**2)/rho2)
            c2 = dsqrt(gamma*p2/rho2)
            s2 = rhou2/rho2 + c2        !# u+c to left of 3-wave
            if (s2 .lt. 0.d0 .and. s3.gt.0.d0 ) then
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
            do m=1,5
               amdq_rot(m) = amdq_rot(m) + sfract*wave_rot(m,3)
            enddo

  200       continue
c
c           # compute the rightgoing flux differences:
c           # df = SUM s*wave   is the total flux difference and apdq = df -
C           # amdq
c
            do m = 1,meqn
               df = 0.d0
               do mws=1,mwaves
                  df = df + s_rot(mws)*wave_rot(m,mws)
               enddo
               apdq_rot(m) = df - amdq_rot(m)
            enddo
         else
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

         endif  !! end of efix

         do mws = 1,mwaves
            call rotate3_tr(rot,wave_rot(2,mws))
            do m = 1,meqn
               wave(i,m,mws) = wave_rot(m,mws)
            enddo
         enddo

         area = auxl(i,locarea)
         do mws = 1,mwaves
            s(i,mws) = area*s_rot(mws)
         enddo

         call rotate3_tr(rot,amdq_rot(2))
         call rotate3_tr(rot,apdq_rot(2))
         do m = 1,meqn
            amdq(i,m) = area*amdq_rot(m)
            apdq(i,m) = area*apdq_rot(m)
         enddo

      enddo  !! end of i loop over 1d sweep array

      return
      end


      subroutine solve_riemann(uvw,enth,delta,wave_rot,s)
      implicit none

      double precision enth, uvw(3), delta(5)
      double precision wave_rot(5,3), s(3)
      double precision u2v2w2, a2, a, g1a2, euv
      double precision a1, a3, a4, a5,u,v,w
      double precision gamma, gamma1
      logical debug

      common /param/ gamma, gamma1


      u = uvw(1)
      v = uvw(2)
      w = uvw(3)

      u2v2w2 = u**2 + v**2 + w**2
      a2 = gamma1*(enth - 0.5d0*u2v2w2)
      if (a2 .lt. 0.d0) then
         write(6,*) 'a2 .lt. 0; ', a2
         stop
      endif
      a = dsqrt(a2)
      g1a2 = gamma1 / a2
      euv = enth - u2v2w2

      a4 = g1a2 * (euv*delta(1)
     &      + u*delta(2) + v*delta(3) + w*delta(4)
     &      - delta(5))
      a2 = delta(3) - v*delta(1)
      a3 = delta(4) - w*delta(1)
      a5 = (delta(2) + (a-u)*delta(1) - a*a4) / (2.d0*a)
      a1 = delta(1) - a4 - a5
c
c     # Compute the waves.
c     # Note that the 2-wave, 3-wave and 4-wave travel at the same speed
c     # and are lumped together in wave(.,.,2).  The 5-wave is then stored
c     # in wave(.,.,3).
c
      wave_rot(1,1)  = a1
      wave_rot(2,1) = a1*(u-a)
      wave_rot(3,1) = a1*v
      wave_rot(4,1) = a1*w
      wave_rot(5,1)  = a1*(enth - u*a)
      s(1) = u - a
c
      wave_rot(1,2)  = a4
      wave_rot(2,2) = a4*u
      wave_rot(3,2) = a4*v + a2
      wave_rot(4,2) = a4*w + a3
      wave_rot(5,2)  = a4*0.5d0*u2v2w2  + a2*v + a3*w
      s(2) = u
c
      wave_rot(1,3)  = a5
      wave_rot(2,3) = a5*(u+a)
      wave_rot(3,3) = a5*v
      wave_rot(4,3) = a5*w
      wave_rot(5,3)  = a5*(enth+u*a)
      s(3) = u + a

      end
