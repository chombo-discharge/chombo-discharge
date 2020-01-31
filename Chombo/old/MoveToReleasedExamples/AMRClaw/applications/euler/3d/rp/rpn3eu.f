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
      implicit none
      integer ixyz,maxm,meqn,mwaves,mbc,mx,maux
c
      double precision wave(1-mbc:maxm+mbc, meqn, mwaves)
      double precision    s(1-mbc:maxm+mbc, mwaves)
      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision amdq(1-mbc:maxm+mbc, meqn)
      double precision apdq(1-mbc:maxm+mbc, meqn)
      double precision auxl(1-mbc:maxm+mbc, maux)
      double precision auxr(1-mbc:maxm+mbc, maux)

      double precision waveb(5,3),sb(3)
      double precision a1,a2,a3,a4,a5
      integer mu,mv,mw, mws, i, m

      integer maxmrp
      parameter (maxmrp = 502)
      double precision u2v2w2(-1:maxmrp), u(-1:maxmrp), v(-1:maxmrp)
      double precision w(-1:maxmrp), enth(-1:maxmrp), a(-1:maxmrp)
      double precision g1a2(-1:maxmrp), euv(-1:maxmrp)
      double precision dtcom, dxcom, dycom, dzcom,tcom,gamma, gamma1
      integer icom, jcom, kcom

      double precision rhsqrtl, rhsqrtr, pl,pr, rhsq2
      double precision rhoim1, pim1, cim1,s0
      double precision rho1,rhou1,rhov1,rhow1,en1,p1,c1,s1,sfract
      double precision rhoi,pi,ci,s3,rho2,rhou2,rhov2,rhow2,en2,p2
      double precision c2,s2, df


      double precision delta(5)
      logical efix

      common /comxyzt/ dtcom,dxcom,dycom,dzcom,tcom,icom,jcom,kcom
      common /comroe/ u2v2w2,u,v,w,enth,a,g1a2,euv
      common /param/ gamma, gamma1
c
      data efix /.true./    !# use entropy fix for transonic rarefactions
c
      if (-1 .gt. 1-mbc .or. maxmrp .lt. maxm+mbc) then
	 write(6,*) 'need to increase maxmrp in rpn3eu'
         write(6,*) 'maxm+mbc=',maxm+mbc
	 stop
      endif

      if (mwaves .ne. 3) then
         write(6,*) '*** Should have mwaves=3 for this Riemann solver'
         stop
      endif

c
c     # set mu to point to  the component of the system that corresponds
c     # to momentum in the direction of this slice, mv and mw to the
c     # orthogonal momentums:
c
      if(ixyz .eq. 1)then
	  mu = 2
	  mv = 3
          mw = 4
      else if(ixyz .eq. 2)then
	  mu = 3
	  mv = 4
          mw = 2
      else
          mu = 4
          mv = 2
          mw = 3
      endif
c
c
c     # note that notation for u,v, and w reflects assumption that the
c     # Riemann problems are in the x-direction with u in the normal
c     # direction and v and w in the orthogonal directions, but with the
c     # above definitions of mu, mv, and mw the routine also works with
c     # ixyz=2 and ixyz = 3
c     # and returns, for example, f0 as the Godunov flux g0 for the
c     # Riemann problems u_t + g(u)_y = 0 in the y-direction.
c
c
c     # Compute the Roe-averaged variables needed in the Roe solver.
c     # These are stored in the common block comroe since they are
c     # later used in routine rpt3eu to do the transverse wave
c     # splitting.
c
      do 10 i = 2-mbc, mx+mbc
	 rhsqrtl = dsqrt(qr(i-1,1))
	 rhsqrtr = dsqrt(ql(i,1))
	 pl = gamma1*(qr(i-1,5) - 0.5d0*(qr(i-1,mu)**2 +
     &		 qr(i-1,mv)**2 + qr(i-1,mw)**2)/qr(i-1,1))
	 pr = gamma1*(ql(i,5) - 0.5d0*(ql(i,mu)**2 +
     &		 ql(i,mv)**2 + ql(i,mw)**2)/ql(i,1))
	 rhsq2 = rhsqrtl + rhsqrtr
	 u(i) = (qr(i-1,mu)/rhsqrtl + ql(i,mu)/rhsqrtr) / rhsq2
	 v(i) = (qr(i-1,mv)/rhsqrtl + ql(i,mv)/rhsqrtr) / rhsq2
	 w(i) = (qr(i-1,mw)/rhsqrtl + ql(i,mw)/rhsqrtr) / rhsq2
	 enth(i) = (((qr(i-1,5)+pl)/rhsqrtl
     &		   + (ql(i,5)+pr)/rhsqrtr)) / rhsq2
	 u2v2w2(i) = u(i)**2 + v(i)**2 + w(i)**2
         a2 = gamma1*(enth(i) - .5d0*u2v2w2(i))
         a(i) = dsqrt(a2)
	 g1a2(i) = gamma1 / a2
	 euv(i) = enth(i) - u2v2w2(i)
   10 continue
c
c
c     # now split the jump in q1d at each interface into waves
c
c     # find a1 thru a5, the coefficients of the 5 eigenvectors:
      do 20 i = 2-mbc, mx+mbc
         delta(1) = ql(i,1) - qr(i-1,1)
         delta(2) = ql(i,mu) - qr(i-1,mu)
         delta(3) = ql(i,mv) - qr(i-1,mv)
         delta(4) = ql(i,mw) - qr(i-1,mw)
         delta(5) = ql(i,5) - qr(i-1,5)
         a4 = g1a2(i) * (euv(i)*delta(1)
     &      + u(i)*delta(2) + v(i)*delta(3) + w(i)*delta(4)
     &      - delta(5))
         a2 = delta(3) - v(i)*delta(1)
         a3 = delta(4) - w(i)*delta(1)
         a5 = (delta(2) + (a(i)-u(i))*delta(1) - a(i)*a4) / (2.d0*a(i))
         a1 = delta(1) - a4 - a5
c
c        # Compute the waves.
c        # Note that the 2-wave, 3-wave and 4-wave travel at the same speed
c        # and are lumped together in wave(.,.,2).  The 5-wave is then stored
c        # in wave(.,.,3).
c
         wave(i,1,1)  = a1
         wave(i,mu,1) = a1*(u(i)-a(i))
         wave(i,mv,1) = a1*v(i)
         wave(i,mw,1) = a1*w(i)
         wave(i,5,1)  = a1*(enth(i) - u(i)*a(i))
         s(i,1) = u(i)-a(i)
c
         wave(i,1,2)  = a4
         wave(i,mu,2) = a4*u(i)
         wave(i,mv,2) = a4*v(i)	 	 + a2
         wave(i,mw,2) = a4*w(i)	 	 + a3
         wave(i,5,2)  = a4*0.5d0*u2v2w2(i)  + a2*v(i) + a3*w(i)
         s(i,2) = u(i)
c
         wave(i,1,3)  = a5
         wave(i,mu,3) = a5*(u(i)+a(i))
         wave(i,mv,3) = a5*v(i)
         wave(i,mw,3) = a5*w(i)
         wave(i,5,3)  = a5*(enth(i)+u(i)*a(i))
         s(i,3) = u(i)+a(i)
   20    continue
c
c
c    # compute flux differences amdq and apdq.
c    ---------------------------------------
c
      if (efix) go to 110
c
c     # no entropy fix
c     ----------------
c
c     # amdq = SUM s*wave   over left-going waves
c     # apdq = SUM s*wave   over right-going waves
c
      do 100 m=1,meqn
         do 100 i=2-mbc, mx+mbc
	    amdq(i,m) = 0.d0
	    apdq(i,m) = 0.d0
	    do 90 mws=1,mwaves
	       if (s(i,mws) .lt. 0.d0) then
		   amdq(i,m) = amdq(i,m) + s(i,mws)*wave(i,m,mws)
		 else
		   apdq(i,m) = apdq(i,m) + s(i,mws)*wave(i,m,mws)
		 endif
   90          continue
  100       continue
      go to 900
c
c-----------------------------------------------------
c
  110 continue
c
c     # With entropy fix
c     ------------------
c
c    # compute flux differences amdq and apdq.
c    # First compute amdq as sum of s*wave for left going waves.
c    # Incorporate entropy fix by adding a modified fraction of wave
c    # if s should change sign.
c
      do 200 i = 2-mbc, mx+mbc
c
c        # check 1-wave:
c        ---------------
c
	 rhoim1 = qr(i-1,1)
	 pim1 = gamma1*(qr(i-1,5) - 0.5d0*(qr(i-1,mu)**2
     &           + qr(i-1,mv)**2 + qr(i-1,mw)**2) / rhoim1)
	 cim1 = dsqrt(gamma*pim1/rhoim1)
	 s0 = qr(i-1,mu)/rhoim1 - cim1     !# u-c in left state (cell i-1)
c
c
c        # check for fully supersonic case:
	 if (s0.ge.0.d0 .and. s(i,1).gt.0.d0)then
c            # everything is right-going
	     do 60 m=1,meqn
		amdq(i,m) = 0.d0
   60           continue
	     go to 200
	     endif
c
         rho1 = qr(i-1,1) + wave(i,1,1)
         rhou1 = qr(i-1,mu) + wave(i,mu,1)
         rhov1 = qr(i-1,mv) + wave(i,mv,1)
         rhow1 = qr(i-1,mw) + wave(i,mw,1)
         en1 = qr(i-1,5) + wave(i,5,1)
         p1 = gamma1*(en1 - 0.5d0*(rhou1**2 + rhov1**2 +
     &                rhow1**2)/rho1)
         c1 = dsqrt(gamma*p1/rho1)
         s1 = rhou1/rho1 - c1  !# u-c to right of 1-wave
         if (s0.lt.0.d0 .and. s1.gt.0.d0) then
c            # transonic rarefaction in the 1-wave
	     sfract = s0 * (s1-s(i,1)) / (s1-s0)
	   else if (s(i,1) .lt. 0.d0) then
c	     # 1-wave is leftgoing
	     sfract = s(i,1)
	   else
c	     # 1-wave is rightgoing
             sfract = 0.d0   !# this shouldn't happen since s0 < 0
	   endif
	 do 120 m=1,meqn
	    amdq(i,m) = sfract*wave(i,m,1)
  120       continue
c
c        # check 2-wave:
c        ---------------
c
         if (s(i,2) .ge. 0.d0) go to 200  !# 2-,3- and 4- waves are rightgoing
	 do 140 m=1,meqn
	    amdq(i,m) = amdq(i,m) + s(i,2)*wave(i,m,2)
  140       continue
c
c        # check 3-wave:
c        ---------------
c
	 rhoi = ql(i,1)
	 pi = gamma1*(ql(i,5) - 0.5d0*(ql(i,mu)**2
     &           + ql(i,mv)**2 + ql(i,mw)**2) / rhoi)
	 ci = dsqrt(gamma*pi/rhoi)
	 s3 = ql(i,mu)/rhoi + ci     !# u+c in right state  (cell i)
c
         rho2 = ql(i,1) - wave(i,1,3)
         rhou2 = ql(i,mu) - wave(i,mu,3)
         rhov2 = ql(i,mv) - wave(i,mv,3)
         rhow2 = ql(i,mw) - wave(i,mw,3)
         en2 = ql(i,5) - wave(i,5,3)
         p2 = gamma1*(en2 - 0.5d0*(rhou2**2 + rhov2**2 +
     &                rhow2**2)/rho2)
         c2 = dsqrt(gamma*p2/rho2)
         s2 = rhou2/rho2 + c2   !# u+c to left of 3-wave
         if (s2 .lt. 0.d0 .and. s3.gt.0.d0 ) then
c            # transonic rarefaction in the 3-wave
	     sfract = s2 * (s3-s(i,3)) / (s3-s2)
	   else if (s(i,3) .lt. 0.d0) then
c            # 3-wave is leftgoing
	     sfract = s(i,3)
	   else
c            # 3-wave is rightgoing
	     go to 200
	   endif
c
	 do 160 m=1,5
	    amdq(i,m) = amdq(i,m) + sfract*wave(i,m,3)
  160       continue
  200    continue
c
c     # compute the rightgoing flux differences:
c     # df = SUM s*wave   is the total flux difference and apdq = df - amdq
c
      do 220 m=1,meqn
	 do 220 i = 2-mbc, mx+mbc
	    df = 0.d0
	    do 210 mws=1,mwaves
	       df = df + s(i,mws)*wave(i,m,mws)
  210          continue
	    apdq(i,m) = df - amdq(i,m)
  220       continue
c
  900 continue
      return
      end
