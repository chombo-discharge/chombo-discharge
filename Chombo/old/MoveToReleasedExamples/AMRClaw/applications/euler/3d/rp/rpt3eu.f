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
      integer ixyz, icoor, maxm, meqn, mwaves,mbc,mx,maux,ilr

      double precision      ql(1-mbc:maxm+mbc, meqn)
      double precision      qr(1-mbc:maxm+mbc, meqn)
      double precision    asdq(1-mbc:maxm+mbc, meqn)
      double precision  bmasdq(1-mbc:maxm+mbc, meqn)
      double precision  bpasdq(1-mbc:maxm+mbc, meqn)
      double precision    aux1(1-mbc:maxm+mbc, maux, 3)
      double precision    aux2(1-mbc:maxm+mbc, maux, 3)
      double precision    aux3(1-mbc:maxm+mbc, maux, 3)
c
      double precision waveb(5,3),sb(3)
      double precision a1,a2,a3,a4,a5
      integer mu,mv,mw, mws, i, m

      integer maxmrp
      parameter (maxmrp = 502)
      double precision u2v2w2(-1:maxmrp), u(-1:maxmrp), v(-1:maxmrp)
      double precision w(-1:maxmrp), enth(-1:maxmrp), a(-1:maxmrp)
      double precision g1a2(-1:maxmrp), euv(-1:maxmrp)

      common /comroe/ u2v2w2,u,v,w,enth,a,g1a2,euv

      if (-3.gt.1-mbc .or. maxmrp .lt. maxm+mbc) then
	 write(6,*) 'need to increase maxmrp in rp3t'
	 stop
	 endif
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
c     # Solve Riemann problem in the second coordinate direction
c
      if( icoor .eq. 2 )then

	 do 20 i = 2-mbc, mx+mbc
            a4 = g1a2(i) * (euv(i)*asdq(i,1)
     &             + u(i)*asdq(i,mu) + v(i)*asdq(i,mv)
     &             + w(i)*asdq(i,mw) - asdq(i,5))
	    a2 = asdq(i,mu) - u(i)*asdq(i,1)
            a3 = asdq(i,mw) - w(i)*asdq(i,1)
	    a5 = (asdq(i,mv) + (a(i)-v(i))*asdq(i,1) - a(i)*a4)
     &              / (2.d0*a(i))
	    a1 = asdq(i,1) - a4 - a5
c
            waveb(1,1)  = a1
            waveb(mu,1) = a1*u(i)
            waveb(mv,1) = a1*(v(i)-a(i))
            waveb(mw,1) = a1*w(i)
            waveb(5,1)  = a1*(enth(i) - v(i)*a(i))
	    sb(1) = v(i) - a(i)
c
            waveb(1,2)  = a4
            waveb(mu,2) = a2 + u(i)*a4
            waveb(mv,2) = v(i)*a4
            waveb(mw,2) = a3 + w(i)*a4
            waveb(5,2)  = a4*0.5d0*u2v2w2(i) + a2*u(i) + a3*w(i)
	    sb(2) = v(i)
c
            waveb(1,3)  = a5
            waveb(mu,3) = a5*u(i)
            waveb(mv,3) = a5*(v(i)+a(i))
            waveb(mw,3) = a5*w(i)
            waveb(5,3)  = a5*(enth(i)+v(i)*a(i))
	    sb(3) = v(i) + a(i)
c
	 do 25 m=1,meqn
	    bmasdq(i,m) = 0.d0
	    bpasdq(i,m) = 0.d0
	    do 25 mws=1,mwaves
	       bmasdq(i,m) = bmasdq(i,m)
     &			   + dmin1(sb(mws), 0.d0) * waveb(m,mws)
	       bpasdq(i,m) = bpasdq(i,m)
     &			   + dmax1(sb(mws), 0.d0) * waveb(m,mws)
 25         continue
c
   20    continue
c
      else
c
c        # Solve Riemann problem in the third coordinate direction
c
	 do 30 i = 2-mbc, mx+mbc
            a4 = g1a2(i) * (euv(i)*asdq(i,1)
     &             + u(i)*asdq(i,mu) + v(i)*asdq(i,mv)
     &             + w(i)*asdq(i,mw) - asdq(i,5))
	    a2 = asdq(i,mu) - u(i)*asdq(i,1)
            a3 = asdq(i,mv) - v(i)*asdq(i,1)
	    a5 = (asdq(i,mw) + (a(i)-w(i))*asdq(i,1) - a(i)*a4)
     &              / (2.d0*a(i))
	    a1 = asdq(i,1) - a4 - a5
c
            waveb(1,1)  = a1
            waveb(mu,1) = a1*u(i)
            waveb(mv,1) = a1*v(i)
            waveb(mw,1) = a1*(w(i) - a(i))
            waveb(5,1)  = a1*(enth(i) - w(i)*a(i))
	    sb(1) = w(i) - a(i)
c
            waveb(1,2)  = a4
            waveb(mu,2) = a2 + u(i)*a4
            waveb(mv,2) = a3 + v(i)*a4
            waveb(mw,2) = w(i)*a4
            waveb(5,2)  = a4*0.5d0*u2v2w2(i) + a2*u(i) + a3*v(i)
	    sb(2) = w(i)
c
            waveb(1,3)  = a5
            waveb(mu,3) = a5*u(i)
            waveb(mv,3) = a5*v(i)
            waveb(mw,3) = a5*(w(i)+a(i))
            waveb(5,3)  = a5*(enth(i)+w(i)*a(i))
	    sb(3) = w(i) + a(i)
c
	 do 35 m=1,meqn
	    bmasdq(i,m) = 0.d0
	    bpasdq(i,m) = 0.d0
	    do 35 mws=1,mwaves
	       bmasdq(i,m) = bmasdq(i,m)
     &			   + dmin1(sb(mws), 0.d0) * waveb(m,mws)
	       bpasdq(i,m) = bpasdq(i,m)
     &			   + dmax1(sb(mws), 0.d0) * waveb(m,mws)
 35         continue
c
 30      continue
c
      endif
c
      return
      end
