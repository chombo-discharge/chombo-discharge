c
c
c =========================================================
       subroutine qinit(maxmx,meqn,mbc,mx,xlower,dx,q,maux,aux)
c =========================================================
c
c     # Set initial conditions for q.
c     # Pulse in pressure, zero velocity
c
c
      implicit double precision (a-h,o-z)
      dimension q(1-mbc:maxmx+mbc, meqn)
      dimension aux(1-mbc:maxmx+mbc, *)

      common /initq/ width
c
      pi = 4.d0*datan(1.d0)

c
      do 150 i=1,mx
	 xcell = xlower + (i-0.5d0)*dx
         r = xcell
c         pressure = exp(-r**2/width**2)/width**2
         if (r .le. 0.5) then
            pressure = 1.d0
         else
            pressure = 0.d0
         endif
	 q(i,1) = pressure
	 q(i,2) = 0.d0
  150    continue
c
      return
      end
