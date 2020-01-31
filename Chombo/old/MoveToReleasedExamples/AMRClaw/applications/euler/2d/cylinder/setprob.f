      subroutine setprob
      implicit none

      double precision gamma, gamma1, pi
      double precision rhol, rhor, pl, pr, ul, ur
      double precision dMs, ar, dMr, S3

      common /param/ gamma,gamma1
      common /compi/ pi
      common /indata/ rhol,rhor,pl,pr,ul,ur
c
      open(unit=7,file='setprob.data',status='old',form='formatted')

      pi = 4.d0*datan(1.d0)

      read(7,*) gamma
      gamma1 = gamma-1.d0

c     read Mach number and state on right hand site from setprob.data
      read(7,*) dMs
      read(7,*) rhor
      read(7,*) ur
      read(7,*) pr
      close(7)
c
c     calculation of the states behind the shock (for a 3-shock wave)
c
      ar  = dsqrt(gamma*pr/rhor)
      dMr = ur/ar
      rhol = rhor*(gamma+1.d0)*((dMr-dMs)**2.d0) /
     &     ((gamma-1.d0)*((dMr-dMs)**2.d0)+2.d0)
      S3 = dMs*ar
      ul = (1.d0-rhor/rhol)*S3 + ur*rhor/rhol
      pl = pr*(2.d0*gamma*(dMr-dMs)**2.d0-(gamma-1.d0))/(gamma+1.d0)

      write(6,100) 'rho_r = ', rhor
      write(6,100) 'u_r   = ', ur
      write(6,100) 'p_r   = ', pr
      write(6,*) ' '
      write(6,100) 'rho_l = ', rhol
      write(6,100) 'u_l   = ', ul
      write(6,100) 'p_l   = ', pl
      write(6,*) 'shockspeed = ', S3

 100  format(A,F16.8)


      return
      end
