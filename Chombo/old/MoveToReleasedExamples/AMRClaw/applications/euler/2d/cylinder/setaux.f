c     ============================================
      subroutine setaux(maxmx,maxmy,mbc,mx,my,xlower,ylower,dx,dy,
     &                  maux,aux)
c     ============================================
      implicit none

      integer maxmx, maxmy, mbc, mx, my, maux
      double precision xlower, ylower, dx, dy
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux)

      integer mcapa

      integer maxlevel, level, refratios(20)
      common /comlevel/maxlevel, level, refratios


      mcapa = 7
      call setquadinfo(maxmx, maxmy, mbc,mx,my,xlower,ylower,
     &      dx,dy, maxlevel,level,refratios, mcapa, maux,aux)

      return
      end
