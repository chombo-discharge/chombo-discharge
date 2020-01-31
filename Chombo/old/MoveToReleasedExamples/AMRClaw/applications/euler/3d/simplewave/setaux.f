c     ==================================================================
      subroutine setaux(maxmx,maxmy,maxmz,mbc,mx,my,mz,xlower,ylower,
     &                  zlower,dx,dy,dz,maux,aux)
c     ==================================================================
c
c     # set auxiliary arrays
c
      implicit none

      integer maxmx, maxmy, maxmz, mx, my, mz, mbc, maux
      integer i,j,k
      double precision xlower, ylower, zlower, dx, dy, dz,t

      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, maux)

      integer ifstore, iface_store(3), mcapa, ivstore

      integer maxlevel, level, refratios(20)
      common /comlevel/ maxlevel, level, refratios

      ifstore = 3  !! number of faces to store geometry for
      ivstore = 2  !! number of vectors at each face to store
      iface_store(1) = 1   !! index of each face to store
      iface_store(2) = 2
      iface_store(3) = 3
      mcapa = 22
      call sethexinfo(ifstore,iface_store,ivstore,
     &      maxmx, maxmy, maxmz,mbc,mx,my,mz,
     &      xlower,ylower,zlower,dx,dy,dz,
     &      maxlevel,level,refratios,mcapa,maux,aux)

      return
      end
