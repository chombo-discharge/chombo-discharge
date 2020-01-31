      subroutine setquadinfo(maxmx, maxmy, mbc,mx,my,xlower,ylower,
     &      dx,dy,maxlevel,level,refratios,mcapa,maux,aux)

      implicit none

      integer maxmx, maxmy, mbc, mx, my, maux ,mcapa
      double precision xlower, ylower, dx, dy
      double precision aux(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc,
     &      maux)
      double precision normals(2,2), quad(0:1,0:1,2), slength(2)

      double precision sum_length(2)
      double precision dxf, dyf, xef, yef, xe, ye, xp, yp
      double precision sum_area, xcorner, ycorner, area
      integer rfactor, i, j, ir, jr, icell, jcell
      integer maxlevel, level

c     refratios should have at least maxlevel entries.
      integer refratios(maxlevel)

      if (maux .lt. mcapa) then
         write(6,*) 'setquadinfo : maux must be >= 7 when using '
         write(6,*) 'mapped grids with routine setquadinfo.f'
         stop
      endif
      if (mcapa .ne. 7) then
         write(6,'(A,A,I2,A)') 'setquadinfo : mcapa should ',
     &         'be set to 7 when using setquadinfo.;  ',mcapa
         stop
      endif

      rfactor = 1
      do ir = level,maxlevel-1
         rfactor = rfactor*refratios(ir+1)
      enddo
      dxf = dx/rfactor
      dyf = dy/rfactor

      do j = 1-mbc,my+mbc
         ye = ylower + (j-1)*dy
         do i = 1-mbc,mx+mbc
            xe = xlower + (i-1)*dx

            sum_area = 0.d0
            sum_length(1) = 0.d0
            sum_length(2) = 0.d0
            do ir = 1,rfactor
               xef = xe + (ir-1)*dxf
               do jr = 1,rfactor
                  yef = ye + (jr - 1)*dyf
                  do icell = 0,1
                     xcorner = xef + icell*dxf
                     do jcell = 0,1
                        ycorner = yef + jcell*dyf
                        call mapc2p(xcorner,ycorner,xp,yp)
                        quad(icell,jcell,1) = xp
                        quad(icell,jcell,2) = yp
                     enddo
                  enddo
c                 # slow version
c                  call compute_area2(quad,slength,area)

c                 # Fast version
                  call compute_info2(quad,normals,slength,area)
                  sum_area = sum_area + area
                  if (ir .eq. 1) then
                     sum_length(1) = sum_length(1) + slength(1)
                  endif
                  if (jr .eq. 1) then
                     sum_length(2) = sum_length(2) + slength(2)
                  endif
               enddo
            enddo

c           # Now get rotation vector at interface
c           # Use coarse grid straight edge approximation
            do icell = 0,1
               xcorner = xe + icell*dx
               do jcell = 0,1
                  ycorner = ye + jcell*dy
                  call mapc2p(xcorner,ycorner,xp,yp)
                  quad(icell,jcell,1) = xp
                  quad(icell,jcell,2) = yp
               enddo
            enddo

c          # Slow version that calls 3d routine
c           call compute_basis2(quad,rot)

c          # Fast version
           call compute_info2(quad,normals,slength,area)

c          # Normal at x face - first row in normals matrix
           aux(i,j,1) = normals(1,1)
           aux(i,j,2) = normals(1,2)
           aux(i,j,3) = sum_length(1)/dy

c          # Normal at y face - second row in normals matrix.
           aux(i,j,4) = normals(2,1)
           aux(i,j,5) = normals(2,2)
           aux(i,j,6) = sum_length(2)/dx

c          Capacity
           if (area .le. 0.d0) then
              write(6,*) 'setquadinfo : capacity <= 0'
              write(6,'(A,A)') 'This usually means that there ',
     &              'is a problem with your mapping function'
              stop
           endif
           aux(i,j,mcapa) = sum_area/(dx*dy)
        enddo
      enddo

      end

      subroutine compute_info2(quad,normals,slength,area)
      implicit none

      double precision normals(2,2),quad(0:1,0:1,2), slength(2)
      double precision v(2), vnorm
      double precision w(2), wnorm
      double precision area, xpcorn(5), ypcorn(5)
      integer ic

c     # compute normals to left edge
      v(1) = quad(0,1,1) - quad(0,0,1)
      v(2) = quad(0,1,2) - quad(0,0,2)
      vnorm = sqrt(v(1)*v(1) + v(2)*v(2))
      slength(1) = vnorm
      if (vnorm .eq. 0.d0) then
         normals(1,1) = 1.d0
         normals(1,2) = 0.d0
      else
         normals(1,1) = v(2)/vnorm  !! Normal should point in positive x dir.
         normals(1,2) = -v(1)/vnorm
      endif
c
      w(1) = quad(1,0,1) - quad(0,0,1)
      w(2) = quad(1,0,2) - quad(0,0,2)
      wnorm = sqrt(w(1)*w(1) + w(2)*w(2))
      slength(2) = wnorm
      if (wnorm .eq. 0.d0) then
         normals(2,1) = 1.d0
         normals(2,2) = 0.d0
      else
         normals(2,1) = -w(2)/wnorm  !! normal should point in positive y dir.
         normals(2,2) = w(1)/wnorm
      endif

      xpcorn(1) = quad(0,0,1)
      xpcorn(2) = quad(0,1,1)
      xpcorn(3) = quad(1,1,1)
      xpcorn(4) = quad(1,0,1)
      xpcorn(5) = quad(0,0,1)

      ypcorn(1) = quad(0,0,2)
      ypcorn(2) = quad(0,1,2)
      ypcorn(3) = quad(1,1,2)
      ypcorn(4) = quad(1,0,2)
      ypcorn(5) = quad(0,0,2)

      area = 0.d0
      do ic=1,4
         area = area + 0.5d0 * (ypcorn(ic)+ypcorn(ic+1)) *
     &         (xpcorn(ic+1)-xpcorn(ic))
      enddo

      end
