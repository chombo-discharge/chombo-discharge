c     # set auxilliary data for hexahedral mesh
c     ==================================================================
      subroutine sethexinfo(ifstore,iface_store,ivstore,
     &      maxmx, maxmy, maxmz,mbc,mx,my,mz,
     &      xlower,ylower,zlower,dx,dy,dz,
     &      maxlevel,level,refratios,mcapa,maux,aux)
c     ==================================================================
c
      implicit none

c     # Input variables
      integer ifstore, iface_store(3), ivstore, maxmx, maxmy, maxmz,
     &      mx, my, mz, mbc, maux, mcapa
      integer maxlevel, level, refratios(maxlevel)

      double precision xlower, ylower, zlower, dx, dy,dz

      double precision aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc,
     &               1-mbc:maxmz+mbc, maux)


c     # Local variables
      double precision volume, area(3), areafact(3), volfact
      double precision xe,ye,ze, xp,yp,zp, xcorner, ycorner, zcorner
      double precision xef,yef,zef, dxf, dyf, dzf
      double precision sum_volume, sum_face_area(3)
      double precision rot(3,3,3), hex(0:1,0:1,0:1,3)
      integer i,j,k, icell, jcell, kcell, ii, jj, i1, n,m
      integer ir, jr, kr, rfactor, mcapa_test

c     # Common block variables.
c     # these are assigned values here, and never changed outside
c     # of this routine.
      integer ifstore_com, iface_store_com(3), ivstore_com, mcapa_com,
     &      iface_store_rev(3)

      common /comstore/ ifstore_com, iface_store_com, ivstore_com,
     &      mcapa_com, iface_store_rev

c     # ivstore is the number of vectors per face to store.  If we store 2, we
C     # can always recover the third (n x t = b)
c     #
c     # ifstore is the number of faces to store.  This will be what the user
c     # decides - depending on the mapping.
c     #

      mcapa_test = 3*(ivstore*ifstore + 1) + 1
      if (maux .lt. mcapa) then
         write(6,'(A,A,I5,A,I2)') 'sethexinfo : maux must ',
     &         'be at least ',mcapa,' when using sethexinfo.'
         stop
      endif
      if (mcapa .ne. mcapa_test) then
         write(6,'(A,A,I2,A)') 'sethexinfo : mcapa should be set to ',
     &         '3*(ivstore*ifstore+1) + 1 when using sethexinfo.'
         stop
      endif

      rfactor = 1
      do ir = level,maxlevel-1
         rfactor = rfactor*refratios(ir+1)
      enddo
      dxf = dx/rfactor
      dyf = dy/rfactor
      dzf = dz/rfactor

c     # Let's just hardwire the number of vectors per face to store
      areafact(1) = 1.d0/(dy*dz)
      areafact(2) = 1.d0/(dx*dz)
      areafact(3) = 1.d0/(dx*dy)
      volfact = 1.d0/(dx*dy*dz)

      do i = 1-mbc,mx+mbc
         xe = xlower + (i - 1)*dx
	 do j = 1-mbc,my+mbc
            ye = ylower + (j - 1)*dy
            do  k = 1-mbc,mz+mbc
               ze = zlower + (k - 1)*dz
               sum_volume = 0.d0
               do jj = 1,3
                  sum_face_area(jj) = 0.d0
               enddo
               do ir = 1,rfactor
                  xef = xe + (ir-1)*dxf
                  do jr = 1,rfactor
                     yef = ye + (jr-1)*dyf
                     do kr = 1,rfactor
                        zef = ze + (kr-1)*dzf

c                       # Get corners of refined mesh cell
                        do icell = 0,1
                           xcorner = xef + icell*dxf
                           do jcell = 0,1
                              ycorner = yef + jcell*dyf
                              do kcell = 0,1
                                 zcorner = zef + kcell*dzf
                                 call mapc2p(xcorner, ycorner, zcorner,
     &                                 xp,yp,zp)
                                 hex(icell,jcell,kcell,1) = xp
                                 hex(icell,jcell,kcell,2) = yp
                                 hex(icell,jcell,kcell,3) = zp
                              enddo
                           enddo
                        enddo !! hex is complete.

                        call compute_volume(hex,volume)
                        if (volume .le. 0) then
                           write(6,*) 'In sethexinfo : volume <= 0'
                           write(6,*) 'level = ', level
                           write(6,'(A,E24.16)') 'volume = ', volume
                           write(6,*) i,j,k
                           write(6,'(3E24.16)') (hex(0,0,0,m),m=1,3)
                           write(6,'(3E24.16)') (hex(1,1,1,m),m=1,3)
                           stop
                        endif
                        sum_volume = sum_volume + volume

                        if (  ir .eq. 1 .or.
     &                        jr .eq. 1 .or.
     &                        kr .eq. 1) then
                           call compute_surf_area(hex,area)
                        endif

                        if (ir .eq. 1) then
                           sum_face_area(1) = sum_face_area(1) + area(1)
                        endif
                        if (jr .eq. 1) then
                           sum_face_area(2) = sum_face_area(2) + area(2)
                        endif
                        if (kr .eq. 1) then
                           sum_face_area(3) = sum_face_area(3) + area(3)
                        endif
                     enddo
                  enddo
               enddo

               do icell = 0,1
                  xcorner = xe + icell*dx
                  do jcell = 0,1
                     ycorner = ye + jcell*dy
                     do kcell = 0,1
                        zcorner = ze + kcell*dz
                        call mapc2p(xcorner,ycorner,zcorner,xp,yp,zp)
                        hex(icell,jcell,kcell,1) = xp
                        hex(icell,jcell,kcell,2) = yp
                        hex(icell,jcell,kcell,3) = zp
                     enddo
                  enddo
               enddo

               call compute_basis(hex,rot)

c              # Set entries 1-9*nstore with nstore 3x3 rotation matrices
               do i1 = 1,ifstore  !! number of faces stored
                  n = iface_store(i1)  !! face stored
c                 # Only store 2 of the three vectors at each face.
                  do ii = 1,2
c                    # Components of each vector
                     do jj = 1,3
                        ir = (i1 - 1)*3*ivstore + (ii - 1)*3 + jj
                        aux(i,j,k,ir) = rot(n,ii,jj)
                     enddo
                  enddo
               enddo

c              # We have to store all of the face info, since even if, say
c              # the z face is the only face that is mapped, all
c              # three faces (xlow, ylow,zlow) could have areas different
c              # from dx*dy, dx*dz or dy*dz
               do jj = 1,3
                  ir = 3*ivstore*ifstore + jj
                  aux(i,j,k,ir) = sum_face_area(jj)*areafact(jj)
               enddo

c              # Set entry 31 to capacity
               aux(i,j,k,mcapa) = sum_volume*volfact
            enddo
         enddo
      enddo

c     # Store these for recall by later routines.  By storing these
c     # things in a common block, these things can be accesses by routines
c     # which are called from the Riemann solvers without having to be
c     # passed explicitly into those routines.  This SHOULD prevent
c     #  problems for the user.

      ifstore_com = ifstore
      do i = 1,ifstore
         iface_store_com(i) = iface_store(i)
      enddo

      mcapa_com = mcapa
      ivstore_com = ivstore

      do i = 1,ifstore
         n = iface_store(i)
         iface_store_rev(n) = i
      enddo

      return
      end


      subroutine compute_auxlr(ixyz,maxm,mbc,mx,meqn,ql,qr,
     &      xlower,ylower,zlower,dx,dy,dz,t,dt,maux,auxl,auxr)
      implicit none

      integer maxm,mbc,mx,meqn,maux,ixyz
      double precision xlower, ylower, zlower, dx,dy,dz,t,dt

      double precision   ql(1-mbc:maxm+mbc, meqn)
      double precision   qr(1-mbc:maxm+mbc, meqn)
      double precision auxl(1-mbc:maxm+mbc, maux)
      double precision auxr(1-mbc:maxm+mbc, maux)

      double precision hex(0:1,0:1,0:1)
      integer i

      do i = 1-mbc,maxm+mbc
      enddo







      end


      subroutine compute_basis(z,rot)
      implicit none

      double precision z(0:1,0:1,0:1,3)

      double precision Jb(3,3), Jinv(3,3), xcv(3), ycv(3), zcv(3)
      double precision eye(3,3), rot(3,3,3), v(3,3), alpha(3),
     &      rhs(3,3)
      integer IPIV(3)

      double precision z000(3),z100(3),z010(3),z001(3), z110(3),
     &      z101(3), z011(3), z111(3)

      double precision a000(3),a100(3),a010(3),a001(3), a110(3),
     &      a101(3), a011(3), a111(3)

      integer i,j, k, n, info, ni
      double precision sum

      logical is_id

      do i = 1,3
         do j = 1,3
            if (i == j) then
               eye(i,j) = 1.d0
            else
               eye(i,j) = 0.d0
            endif
         enddo
      enddo


c     Get centers of faces so we can find basis vectors at each face.
      do i = 1,3
         xcv(i) = 0.5d0
         ycv(i) = 0.5d0
         zcv(i) = 0.5d0
      enddo

      xcv(1) = 0.d0
      ycv(2) = 0.d0
      zcv(3) = 0.d0

c     # Make notation easy..
      do i = 1,3
         z000(i) = z(0,0,0,i)
         z100(i) = z(1,0,0,i)
         z010(i) = z(0,1,0,i)
         z001(i) = z(0,0,1,i)
         z110(i) = z(1,1,0,i)
         z101(i) = z(1,0,1,i)
         z011(i) = z(0,1,1,i)
         z111(i) = z(1,1,1,i)
      enddo

c     # Get coefficients for trilinear map.
      do i = 1,3
         a000(i) = z000(i)
         a001(i) = z001(i) - z000(i)
         a010(i) = z010(i) - z000(i)
         a100(i) = z100(i) - z000(i)

         a111(i) = z111(i) - z110(i) - z101(i) -
     &         z011(i) + z100(i) + z010(i) + z001(i) - z000(i)
         a011(i) = z011(i) - z010(i) - z001(i) + z000(i)
         a101(i) = z101(i) - z100(i) - z001(i) + z000(i)
         a110(i) = z110(i) - z100(i) - z010(i) + z000(i)
      enddo

c     # Start computing basis vectors.
      do n = 1,3  !! Loop over each face.

         do i = 1,3 !! Get three rows of Jacobian
            Jb(i,1) = a100(i) + a110(i)*ycv(n)  + a101(i)*zcv(n) +
     &            a111(i)*ycv(n)*zcv(n)
            Jb(i,2) = a010(i) + a110(i)*xcv(n)  + a011(i)*zcv(n) +
     &            a111(i)*xcv(n)*zcv(n)
            Jb(i,3) = a001(i) + a101(i)*xcv(n)  + a011(i)*ycv(n) +
     &            a111(i)*xcv(n)*ycv(n)
         enddo

         do i = 1,3
            IPIV(i) = 0
            do j = 1,3
               rhs(i,j) = eye(i,j)
            enddo
         enddo

c        # Compute inverse of Jacobian.
c         call dgesv(3,3,Jb,3,IPIV,rhs,3,info)
         call compute_Jinv(jb,rhs)

c        if (info .ne. 0) then
c           # The Jacobian is singular;  set rotation to identity and hope
c           # that it is never really needed in other clawpack calculations...
c           do i = 1,3
c               ni = mod(i+n-2,3) + 1
c              ni = i
c              do j = 1,3
c                 rot(n,i,j) = eye(ni,j)
c              enddo
c           enddo
c
c           # Skip to next face (cycle through n loop)
c           cycle
c
c        endif

c        # Now we have to permute Jinv to get normal to face n in first row.
C        # We then orthogonalize other rows against this vector.
         do i = 1,3
            ni = mod(i+n-2,3) + 1
            do j = 1,3
               Jinv(i,j) = rhs(ni,j)
            enddo
         enddo

c        # Do Gram-Schmidt on rows of Jinv, to get orthogonal basis.
         do i = 1,3

c           # Project onto each previous vector v(i,:)
            alpha(1) = 0
            do k = 1,i-1
               alpha(k) = 0
               do j = 1,3
                  alpha(k) = alpha(k) + Jinv(i,j)*v(k,j)
               enddo
            enddo

            sum = 0
            do j = 1,3
               v(i,j) = Jinv(i,j)
               do k = 1,i-1
                  v(i,j) = v(i,j) - alpha(k)*v(k,j)
               enddo
               sum = sum + v(i,j)*v(i,j)
            enddo

            do j = 1,3
               v(i,j) = v(i,j)/sqrt(sum)
            enddo
         enddo


c        # We now have our three basis vectors for face n. Now, permute them
C        back so that row 1 is always normal in direction xi, row 2 in
C        direction eta, and row3 in direction zeta.

c        # If we don't do this, then u is always the direction normal to the
C        current face, v is always the first transverse, w the second.  The
C        Riemann solvers as they are currently written assume that this is the
C        case.

c        ## Test to see if we are storing the identity matrix;  if so, we
C        shouldn't store it above.
         is_id = .true.
         do i = 1,3
c            ni = mod(i+n-2,3) + 1
            ni = i
            do j = 1,3
               if (abs(v(i,j) - eye(i,j)) > 1e-12) then
                  is_id = .false.
               endif
               rot(n,ni,j) = v(i,j)
            enddo
         enddo
      enddo

      end


      subroutine compute_volume(cube,volume)
      implicit  none


      double precision dot_cross
      double precision cube(0:1,0:1,0:1,3)
      double precision u(3),v(3),w(3), p(3), q(3), r(3)

      integer i,j,k, ip,jp, kp, m
      double precision volume

      do j = 1,3
         u(j) = cube(0,1,1,j) - cube(0,0,0,j)
         v(j) = cube(1,0,1,j) - cube(0,0,0,j)
         w(j) = cube(1,1,0,j) - cube(0,0,0,j)
      enddo
      volume = dabs(dot_cross(u,v,w))/6.d0

      do j = 1,3
         u(j) = cube(0,0,1,j) - cube(1,1,1,j)
         v(j) = cube(0,1,0,j) - cube(1,1,1,j)
         w(j) = cube(1,0,0,j) - cube(1,1,1,j)
      enddo
      volume = volume + dabs(dot_cross(u,v,w))/6.d0

      do i = 0,1
         do j = 0,1
            do k = 0,1
               ip = mod(i+1,2)
               jp = mod(j+1,2)
               kp = mod(k+1,2)

               do m = 1,3
                  p(m) = cube(ip,j,k,m) - cube(i,j,k,m)
                  q(m) = cube(i,jp,k,m) - cube(i,j,k,m)
                  r(m) = cube(i,j,kp,m) - cube(i,j,k,m)
               enddo

               volume = volume + dabs(dot_cross(p,q,r))/6.d0
            enddo
         enddo
      enddo

      volume = volume/2.d0

      end

      double precision function dot_cross(u,v,w)
      implicit none
      double precision u(3), v(3), w(3)

      dot_cross = u(1)*(v(2)*w(3) - v(3)*w(2)) -
     &            u(2)*(v(1)*w(3) - v(3)*w(1)) +
     &            u(3)*(v(1)*w(2) - v(2)*w(1))

      return
      end


c     # This will give us the exact surface area of each face in the case
c     # the surface is flat.  Otherwise, this seems like it should give a
c     # reasonable answer, but this is merely a conjecture!

      subroutine compute_surf_area(cube,area)
      implicit none

      double precision cube(0:1,0:1,0:1,3)
      double precision w(3), p(3), q(3), area(3)
      double precision quad(0:1,0:1,3)

      integer i,j,m,n, ip,jp
      double precision d

      do n = 1,3
         do i = 0,1
            do j = 0,1
               do m = 1,3
                  if (n == 1) then
                     quad(i,j,m) = cube(0,i,j,m)
                  elseif (n == 2) then
                     quad(i,j,m) = cube(i,0,j,m)
                  else
                     quad(i,j,m) = cube(i,j,0,m)
                  endif
               enddo
            enddo
         enddo

         area(n) = 0
         do i = 0,1
            do j = 0,1
               ip = mod(i+1,2)
               jp = mod(j+1,2)

               do m = 1,3
                  p(m) = quad(ip,j,m) - quad(i,j,m)
                  q(m) = quad(i,jp,m) - quad(i,j,m)
               enddo
               w(1) =   p(2)*q(3) - p(3)*q(2)
               w(2) = -(p(1)*q(3) - p(3)*q(1))
               w(3) =   p(1)*q(2) - p(2)*q(1)

               d = w(1)*w(1) + w(2)*w(2) + w(3)*w(3)
               area(n) = area(n) + dsqrt(d)/2.d0
            enddo
         enddo
         area(n) = area(n)/2.d0

      enddo !! face loop

      end

c     # This subroutine uses Cramer's Rule to compute the inverse of the
c     # Jacobian
      subroutine compute_Jinv(Jac,Jinv)
      implicit none

      double precision Jac(3,3), Jinv(3,3)
      double precision detJ, Jadj(2,2)
      double precision dot_cross, detJtmp, s
      integer k1(3,2), k2(3,2), i,j,k,ii,jj

      data k1 /2, 1, 1, 3, 3, 2/
      data k2 /2, 1, 1, 3, 3, 2/

      detJ  = dot_cross(jac(1,1),jac(1,2),jac(1,3))

      s = 1.d0/detJ
      do j = 1,3
         do i = 1,3
            Jinv(i,j) = s*(Jac(k1(j,1),k2(i,1))*Jac(k1(j,2),k2(i,2)) -
     &            Jac(k1(j,1),k2(i,2))*Jac(k1(j,2),k2(i,1)))
            s = -s
         enddo
      enddo

      end
