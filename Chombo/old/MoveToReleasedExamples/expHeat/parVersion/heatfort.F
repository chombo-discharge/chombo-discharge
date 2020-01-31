
      subroutine heatsub2d(
     &     phi,
     &     nlphi0, nhphi0,
     &     nlphi1, nhphi1,
     &     lph,
     &     nllph0, nhlph0,
     &     nllph1, nhlph1,
     &     nlreg, nhreg,
     &     nldom, nhdom,
     &     dt, dx, nu)

      implicit none
      integer nlphi0, nhphi0
      integer nlphi1, nhphi1
      integer nllph0, nhlph0
      integer nllph1, nhlph1
      integer nlreg(2), nhreg(2)
      integer nldom(2), nhdom(2)

      real*8 phi(nlphi0:nhphi0, nlphi1:nhphi1)
      real*8 lph(nllph0:nhlph0, nllph1:nhlph1)
      real*8 dt, dx, nu, lapphi
      integer i,j

c      enforce boundary conditions
c     phi(ghost) = -phi(interior)
      if(nlreg(1) .eq. nldom(1)) then
         do j = nlreg(2), nhreg(2)
            phi(nlreg(1)-1,j) = -phi(nlreg(1),j)
         enddo
      endif
      if(nhreg(1) .eq. nhdom(1)) then
         do j = nlreg(2), nhreg(2)
            phi(nhreg(1)+1,j) = -phi(nhreg(1),j)
         enddo
      endif

      if(nlreg(2) .eq. nldom(2)) then
         do i = nlreg(1), nhreg(1)
            phi(i,nlreg(2)-1) = -phi(i, nlreg(2))
         enddo
      endif
      if(nhreg(2) .eq. nhdom(2)) then
         do i = nlreg(1), nhreg(1)
            phi(i, nhreg(2)+1) = -phi(i, nhreg(2))
         enddo
      endif

c      advance solution
      do j = nlreg(2), nhreg(2)
         do i = nlreg(1), nhreg(1)
            lapphi =
     &           (phi(i+1,j)+phi(i,j+1)
     &           +phi(i-1,j)+phi(i,j-1)
     &           -4.0d0*phi(i,j))/(dx*dx)

            lph(i,j) = lapphi
         enddo
      enddo

      do j = nlreg(2), nhreg(2)
         do i = nlreg(1), nhreg(1)
            phi(i,j) = phi(i,j) + nu*dt*lph(i,j)
         enddo
      enddo

      return
      end

      subroutine heatsub3d(
     &     phi,
     &     nlphi0, nhphi0,
     &     nlphi1, nhphi1,
     &     nlphi2, nhphi2,
     &     lph,
     &     nllph0, nhlph0,
     &     nllph1, nhlph1,
     &     nllph2, nhlph2,
     &     nlreg, nhreg,
     &     nldom, nhdom,
     &     dt, dx, nu)

      implicit none
      integer nlphi0, nhphi0
      integer nlphi1, nhphi1
      integer nlphi2, nhphi2
      integer nllph0, nhlph0
      integer nllph1, nhlph1
      integer nllph2, nhlph2
      integer nlreg(3), nhreg(3)
      integer nldom(3), nhdom(3)
      integer i,j,k
      real*8 phi(nlphi0:nhphi0,nlphi1:nhphi1,nlphi2:nhphi2)
      real*8 lph(nllph0:nhlph0,nllph1:nhlph1,nllph2:nhlph2)
      real*8 dt, dx, nu, lapphi

c     enforce boundary conditions
c     phi(ghost) = -phi(interior)
      if(nlreg(1) .eq. nldom(1)) then
         do k = nlreg(3), nhreg(3)
            do j = nlreg(2), nhreg(2)
               phi(nlreg(1)-1,j,k) = -phi(nlreg(1),j,k)
            enddo
         enddo
      endif

      if(nhreg(1) .eq. nhdom(1)) then
         do k = nlreg(3), nhreg(3)
            do j = nlreg(2), nhreg(2)
               phi(nhreg(1)+1,j,k) = -phi(nhreg(1),j,k)
            enddo
         enddo
      endif

      if(nlreg(2) .eq. nldom(2)) then
         do k = nlreg(3), nhreg(3)
            do i = nlreg(1), nhreg(1)
               phi(i,nlreg(2)-1,k) = -phi(i, nlreg(2),k)
            enddo
         enddo
      endif

      if(nhreg(2) .eq. nhdom(2)) then
         do k = nlreg(3), nhreg(3)
            do i = nlreg(1), nhreg(1)
               phi(i, nhreg(2)+1,k) = -phi(i, nhreg(2),k)
            enddo
         enddo
      endif

      if(nlreg(3) .eq. nldom(3)) then
         do j = nlreg(2), nhreg(2)
            do i = nlreg(1), nhreg(1)
               phi(i,j,nlreg(3)-1) = -phi(i,j,nlreg(3))
            enddo
         enddo
      endif

      if(nhreg(3) .eq. nhdom(3)) then
         do j = nlreg(2), nhreg(2)
            do i = nlreg(1), nhreg(1)
               phi(i,j,nhreg(3)+1) = -phi(i,j,nhreg(3))
            enddo
         enddo
      endif

c      advance solution
      do k = nlreg(3), nhreg(3)
         do j = nlreg(2), nhreg(2)
            do i = nlreg(1), nhreg(1)

               lapphi =
     &              (phi(i+1,j  ,k  )+phi(i-1,j  ,k  )
     &              +phi(i  ,j+1,k  )+phi(i  ,j-1,k  )
     &              +phi(i  ,j  ,k+1)+phi(i  ,j  ,k-1)
     &              -6.0d0*phi(i,j,k))/(dx*dx)

               lph(i,j,k) = lapphi
            enddo
         enddo
      enddo


      do k = nlreg(3), nhreg(3)
         do j = nlreg(2), nhreg(2)
            do i = nlreg(1), nhreg(1)
               phi(i,j,k) = phi(i,j,k) + nu*dt*lph(i,j,k)
            enddo
         enddo
      enddo

      return
      end


