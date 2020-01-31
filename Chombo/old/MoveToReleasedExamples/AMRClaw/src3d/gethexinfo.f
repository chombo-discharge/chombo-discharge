c     # This gets the aux locations needed for the normal Riemann
c     # solve
      subroutine get_aux_locations_n(ixyz,mcapa,locrot,locarea)
      implicit none

      integer ixyz, mcapa,locrot, locarea

      integer ifstore, iface_store(3), ivstore, mcapa_com,
     &      iface_store_rev(3), n

      common /comstore/ ifstore, iface_store, ivstore, mcapa_com,
     &      iface_store_rev

      mcapa = mcapa_com
      n = iface_store_rev(ixyz)
      locarea = 3*ifstore*ivstore + ixyz
      locrot  = 3*(n-1)*ivstore + 1
      end

c     # This gets the aux locations needed for the transverse
c     # and double transverse Riemann solves
      subroutine get_aux_locations(ixyz,icoor,mcapa,
     &      locrot,locarea)
      implicit none

      integer ixyz, icoor, locrot, locarea, impt
      integer ifstore, iface_store(3), ivstore,mcapa_com,
     &      iface_store_rev(3), n, nxyz, mcapa

      common /comstore/ ifstore, iface_store,ivstore, mcapa_com,
     &      iface_store_rev


      mcapa = mcapa_com
      if (ixyz .eq. 1) then
         if (icoor .eq. 2) then
c           # y-like direction is a y-face
            nxyz = 2
         elseif (icoor .eq. 3) then
c           # z-like direction is a z-face
            nxyz = 3
         endif
      elseif (ixyz .eq. 2) then
         if (icoor .eq. 2) then
c           # y-like direction is a z-face
            nxyz = 3
         else
c           # y-like direction is an x-face
            nxyz = 1
         endif
      elseif (ixyz .eq. 3) then
         if (icoor .eq. 2) then
c           # y-like direction is an x-face
            nxyz = 1
         elseif (icoor .eq. 3) then
c           # z-like direction is a y-face
            nxyz = 2
         endif
      endif

      n = iface_store_rev(nxyz)
      locarea = 3*ivstore*ifstore + nxyz
      locrot = 3*(n-1)*ivstore + 1

      end

      subroutine rotate3(rot,velcomps)
      implicit  none

      double precision velcomps(3), rot(9)
      double precision v1, v2, v3

      v1 = velcomps(1)
      v2 = velcomps(2)
      v3 = velcomps(3)

      velcomps(1) = rot(1)*v1 + rot(2)*v2 + rot(3)*v3
      velcomps(2) = rot(4)*v1 + rot(5)*v2 + rot(6)*v3
      velcomps(3) = rot(7)*v1 + rot(8)*v2 + rot(9)*v3
      end

      subroutine rotate3_tr(rot,velcomps)
      implicit none

      double precision velcomps(3),rot(9)
      double precision v1, v2, v3

      v1 = velcomps(1)
      v2 = velcomps(2)
      v3 = velcomps(3)

      velcomps(1) = rot(1)*v1 + rot(4)*v2 + rot(7)*v3
      velcomps(2) = rot(2)*v1 + rot(5)*v2 + rot(8)*v3
      velcomps(3) = rot(3)*v1 + rot(6)*v2 + rot(9)*v3
      end


      subroutine compute_binormal(rot)
      implicit none

      double precision rot(9)
      double precision w1,w2,w3,r

      w1 =  rot(2)*rot(6) - rot(3)*rot(5)
      w2 = -rot(1)*rot(6) + rot(3)*rot(4)
      w3 =  rot(1)*rot(5) - rot(2)*rot(4)
      r = dsqrt(w1*w1 + w2*w2 + w3*w3)

      rot(7) = w1/r
      rot(8) = w2/r
      rot(9) = w3/r

      end
