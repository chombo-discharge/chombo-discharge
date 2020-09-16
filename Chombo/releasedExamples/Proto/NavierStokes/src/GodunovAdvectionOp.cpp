#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include "GodunovAdvectionOp.H"
#include "Proto.H"
#include "Proto_Timer.H"

using std::cout;

namespace Proto
{

  void
  GodunovAdvectionOp::
  divUPhi(BoxData<double,1>          & a_div,
          const BoxData<double,   1>   a_vel[DIM],       // 
          const BoxData<double,   1>   a_phi[DIM],
          const Box                   & a_region) const
  {
    PR_TIME("gao::divuphi");
    a_div.setVal(0.);
    for(int faceDir = 0; faceDir < DIM; faceDir++)
    {
      Box Bface = a_region.extrude(faceDir);
      BoxData<double, 1> flux(Bface);
      a_vel[faceDir].copyTo(flux);
      flux *= a_phi[faceDir];
      //increment the divergence with the flux difference
      a_div += m_fluxDivIncr[faceDir](flux, a_region);
    }

  }

  GodunovAdvectionOp::
  GodunovAdvectionOp(double a_dx)
  {
    m_dx = a_dx;
    PR_TIME("gao::gao");
    for (int dir = 0; dir < DIM; dir++)
    {
      Point hidir =  Point::Basis(dir,  1.0);
      Point lodir =  Point::Basis(dir, -1.0);
      Point zedir =  Point::Zeros();
    
      m_slopeLo[dir] =  (1.0)*Shift(zedir) + (-1.0)*Shift(lodir);
      m_slopeHi[dir] =  (1.0)*Shift(hidir) + (-1.0)*Shift(zedir);
      m_faceToCell[dir] =(0.5)*(1.0*Shift(lodir) + 1.0*Shift(zedir));
      m_fluxDivIncr[dir] = (1.0/m_dx)*Shift(hidir) + (-1.0/m_dx)*Shift(zedir);
    }
  }



  PROTO_KERNEL_START
  unsigned int 
  upwindSlopeF(Scalar& a_out,
               const Scalar& a_lo,
               const Scalar& a_hi,
               const Pvector& a_vel,
               int            a_idir)
  {
    if(a_vel(a_idir) > 0) 
    {
      a_out(0) = a_lo(0);
    }
    else if (a_vel(a_idir) < 0) 
    {
      a_out(0) = a_hi(0);
    }
    else
    {
      a_out(0) = 0.0;
    }
//  a_out(0) = 0.0;
    return 0;
  }
  PROTO_KERNEL_END(upwindSlopeF, upwindSlope)

#ifdef PROTO_CUDA
  __device__
#endif
  double minVal(double vala, double valb)
  {
    double minval;
    if(vala < valb)
    {
      minval = vala;
    }
    else
    {
      minval = valb;
    }
    return minval;
  }

#ifdef PROTO_CUDA
  __device__
#endif
  double absVal(double vala)
  {
    double absval;
    if(vala < 0.0)
    {
      absval = -vala;
    }
    else
    {
      absval =  vala;
    }
    return absval;
  }


#ifdef PROTO_CUDA
  __device__
#endif
  double cpySign(double vala, double signval)
  {
    double retval;
    if(signval < 0.0)
    {
      retval = -vala;
    }
    else
    {
      retval =  vala;
    }
    return retval;
  }

  PROTO_KERNEL_START
  unsigned int
  vanLeerSlopeF(Scalar& a_slope,
                const Scalar& a_lo,
                const Scalar& a_hi)
  {
    double mono = a_lo(0)*a_hi(0);
    if(mono < 0) 
    {
      a_slope(0) = 0.0;
    }
    else
    {
      double cen = 0.5*(a_lo(0) + a_hi(0));
     
      a_slope(0) = minVal(absVal(cen)       , absVal(2.*a_lo(0)));
      a_slope(0) = minVal(absVal(a_slope(0)), absVal(2.*a_hi(0)));
      a_slope(0) = cpySign(a_slope(0),cen);
//    a_slope(0) = std::min(std::abs(cen)       , std::abs(2.*a_lo(0)));
//    a_slope(0) = std::min(std::abs(a_slope(0)), std::abs(2.*a_hi(0)));
//    a_slope(0) = std::copysign(a_slope(0),cen);


    }
//  a_slope(0) = 0.0;
    return 0;
  }
  PROTO_KERNEL_END(vanLeerSlopeF, vanLeerSlope)

  PROTO_KERNEL_START
  unsigned int
  addNormalDerivativesF(Scalar        & a_phiFace,
                        const Scalar  & a_vanLeerSlope,
                        const Pvector  & a_vel,
                        double   a_dx,
                        double   a_dt,
                        int      a_slopeSign,
                        int      a_idir)
  {
    double extrap = 0.5*(a_slopeSign - (a_dt/a_dx)*a_vel(a_idir))*a_vanLeerSlope(0);
    a_phiFace(0) += extrap;
    return 0;
  }
  PROTO_KERNEL_END(addNormalDerivativesF, addNormalDerivatives)


  PROTO_KERNEL_START
  unsigned int
  addTangentialDerivativesF(Scalar        & a_phiFace,
                            const Scalar  & a_upwindSlope,
                            const Pvector  & a_vel,
                            double   a_dx,
                            double   a_dt,
                            int      a_idir)
  {
    double extrap = 0.5*a_vel(a_idir)*(-a_dt/a_dx)*a_upwindSlope(0);
    a_phiFace(0) += extrap;
    return 0;
  }
  PROTO_KERNEL_END(addTangentialDerivativesF, addTangentialDerivatives)

  PROTO_KERNEL_START
  unsigned int 
  getPhiHF(Scalar       & a_phiH, 
           const Scalar & a_phiLeft, 
           const Scalar & a_phiRigh, 
           const Scalar & a_vel,
           bool  a_doingVel)
  {
    //choose upwind state
    double phiup ;
    const double& phil = a_phiLeft(0);
    const double& phir = a_phiRigh(0);

    if (a_vel(0) > 0) 
    {
      phiup = phil;
    }
    else if (a_vel(0) < 0) 
    {
      phiup = phir;
    }
    else if(!a_doingVel)
    {
      phiup = 0.5*(phil + phir);
    }
    else
    {
      phiup = 0.0;
    }

    a_phiH(0) = phiup;
    return 0;
  }
  PROTO_KERNEL_END(getPhiHF, getPhiH)


  PROTO_KERNEL_START
  unsigned int 
  getFluxF(Scalar       & a_F, 
           const Scalar & a_phiLeft, 
           const Scalar & a_phiRigh, 
           const Scalar & a_vel,
           bool    a_doingVel)
  {
    //choose upwind state
    double phiup ;
    const double& phil = a_phiLeft(0);
    const double& phir = a_phiRigh(0);
    const double& velo = a_vel(0);
    if (a_vel(0) > 0) 
    {
      phiup = phil;
    }
    else if (a_vel(0) < 0) 
    {
      phiup = phir;
    }
    else if(!a_doingVel)
    {
      phiup = 0.5*(phil + phir);
    }
    else
    {
      phiup = 0.0;
    }
    a_F(0) = phiup*velo;
//  cout << "phiup = " << phiup << ", velo = " << velo << ", flux = "<< a_F(0)<< std::endl;
    return 0;
  }
  PROTO_KERNEL_END(getFluxF, getFlux)

  void
  GodunovAdvectionOp::
  divFluxNPH(BoxData<double,1>          & a_div,
             const BoxData<double,   1>   a_vel[DIM],
             const BoxData<double,   1> & a_phi,
             const BoxData<double,   1> & a_src,
             const BoxData<double, DIM> & a_velCell,
             const Box                  & a_dbx0,
             const int                  & a_doingVel,
             const double               & a_dt) const
  {
    PR_TIME("GodunovAdvectionOp::operator");

    bool doingVel = (a_doingVel == 1);
    double dx = m_dx;
    double dt = a_dt;
    Box B_0 = a_dbx0;
    Box B_1 = B_0.grow(1);
    Box B_2 = B_0.grow(2);
    //face centered boxes
    Box Bface[DIM];
    for(int idir = 0; idir < DIM; idir++)
    {
      Bface[idir] = B_0.extrude(idir);
    }

    const BoxData<double, DIM>& cellVel = a_velCell;

    //get upwind and van leer slopes from one sided slopes.
    BoxData<double,1> slopeUpwind[DIM];
    BoxData<double,1> slopeVLimit[DIM];
    for(int idir = 0; idir < DIM; idir++)
    {
      BoxData<double,1> slopeLo(B_2);
      BoxData<double,1> slopeHi(B_2);


      slopeLo |= m_slopeLo[idir](a_phi, B_2);
      slopeHi |= m_slopeHi[idir](a_phi, B_2);

      slopeUpwind[idir].define(B_1);
      slopeVLimit[idir].define(B_1);

      forallInPlace(upwindSlope, B_1, slopeUpwind[idir], slopeLo, slopeHi, cellVel, idir);
      forallInPlace(vanLeerSlope, B_1, slopeVLimit[idir], slopeLo, slopeHi);
    }


    //to avoid too many temporaries, this is a big loop
    int signLo = -1;
    int signHi =  1;
    //computing the divergence  by increment so we have to set it to zero
    a_div.setVal(0.0);
    for(int faceDir = 0; faceDir < DIM; faceDir++)
    {
      BoxData<double, 1> flux(Bface[faceDir]);
      BoxData<double, 1> phiLo(B_1);
      BoxData<double, 1> phiHi(B_1);

      a_phi.copyTo(phiHi, B_1);
      a_phi.copyTo(phiLo, B_1);
      for(int jdir = 0; jdir < DIM; jdir++)
      {
        if(jdir == faceDir)
        {
          forallInPlace(addNormalDerivatives,B_1,
                        phiLo,slopeVLimit[jdir],cellVel,dx,dt,signLo,jdir);
          forallInPlace(addNormalDerivatives,B_1,
                        phiHi,slopeVLimit[jdir],cellVel,dx,dt,signHi,jdir);
        }
        else
        {
          forallInPlace(addTangentialDerivatives,B_1,
                        phiLo, slopeUpwind[jdir], cellVel, dx, dt, jdir);
          forallInPlace(addTangentialDerivatives,B_1,
                        phiHi, slopeUpwind[jdir], cellVel, dx, dt, jdir);
        }
      }
      //now shift phi lo and high to make left and right states on the face.
      //philo was extrapolated in the low  direction so it becomes the right state (and its index space is correct)
      //phihi was extrapolated in the high direction so it becomes the left  state (and its index space needs to be shifted up)

      BoxData<double, 1>& phiRigh = phiLo;
      BoxData<double, 1>& phiLeft = phiHi;
      phiLeft.shift(Point::Basis(faceDir));

      //now we can solve the riemann problem and get the flux all at once.
      forallInPlace(getFlux, Bface[faceDir], flux, phiLeft, phiRigh, a_vel[faceDir], doingVel);
      /*
        forall([this, doingVel](Scalar& flux, const Scalar&  phil, const Scalar& phir, const Scalar& vel)
        {return getFlux(flux, phil, phir, vel, doingVel);}, 
        Bface[faceDir], flux, phiLeft, phiRigh, a_vel[faceDir]);
      */
    
      //cout << faceDir << ", flux max  = " << flux.max() << ", flux min = " << flux.min() << std::endl;
      //cout << faceDir << ", phil max  = " << phiLeft.max() << ", phil min = " << phiLeft.min() << std::endl;
      //cout << faceDir << ", phir max  = " << phiRigh.max() << ", phir min = " << phiRigh.min() << std::endl;
      //cout << faceDir << ", vel  max  = " << a_vel[faceDir].max() << ", vel  min = " << a_vel[faceDir].min() << std::endl;
      //increment the divergence with the flux difference
      a_div += m_fluxDivIncr[faceDir](flux, B_0);
    }

    //finally we need to take the divergence of the flux.
    return;
  }

//same as above but without the divergence
  void
  GodunovAdvectionOp::
  getFluxNPH(BoxData<double,1>            a_flux[DIM],
             const BoxData<double,   1>   a_vel[DIM],
             const BoxData<double,   1> & a_phi,
             const BoxData<double,   1> & a_src,
             const BoxData<double, DIM> & a_velCell,
             const Box                  & a_dbx0,
             const int                  & a_doingVel,
             const double               & a_dt) const
  {
    PR_TIME("gao::getfluxnph");

    bool doingVel = (a_doingVel == 1);
    double dx = m_dx;
    double dt = a_dt;
    Box B_0 = a_dbx0;
    Box B_1 = B_0.grow(1);
    Box B_2 = B_0.grow(2);
    //face centered boxes
    Box Bface[DIM];
    for(int idir = 0; idir < DIM; idir++)
    {
      Bface[idir] = B_0.extrude(idir);
    }

    const BoxData<double, DIM>& cellVel = a_velCell;

    //get upwind and van leer slopes from one sided slopes.
    BoxData<double,1> slopeUpwind[DIM];
    BoxData<double,1> slopeVLimit[DIM];
    for(int idir = 0; idir < DIM; idir++)
    {
      BoxData<double,1> slopeLo(B_2);
      BoxData<double,1> slopeHi(B_2);


      slopeLo |= m_slopeLo[idir](a_phi, B_2);
      slopeHi |= m_slopeHi[idir](a_phi, B_2);

      slopeUpwind[idir].define(B_1);
      slopeVLimit[idir].define(B_1);

      forallInPlace(upwindSlope, B_1,
                    slopeUpwind[idir], slopeLo, slopeHi, cellVel, idir);
      forallInPlace(vanLeerSlope, B_1,
                    slopeVLimit[idir], slopeLo, slopeHi);
    }


    //to avoid too many temporaries, this is a big loop
    int signLo = -1;
    int signHi =  1;

    for(int faceDir = 0; faceDir < DIM; faceDir++)
    {
      BoxData<double, 1>& flux = a_flux[faceDir];
      BoxData<double, 1> phiLo(B_1);
      BoxData<double, 1> phiHi(B_1);

      a_phi.copyTo(phiHi, B_1);
      a_phi.copyTo(phiLo, B_1);
      for(int jdir = 0; jdir < DIM; jdir++)
      {
        if(jdir == faceDir)
        {
          forallInPlace(addNormalDerivatives, B_1,
                        phiLo, slopeVLimit[jdir], cellVel, dx, dt, signLo, jdir);
          forallInPlace(addNormalDerivatives, B_1,
                        phiHi, slopeVLimit[jdir], cellVel, dx, dt, signHi, jdir);
        }
        else
        {
          forallInPlace(addTangentialDerivatives, B_1,
                        phiLo, slopeUpwind[jdir], cellVel, dx, dt, jdir);
          forallInPlace(addTangentialDerivatives, B_1,
                        phiHi, slopeUpwind[jdir], cellVel, dx, dt, jdir);
        }
      }
      //now shift phi lo and high to make left and right states on the face.
      //philo was extrapolated in the low  direction so it becomes the right state (and its index space is correct)
      //phihi was extrapolated in the high direction so it becomes the left  state (and its index space needs to be shifted up)

      BoxData<double, 1>& phiRigh = phiLo;
      BoxData<double, 1>& phiLeft = phiHi;
      phiLeft.shift(Point::Basis(faceDir));
    
      //now we can solve the riemann problem and get the flux all at once.
    
      forallInPlace(getFlux, Bface[faceDir],
                    flux, phiLeft, phiRigh, a_vel[faceDir], doingVel);
      /*
        forall([this, doingVel](Scalar& flux, const Scalar&  phil, const Scalar& phir, const Scalar& vel)
        {return getFlux(flux, phil, phir, vel, doingVel);}, 
        Bface[faceDir], flux, phiLeft, phiRigh, a_vel[faceDir]);

      */
    }

    //finally we need to take the divergence of the flux.
    return;
  }



//same as above but without multiplying in the velocity
  void
  GodunovAdvectionOp::
  advectToFaces(BoxData<double,1>            a_phiHalf[DIM],
                const BoxData<double,   1>   a_vel[DIM],
                const BoxData<double,   1> & a_phi,
                const BoxData<double,   1> & a_src,
                const BoxData<double, DIM> & a_velCell,
                const Box                  & a_dbx0,
                const int                  & a_doingVel,
                const double               & a_dt) const
  {
    PR_TIME("gao::advecttofaces");

    bool doingVel = (a_doingVel == 1);

    double dx = m_dx;
    double dt = a_dt;
    Box B_0 = a_dbx0;
    Box B_1 = B_0.grow(1);
    Box B_2 = B_0.grow(2);
    //face centered boxes
    Box Bface[DIM];
    for(int idir = 0; idir < DIM; idir++)
    {
      Bface[idir] = B_0.extrude(idir);
    }

    const BoxData<double, DIM>& cellVel = a_velCell;

    //get upwind and van leer slopes from one sided slopes.
    BoxData<double,1> slopeUpwind[DIM];
    BoxData<double,1> slopeVLimit[DIM];
    for(int idir = 0; idir < DIM; idir++)
    {
      BoxData<double,1> slopeLo(B_2);
      BoxData<double,1> slopeHi(B_2);


      slopeLo |= m_slopeLo[idir](a_phi, B_2);
      slopeHi |= m_slopeHi[idir](a_phi, B_2);

      slopeUpwind[idir].define(B_1);
      slopeVLimit[idir].define(B_1);

      forallInPlace(upwindSlope, B_1,
                    slopeUpwind[idir], slopeLo, slopeHi, cellVel, idir);
      forallInPlace(vanLeerSlope, B_1,
                    slopeVLimit[idir], slopeLo, slopeHi);
      /*
        forall([this,idir](Scalar& upws, const Scalar& losl, const Scalar& hisl, const Pvector& vel)
        {return upwindSlope(upws, losl, hisl, vel, idir);}, 
        B_1, slopeUpwind[idir], slopeLo, slopeHi, cellVel);

        forall([this,idir](Scalar& vlsl, const Scalar& losl, const Scalar& hisl)
        {return vanLeerSlope(vlsl, losl, hisl);}, 
        B_1, slopeVLimit[idir], slopeLo, slopeHi);
      */
    }


    //to avoid too many temporaries, this is a big loop
    int signLo = -1;
    int signHi =  1;

    for(int faceDir = 0; faceDir < DIM; faceDir++)
    {
      BoxData<double, 1> phiLo(B_1);
      BoxData<double, 1> phiHi(B_1);

      a_phi.copyTo(phiHi, B_1);
      a_phi.copyTo(phiLo, B_1);
      for(int jdir = 0; jdir < DIM; jdir++)
      {
        if(jdir == faceDir)
        {
          forallInPlace(addNormalDerivatives,B_1,
                        phiLo,slopeVLimit[jdir],cellVel,dx,dt,signLo,jdir);
          forallInPlace(addNormalDerivatives,B_1,
                        phiHi,slopeVLimit[jdir],cellVel,dx,dt,signHi,jdir);
        }
        else
        {
          forallInPlace(addTangentialDerivatives,B_1,
                        phiLo, slopeUpwind[jdir], cellVel, dx, dt, jdir);
          forallInPlace(addTangentialDerivatives,B_1,
                        phiHi, slopeUpwind[jdir], cellVel, dx, dt, jdir);
        }
      }
      //now shift phi lo and high to make left and right states on the face.
      //philo was extrapolated in the low  direction so it becomes the right state (and its index space is correct)
      //phihi was extrapolated in the high direction so it becomes the left  state (and its index space needs to be shifted up)

      BoxData<double, 1>& phiRigh = phiLo;
      BoxData<double, 1>& phiLeft = phiHi;
      phiLeft.shift(Point::Basis(faceDir));

      BoxData<double, 1>& phiHalf = a_phiHalf[faceDir];
      //now we can solve the riemann problem and get the flux all at once.
      forallInPlace(getPhiH, Bface[faceDir],
                    phiHalf, phiLeft, phiRigh, a_vel[faceDir], doingVel);
    }

    //finally we need to take the divergence of the flux.
    return;
  }
}
