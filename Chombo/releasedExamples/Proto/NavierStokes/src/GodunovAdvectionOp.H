#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _GODUNOVADVECTIONOP_H_
#define _GODUNOVADVECTIONOP_H_

#include "Proto.H"

namespace Proto
{
  typedef Var<double,   1> Scalar;
  typedef Var<double,DIM>  Pvector;

///
/**
 */
  class GodunovAdvectionOp
  {
  public:
    /// Default constructor.
    GodunovAdvectionOp(double m_dx);



    ///  Destructor.
    ~GodunovAdvectionOp(){};

    /// computes u*phi at faces at n+1/2
    void   getFluxNPH(BoxData<double,1>            a_flux[DIM],
                      const BoxData<double,   1>   a_vel[DIM],
                      const BoxData<double,   1> & a_phi,
                      const BoxData<double,   1> & a_src,
                      const BoxData<double, DIM> & a_velCell,
                      const Proto::Box                  & a_dbx0,
                      const int                  & a_doingVel,
                      const double               & a_dt) const;


    /// imitates the patch integrator interface (just gets phi n+1/2, i+1/2
    void   advectToFaces(BoxData<double,1>            a_vface[DIM],
                         const BoxData<double,   1>   a_macVel[DIM],
                         const BoxData<double,   1> & a_phi,
                         const BoxData<double,   1> & a_src,
                         const BoxData<double, DIM> & a_velCell,
                         const Proto::Box                  & a_dbx0,
                         const int                  & a_doingVel,
                         const double               & a_dt) const;

    /**
       computes the divergence of the flux at time n+1/2 
       Given velocity and phi at time n.
    */
    void   divFluxNPH(BoxData<double,1>          & a_div,
                      const BoxData<double,   1>   a_vel[DIM],
                      const BoxData<double,   1> & a_phi,
                      const BoxData<double,   1> & a_src,
                      const BoxData<double, DIM> & a_velCell,
                      const Proto::Box                  & a_dbx0,
                      const int                  & a_doingVel,
                      const double               & a_dt) const;


    /**
       calculates the flux (= u*phi) and takes the divergence
    */
    void   divUPhi(BoxData<double,1>          & a_div,
                   const BoxData<double,   1>   a_vel[DIM],
                   const BoxData<double,   1>   a_phi[DIM],
                   const Box                  & a_region) const;

  private:
    double m_dx;

    /// Default constructor disallowed.
    GodunovAdvectionOp();
    Stencil<double> m_slopeLo[DIM];
    Stencil<double> m_slopeHi[DIM];
    Stencil<double> m_faceToCell[DIM];
    Stencil<double> m_fluxDivIncr[DIM];

  };
}
#endif //end include guard
