#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#include <cmath>

#include "LoHiSide.H"
#include "BoxIterator.H"
#include "FourthOrderUtil.H"
#include "PolynomialPatchIBC.H"

// Null constructor
PolynomialPatchIBC::PolynomialPatchIBC()
{
    m_params_are_set = false;
}
// Sets parameters used by initial conditions
//     a_r0            - Full width of PolynomialPatch.
void PolynomialPatchIBC::setParams(Real a_rad,
                                   RealVect a_center,
                                   Real a_mag)
{
    CH_assert(m_params_are_set == false);

    m_rad = a_rad;
    m_center = a_center;
    m_mag = a_mag;

    // coefficients follow from mag
    m_A = m_mag;
    m_B = 0.0;
    m_C = -10*m_mag;
    m_D = 20*m_mag;
    m_E = -15.0*m_mag;
    m_F = 4*m_mag;

    m_params_are_set = true;
}


/// set uniform velocity field
void
PolynomialPatchIBC::setUniformVel(const RealVect& a_vel)
{
  m_velType = uniform;
  m_uniformVel = a_vel;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new BasicIBC object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
BasicIBC* PolynomialPatchIBC::new_basicIBC()
{
  PolynomialPatchIBC* retval = new PolynomialPatchIBC();
  retval -> setParams(m_rad,m_center, m_mag);
  //[NOTE: do this even though setParams() will set it true
  //       because it might not have been set in the source
  //       object (this).  Of course, that would be a bad idea
  //       because then any object created by the factor wont
  //       be usable.  Caveat usor.  -dbs]

  retval ->m_velType = m_velType;
  retval->m_uniformVel = m_uniformVel;
  retval->m_params_are_set = m_params_are_set ;

  return static_cast<BasicIBC*>(retval);
}

// Set up initial conditions
void PolynomialPatchIBC::initialize(LevelData<FArrayBox>& a_phi,
                                    const ProblemDomain& a_domain,
                                    const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                                    Real a_dx,
                                    Real a_time)
{
  CH_assert(m_params_are_set == true);

  DataIterator dit = a_phi.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& phi = a_phi[dit()];


    // Box of current grid
    Box phiBox = phi.box();

    pointVal(phi, a_domain, a_coordSys,
             phiBox, true, a_dx, a_time);
  }

  // convert point values into 4th-order cell averages
  fourthOrderAverage(a_phi);
}


void
PolynomialPatchIBC::advVel(LevelData<FArrayBox>& a_advVel,
                           const ProblemDomain& a_domain,
                           const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                           Real a_dt, Real a_time)
{
   DataIterator dit = a_advVel.dataIterator();
   for (dit.begin(); dit.ok(); ++dit)
   {
      for (int dir=0; dir<SpaceDim; dir++)
      {
         if (m_velType == uniform)
         {
            a_advVel[dit].setVal(m_uniformVel[dir],dir);
         }
         else
         {
            MayDay::Error("PolynomialPatchBC::advVel -- bad velType");
         }
      }
   }
  fourthOrderAverage(a_advVel);
}


/// fill ghost cell values at domain boundaries
void
PolynomialPatchIBC::ghostCellBC(LevelData<FArrayBox>& a_phi,
                                const ProblemDomain& a_domain,
                                const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                                Real a_dx,
                                Real a_time)
{
  MayDay::Error("PolynomialPatchIBC::getGhostBC not implemented yet");
}



/// compute exact solution
void
PolynomialPatchIBC::exactSoln(LevelData<FArrayBox>& a_phi,
                              const ProblemDomain& a_domain,
                              const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                              Real a_dx,
                              Real a_time)
{
//  initialize(a_phi, a_domain, a_dx, a_time);
}


/// solution at a given time
/**
 */
void
PolynomialPatchIBC::pointVal(FArrayBox& a_phi,
                             const ProblemDomain& a_domain,
                             const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                             const Box& a_box,
                             bool a_includeJ,
                             Real a_dx,
                             Real a_time)
{

  // compute initial values

  // compute current center
  RealVect center = m_center;
  if (m_velType == uniform)
    {
      center += a_time*m_uniformVel;
    }

  // 0.5 for cell-centered, 0 for node-centering
  RealVect meshOffset(0.5*RealVect::Unit);
  for (int dir=0; dir<SpaceDim; dir++)
    {
      if (a_box.type(dir) == IndexType::NODE) meshOffset[dir] = 0.0;
    }


  BoxIterator bit(a_box);
  for (bit.begin();bit.ok();++bit)
    {
      IntVect iv = bit();
      Real dist = 0.;
      RealVect Xi;
      for (int idir = 0;idir < SpaceDim; idir++)
        {
          Xi[idir] = ((meshOffset[idir]+iv[idir])*a_dx);
        }

      RealVect X(a_coordSys.realCoord(Xi));


      for (int idir=0; idir<SpaceDim; idir++)
        {
          Real loc = X[idir] - center[idir];
          if (a_domain.isPeriodic(idir)  )
            {
              loc = min(abs(loc), abs(loc + 1.0));
              loc = min(abs(loc), abs(loc - 1.0));
            }

          dist = dist + loc*loc;
        }
      dist = dist/m_rad;

      for (int icomp = 0;icomp < a_phi.nComp();icomp++)
        {
          Real phiVal = m_A
            +m_B*dist
            +m_C*dist*dist
            +m_D*dist*dist*dist
            +m_E*dist*dist*dist*dist
            +m_F*dist*dist*dist*dist*dist;

          if (a_includeJ)
            {
              phiVal *= a_coordSys.pointwiseJ(X);
            }
          a_phi(iv,icomp) = phiVal;
        }
    }

}
