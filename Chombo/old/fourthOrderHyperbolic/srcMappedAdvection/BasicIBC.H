#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

//  BasicIBC.H
// ============
//
// Virtual base class through which a user specifies the initial and boundary
// conditions for a hyperbolic system of PDEs.
//

#ifndef _BASICIBC_H_
#define _BASICIBC_H_

#include "FArrayBox.H"
#include "FluxBox.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "UsingNamespace.H"
#include "CoordSys.H"
//#include "SecondOrderCoordSys.H"
#include "FourthOrderCoordSys.H"
#include "MappedPhysIBC.H"

class MappedPhysIBC;

/// physical/domain initial and boundary conditions
/**
   Virtual base class through which a user specifies the initial and
   boundary conditions for a hyperbolic system of PDEs.
 */
class BasicIBC
{
  friend class MappedPhysIBC;

public:
  /// Constructor
  /**
   */
  BasicIBC()
  {
    m_exactOnPhysicalDomain=true;
  }

  /// Destructor
  /**
   */
  virtual ~BasicIBC()
  {
  }

  /// Factory method - this object is its own factory
  /**
     Return a point to a new BasicIBC object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  virtual BasicIBC* new_basicIBC() = 0;

  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_phi,
                          const ProblemDomain& a_domain,
                          const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                          Real a_dx,
                          Real a_time=0.0) = 0;


  virtual void advVel(LevelData<FArrayBox>& a_advVel,
                      const ProblemDomain& a_domain,
                      const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                      Real a_dx,
                      Real a_time=0.0) = 0;


  /// fill ghost cell values at domain boundaries
  virtual void ghostCellBC(LevelData<FArrayBox>& a_phi,
                           const ProblemDomain& a_domain,
                           const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                           Real a_dx,
                           Real a_time) = 0;

  /// compute transverse flux boundary conditions in ghost-cell faces
  /** In the mapped-grid case, where the transverse derivative dN/dx_j
      doesn't vanish at the boundary, we need values for the fluxes on
      ghost-cell faces transverse to the boundary, in order to be able
      to compute dFlux/dx_j to dot with the metric-term gradients.
      In other words, in a (0:N-1,0:N-1) domain, fluxes on x-faces will
      need to be set for all faces where j=-1 and N.
      I think we can generally use second-order extrapolation, which
      produces O(h) gradients.  (dropping one order of accuracy at the
      boundary)
  */
  virtual void fluxBC(LevelData<FluxBox>& a_flux,
                      const ProblemDomain& a_domain,
                      const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                      Real a_dx,
                      Real a_time);

  /// if possible, compute exact solution
  /**
   */
  virtual void exactSoln(LevelData<FArrayBox>& a_phi,
                         const ProblemDomain& a_domain,
                         const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                         Real a_dx,
                         Real a_time=0.0) = 0;

   /**
      If true, exactSoln returns the "exact" cell average on the physical
      domain; otherwise, it returns the "exact" cell average of phi*J on
      the computational domain
   */
  virtual void setExactOnPhysicalDomain( const bool a_value=false )
      { m_exactOnPhysicalDomain=a_value; }


protected:
  // define() has been called
  bool          m_isDefined;
  bool          m_exactOnPhysicalDomain;

  /// solution at a given time
  /**
   */
  virtual void pointVal(FArrayBox& a_phi,
                        const ProblemDomain& a_domain,
                        const CoordSys<FArrayBox,FluxBox>& a_coordSys,
                        const Box& a_box,
                        bool a_includeJ,
                        Real a_dx,
                        Real time) = 0;


private:
  // Disallowed for all the usual reasons
  void operator=(const BasicIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  BasicIBC(const BasicIBC& a_input)
  {
    MayDay::Error("invalid operator");
  }

};

#endif
