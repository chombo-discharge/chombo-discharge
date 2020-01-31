#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MOLANALYTICLONLATADVECTIONPHYSICS_H_
#define _MOLANALYTICLONLATADVECTIONPHYSICS_H_

#include <string>
using std::string;

#include "CH_HDF5.H"

#include "MOLPhysics.H"
#include "NewFourthOrderCoordSys.H"
#include "RefCountedPtr.H"
#include "VectorFunction.H"
#include "ScalarFunction.H"

#include "NamespaceHeader.H"

//! \class MOLAnalyticLonLatAdvectionPhysics
//! This class implements the advection equation for a given
//! vector-valued function V(x, t).
class MOLAnalyticLonLatAdvectionPhysics: public MOLPhysics
{
  public:

  //! Constructs advection physics for a given (vector) velocity function.
  explicit MOLAnalyticLonLatAdvectionPhysics(RefCountedPtr<VectorFunction> a_advectionVelocity);

  //! Constructs advection physics for a given stream function.
  explicit MOLAnalyticLonLatAdvectionPhysics(RefCountedPtr<ScalarFunction> a_streamFunc);

  /// Destructor
  /**
   */
  ~MOLAnalyticLonLatAdvectionPhysics();

  // Overridden methods.

  Real getMaxWaveSpeed(const FArrayBox& a_U,
                       const Box&       a_box);

  void soundSpeed(FArrayBox& a_speed,
                  const FArrayBox& a_U,
                  const Box&       a_box);

  virtual MOLPhysics* new_molPhysics() const;

  int numConserved() const;

  Vector<string> stateNames();

  int numFluxes() const;

  void getFlux(FArrayBox&       a_flux,
               const FArrayBox& a_whalf,
               const int&       a_dir,
               const Box&       a_box);

  // return velocity in equiangular coordinates
  void getVelocity(FArrayBox&       a_velFab,
                   int              a_dir,
                   const Box&       a_box);

  int numPrimitives() const;

  void riemann(FArrayBox&       a_WStar,
               const FArrayBox& a_WLeft,
               const FArrayBox& a_WRight,
               const FArrayBox& a_W,
               const Real&      a_time,
               const int&       a_dir,
               const Box&       a_box);

  void consToPrim(FArrayBox&       a_W,
                  const FArrayBox& a_U,
                  const Box&       a_box);

  virtual Interval velocityInterval() const;

  virtual int pressureIndex();

  virtual Real smallPressure();

  virtual int bulkModulusIndex();

  void setCurrentTime(Real a_time);

  void setCurrentCoordSys(const NewCoordSys* a_coordSys);

#ifdef CH_USE_HDF5
  virtual void expressions(HDF5HeaderData& a_holder) const;
#endif

  private:

  RefCountedPtr<VectorFunction> m_velocityFunc;
  RefCountedPtr<ScalarFunction> m_streamFunc;
  NewFourthOrderCoordSys* m_coordSys;
  Real m_time;

  // Disallowed for all the usual reasons
  MOLAnalyticLonLatAdvectionPhysics& operator=(const MOLAnalyticLonLatAdvectionPhysics& a_input);
  MOLAnalyticLonLatAdvectionPhysics(const MOLAnalyticLonLatAdvectionPhysics& a_input);
};

#include "NamespaceFooter.H"

#endif
