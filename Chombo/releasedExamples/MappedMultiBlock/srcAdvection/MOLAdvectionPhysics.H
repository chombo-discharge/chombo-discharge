#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MOLADVECTIONPHYSICS_H_
#define _MOLADVECTIONPHYSICS_H_

#include <string>
using std::string;

#include "CH_HDF5.H"

#include "MOLMappedPhysics.H"

#include "NamespaceHeader.H"

class MOLAdvectionPhysics: public MOLMappedPhysics
{
public:
  /// Constructor
  /**
   */
  MOLAdvectionPhysics();

  /// Destructor
  /**
   */
  ~MOLAdvectionPhysics();

  /// Compute the maximum wave speed
  /**
   */
  Real getMaxWaveSpeed(const FArrayBox& a_U,
                       const Box&       a_box);

  /// For advection, set a_speed = 0.
  /**
   */
  void soundSpeed(FArrayBox& a_speed,
                  const FArrayBox& a_U,
                  const Box&       a_box);

  /// Object factory for this class
  /**
   */
  virtual MOLPhysics* new_molPhysics() const;

  /// Number of conserved variables
  /**
     Return the number of conserved variables.
  */
  int numConserved() const;

  /// Names of the conserved variables
  /**
     Return the names of the conserved variables.  A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numConserved()-1.
  */
  Vector<string> stateNames();

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
  */
  int numFluxes() const;

  /// Compute a flux from primitive variable values on a face
  /**
   */
  void getFlux(FArrayBox&       a_flux,
               const FArrayBox& a_whalf,
               const int&       a_dir,
               const Box&       a_box);

  /// Number of primitive variables
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
  */
  int numPrimitives() const;

  /// Compute the solution to the Riemann problem.
  /**
     Given input left and right states in a direction, a_dir, compute a
     Riemann problem and generate fluxes at the faces within a_box.
  */
  void riemann(FArrayBox&       a_WStar,
               const FArrayBox& a_WLeft,
               const FArrayBox& a_WRight,
               const FArrayBox& a_W,
               const Real&      a_time,
               const int&       a_dir,
               const Box&       a_box);

  /// Compute primitive variables from conserved variables.
  /**
   */
  void consToPrim(FArrayBox&       a_W,
                  const FArrayBox& a_U,
                  const Box&       a_box);

  /// Interval within the primitive variables corresponding to the velocities
  /**
     Return the interval of component indices within the primitive variable
     of the velocities.  Used for slope flattening (slope computation) and
     computing the divergence of the velocity (artificial viscosity).
   */
  virtual Interval velocityInterval() const;

  /// Component index within the primitive variables of the pressure
  /**
     Return the component index withn the primitive variables for the
     pressure.  Used for slope flattening (slope computation).
   */
  virtual int pressureIndex();

  /// Used to limit the absolute value of a "pressure" difference (away from zero)
  /**
     Return a value that is used by slope flattening to limit (away from
     zero) the absolute value of a slope in the pressureIndex() component
     (slope computation).
   */
  virtual Real smallPressure();

  /// Component index within the primitive variables of the bulk modulus
  /**
     Return the component index withn the primitive variables for the
     bulk modulus.  Used for slope flattening (slope computation) used
     as a normalization to measure shock strength.
   */
  virtual int bulkModulusIndex();

  virtual void setCurrentBox(const Box& a_currentBox);

  virtual void setCurrentCoordSys(const NewCoordSys* a_coordSys);

  virtual const FArrayBox& velFab();

  virtual const FluxBox& velFlux();

#ifdef CH_USE_HDF5
  virtual void expressions(HDF5HeaderData& a_holder) const;
#endif

protected:

  Box m_currentBox;

  bool m_haveAdvVelFab;

  bool m_haveAdvVelFlux;

  /// face-centered advection velocity
  mutable FArrayBox m_advVelFab;

  mutable FluxBox m_advVelFlux;
};

#include "NamespaceFooter.H"

#endif
