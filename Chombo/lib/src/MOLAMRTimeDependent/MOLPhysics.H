#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MOLPHYSICS_H_
#define _MOLPHYSICS_H_

#include <string>
using std::string;

#include "FArrayBox.H"
#include "Vector.H"
#include "PhysIBC.H"
#include "MOLUtilities.H"

#include "NamespaceHeader.H"

class HDF5HeaderData;

///
/**
   The base class MOLPhysics provides the physics-dependent components
   for a higher-order method for a single patch: Riemann solver and
   transformations between conserved, primitive, and flux variables.
   This class is essentially pure virtual, i.e. all of its member functions
   are virtual, and the ones that have default implementations are ones
   which are optionally defined, i.e. the default definition is to send
   an error message. Physics-dependent versions of this class that are
   required in real applications are derived from this class by inheritance.
*/
class MOLPhysics
{
public:
  /// Constructor
  /**
   */
  MOLPhysics();

  /// Destructor
  /**
   */
  virtual ~MOLPhysics();

  /// Get the initial and boundary condition object
  /**
   */
  PhysIBC* getPhysIBC() const;

  /// Set the initial and boundary condition object
  /**
   */
  void setPhysIBC(PhysIBC* a_bc);

  /// Define the object
  /**
   */
  virtual void define(const ProblemDomain& a_domain,
                      const Real&          a_dx);

  /// Set the current box (default implementation - do nothing)
  /**
   */
  virtual void setCurrentBox(const Box& a_currentBox);

  /// Returns true if 4th-order artificial viscosity is defined.
  /**
   */
  virtual bool fourthOrderArtificialViscosityIsDefined() const;

  /// Defines fourth-order artifical viscosity strong shock threshold.
  /**
   */
  virtual void setFourthOrderArtificialViscosityParameter(const Real& M0sq);
 /// Returns fourth-order artifical viscosity strong shock threshold.
  /**
   */
  virtual Real getFourthOrderArtificialViscosityParameter() const;

  /// Compute the maximum wave speed.
  /**
   */
  virtual Real getMaxWaveSpeed(const FArrayBox& a_U,
                               const Box&       a_box) = 0;

  /// Compute the maximum wave speed
  /**
   */
  virtual void soundSpeed(FArrayBox& a_speed,
                          const FArrayBox& a_U,
                          const Box&       a_box);

  /// Object factory for this class
  /**
   */
  virtual MOLPhysics* new_molPhysics() const = 0;

  /// Copy our values into one of us
  virtual void copyTo(MOLPhysics* a_mol) const;

  /// Number of conserved variables
  /**
     Return the number of conserved variables.
  */
  virtual int numConserved() const = 0;

  /// Names of the conserved variables
  /**
     Return the names of the conserved variables.
  */
  virtual Vector<string> stateNames() = 0;

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
  */
  virtual int numFluxes() const = 0;

  /// Compute the fluxes from primitive variable values on a face
  /**
     This has a default implementation which throws an error.  The method is
     here so that the default implementation of "computeUpdate" can use it
     and the user can supply it.  It has an implementation so if the user
     redefines "computeUpdate" they aren't force to implement "getFlux" -
     which is only used by the default implementation of "computeUpdate".
   */
  virtual void getFlux(FArrayBox&       a_flux,
                       const FArrayBox& a_WHalf,
                       const int&       a_dir,
                       const Box&       a_box);

  /// Is the object completely defined
  /**
     Return true if the object is completely defined.
  */
  virtual bool isDefined() const;

  /// Number of primitive variables
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
  */
  virtual int numPrimitives() const = 0;

  /// Compute the solution to the Riemann problem.
  /**
     Given input left and right states in a direction, a_dir, compute a
     Riemann problem and generate fluxes at the faces within a_box.
  */
  virtual void riemann(/// face-centered solution to Riemann problem
                       FArrayBox&       a_WStar,
                       /// left state, on cells to left of each face
                       const FArrayBox& a_WLeft,
                       /// right state, on cells to right of each face
                       const FArrayBox& a_WRight,
                       /// state on cells, used to set boundary conditions
                       const FArrayBox& a_W,
                       /// current time
                       const Real&      a_time,
                       /// direction of faces
                       const int&       a_dir,
                       /// face-centered box on which to set a_WStar
                       const Box&       a_box) = 0;

  /// Compute primitive variables from conserved variables.
  /**
   */
  virtual void consToPrim(FArrayBox&       a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box) = 0;

  /// Interval within the primitive variables corresponding to the velocities
  /**
     Return the interval of component indices within the primitive variable
     of the velocities.  Used for slope flattening (slope computation) and
     computing the divergence of the velocity (artificial viscosity).
   */
  virtual Interval velocityInterval() const = 0;

  /// Interval within the flux variables corresponding to vector flux
  virtual Interval vectorFluxInterval() const;

  /// Component index within the primitive variables of the pressure
  /**
     Return the component index withn the primitive variables for the
     pressure.  Used for slope flattening (slope computation).
   */
  virtual int pressureIndex() = 0;

  /// Used to limit the absolute value of a "pressure" difference
  /**
     Return a value that is used by slope flattening to limit (away from
     zero) the absolute value of a slope in the pressureIndex() component
     (slope computation).
   */
  virtual Real smallPressure() = 0;

  /// Component index within the primitive variables of the bulk modulus
  /**
     Return the component index withn the primitive variables for the
     bulk modulus.  Used for slope flattening (slope computation) used
     as a normalization to measure shock strength.
   */
  virtual int bulkModulusIndex() = 0;

 /// Component index within the primitive variables of the density.
  /**
     Return the component index within the primitive variables for the
     density.  Used for fourth-order accurate artificial viscosity.
   */
  virtual int densityIndex();

  /// Compute the artificial viscosity contribution to the flux
  /**
     Compute the artificial viscosity contribution to the flux.  This has
     a default implementation but this can be overridded as needed.

     a_box is the valid cell-centered box.

     Increment a_F
     on all a_dir-faces of a_box
     using a_U
     on all cells of bx1inDomain = grow(a_box, 1) & m_domain.
   */
  virtual void artVisc(FArrayBox&       a_F,
                       const FArrayBox& a_U,
                       const Real&      a_artificialViscosity,
                       const Real&      a_currentTime,
                       const int&       a_dir,
                       const Box&       a_box);

#ifdef CH_USE_HDF5
  virtual void expressions(HDF5HeaderData& a_holder) const
  {
  }
#endif


protected:

  // Has this object been defined
  bool m_isDefined;

  // Problem domain and grid spacing
  ProblemDomain m_domain;
  Real          m_dx;

  // Object containing various methods for the MOL computation
  MOLUtilities m_util;

  // Flag to use fourth-order accurate artificial viscosity, strong shock
  // threshhold M0sq (looks like the shock Mach number squared).
  bool m_useFourthOrderArtificialViscosity;
  Real m_M0sq;

  // Initial and boundary condition object and has it been set
  PhysIBC* m_bc;
  bool     m_isBCSet;

private:
  // Disallowed for all the usual reasons
  MOLPhysics& operator=(const MOLPhysics& a_input);
  MOLPhysics(const MOLPhysics& a_input);
};

#include "NamespaceFooter.H"
#endif
