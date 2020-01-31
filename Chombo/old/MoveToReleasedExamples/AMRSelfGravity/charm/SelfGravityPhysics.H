#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _SELFGRAVITYPHYSICS_H_
#define _SELFGRAVITYPHYSICS_H_

#include <string>
using std::string;

#include "GodunovPhysics.H"
#include "FluxBox.H"
#include "LGintegrator.H"

/// GodonovPhysics-derived class encapsulating the self-gravity problem
class SelfGravityPhysics : public GodunovPhysics
{
public:
  /// Constructor
  SelfGravityPhysics();

  /// Destructor
  virtual ~SelfGravityPhysics();

  /// Compute the maximum wave speed
  virtual Real getMaxWaveSpeed(const FArrayBox& a_U,
                               const Box&       a_box);

  /// Compute the maximum wave speed with gravity term
  virtual Real getMaxWaveSpeedWithSource(const FArrayBox& a_U,
                                         const FArrayBox& a_S,
                                         const Real&      a_dx,
                                         const Box&       a_box);

  /// Object factory for this class
  virtual GodunovPhysics* new_godunovPhysics() const;

  /// Names of the conserved variables
  virtual Vector<string> stateNames();

  /// Names of the plot variables
  virtual Vector<string> plotNames();

  /// Number of conserved variables
  virtual int numConserved();

  /// Number of flux variables
  /**
     This can be greater than the number of conserved variables if
     addition fluxes/face-centered quantities are computed.
   */
  virtual int numFluxes();

  virtual bool isDefined();

  /// Number of primitive variables
  virtual int numPrimitives();

  ///  Number of primitive variables for which slopes are computed
  virtual int numSlopes();

  /// Compute the primitive variables from the conserved variables within a_box
  virtual void consToPrim(FArrayBox&       a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box);

  /// Compute the fluxes from primitive variable values on a face
  virtual void getFlux(FArrayBox&       a_flux,
                       const FArrayBox& a_whalf,
                       const int&       a_dir,
                       const Box&       a_box);

  /// Transform a_dW from primitive to characteristic variables
  /**
     On input, a_dW contains the increments of the primitive variables. On
     output, it contains the increments in the characteristic variables.

     IMPORTANT NOTE: It is assumed that the characteristic analysis puts the
     eigenvalues in order from smallest to largest and orders the eigenvalues
     and characteristic variables accordingly.
  */
  virtual void charAnalysis(FArrayBox&       a_dW,
                            const FArrayBox& a_W,
                            const int&       a_dir,
                            const Box&       a_box);

  /// Transform a_dW from characteristic to primitive variables
  virtual void charSynthesis(FArrayBox&       a_dW,
                             const FArrayBox& a_W,
                             const int&       a_dir,
                             const Box&       a_box);

  /// Compute the characteristic values (eigenvalues).
  virtual void charValues(FArrayBox&       a_lambda,
                          const FArrayBox& a_W,
                          const int&       a_dir,
                          const Box&       a_box);

  /// Add to (increment) the source terms given the current state
  virtual void incrementSource(FArrayBox&       a_localS,
                               const FArrayBox& a_S,
                               const FArrayBox& a_W,
                               const Real&      a_time,
                               const Real&      a_dt,
                               const Box&       a_box);
  ///
  virtual void incrementSource(FArrayBox&       a_S,
                               const FArrayBox& a_W,
                               const Box&       a_box)
  {
  }

  /// Compute the solution to the Riemann problem.
  virtual void riemann(FArrayBox&       a_WStar,
                       const FArrayBox& a_WLeft,
                       const FArrayBox& a_WRight,
                       const FArrayBox& a_W,
                       const Real&      a_time,
                       const int&       a_dir,
                       const Box&       a_box);

  /// Post-normal predictor calculation.
  virtual void postNormalPred(FArrayBox&       a_dWMinus,
                              FArrayBox&       a_dWPlus,
                              const FArrayBox& a_W,
                              const Real&      a_dt,
                              const Real&      a_dx,
                              const int&       a_dir,
                              const Box&       a_box);

  /// Compute the quasilinear update A*dW/dx.
  virtual void quasilinearUpdate(FArrayBox&       a_AdWdx,
                                 const FArrayBox& a_wHalf,
                                 const FArrayBox& a_W,
                                 const Real&      a_scale,
                                 const int&       a_dir,
                                 const Box&       a_box);

  /// synchronize entropy and pressure
  virtual void synchronize(FArrayBox&       a_U,
                           const FArrayBox& a_UOld,
                           const Box&       a_box);

  ///
  virtual void setupSourceTerm(FArrayBox&       a_S,
                               const FArrayBox& a_force,
                               const FArrayBox& a_U,
                               const Real&      a_time,
                               const Real&      a_dt,
                               const Box&       a_box);

  /// apply source and returns max dt based on CFL condition
  virtual Real applySource(FArrayBox&       a_U,
                           const FArrayBox& a_UOld,
                           const FArrayBox& a_force,
                           const Real&      a_time,
                           const Real&      a_dt,
                           const Box&       a_box);

  ///
  virtual void forceCorrection(FArrayBox&       a_U,
                               const FArrayBox& a_force,
                               const Real&      a_time,
                               const Real&      a_dt,
                               const Box&       a_box);

  ///
  virtual void setPressureToEntropy(FArrayBox& a_U,
                                    const Box& a_box);

  /// Interval within the primitive variables corresponding to the velocities
  virtual Interval velocityInterval();

  /// Component index within the primitive variables of the pressure
  virtual int pressureIndex();

  /// Component index within the primitive variables of the pressure
  virtual int densityIndex();

  /// Component index within the primitive variables of the bulk modulus
  virtual int bulkModulusIndex();

  /// Component index within the primitive variables of the entropy
  virtual int entropyIndex();

  /// Used to limit the absolute value of a "pressure" difference
  virtual Real smallPressure();

  /// set smallPressure parameter
  virtual void setSmallPressure(const Real& a_smallPressure);

protected:
  /// Used to limit the absolute value of a "pressure" difference
  Real m_smallPressure;

private:
  /// Disallowed for all the usual reasons
  void operator=(const SelfGravityPhysics& a_input)
  {
    MayDay::Error("invalid operator");
  }

  /// Disallowed for all the usual reasons
  SelfGravityPhysics(const SelfGravityPhysics& a_input)
  {
    MayDay::Error("invalid operator");
  }
};


#endif
