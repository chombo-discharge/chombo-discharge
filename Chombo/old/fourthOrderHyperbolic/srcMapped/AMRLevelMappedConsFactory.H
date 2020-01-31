#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELMAPPEDCONSFACTORY_H_
#define _AMRLEVELMAPPEDCONSFACTORY_H_

#include "AMRLevelFactory.H"
// #include "PhysIBC.H"
#include "AMRLevelMappedCons.H"
#include "UsingNamespace.H"

/// AMR Cons Equation factory
/**
 */
class AMRLevelMappedConsFactory : public AMRLevelFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelMappedConsFactory();

  /// Virtual constructor
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  /**
   */
  virtual ~AMRLevelMappedConsFactory();

  /// CFL number
  /**
   */
  virtual void CFL(Real a_cfl);

  /// spatial order of accuracy (default is 4)
  virtual void spaceOrder(int a_spaceOrder);

  /// sets whether to limit face values in advection (default is false)
  virtual void limitFaceValues(bool a_limitFaceValues);

  /// sets whether initial data is average; default false
  virtual void initialAverage(bool a_initialAverage);

  /// sets whether to flatten extrapolations to faces
  virtual void useFlattening(bool a_useFlattening);

  /// sets whether to turn PPM off
  virtual void noPPM(bool a_noPPM);

  /// sets m_doDeconvolution
  virtual void doDeconvolution(bool a_doDeconvolution);

  /// sets m_doFaceDeconvolution
  virtual void doFaceDeconvolution(bool a_doFaceDeconvolution);

  /// sets whether to use artificial viscosity (from divergence)
  virtual void useArtificialViscosity(bool a_useArtificialViscosity);

  /// sets coefficient of artificial viscosity (from divergence)
  virtual void artificialViscosity(Real a_artificialViscosity);

  /// sets whether to use artificial viscosity
  virtual void useArtVisc(bool a_useArtVisc);

  /// sets coefficient of artificial viscosity
  virtual void ratioArtVisc(Real a_ratioArtVisc);

  /// sets whether forward Euler is used (instead of RK4)
  virtual void forwardEuler(bool a_forwardEuler);

  /// sets whether to enforce a min value in advection, along with valeu
  virtual void enforceMinVal(bool a_enforceMinVal, Real a_minVal);

  /// Physical dimension of the longest side of the domain
  /**
   */
  virtual void domainLength(Real a_domainLength);

  void x0(const Real& a_x0);
  //  void IBC(PhysIBC* a_w0)
  //  {
  //    m_cons_ibc = dynamic_cast<PhysIBC*>(a_w0->new_physIBC());
  //  }

  ///
  /**
   */
  virtual void verbosity(const int& verbosity);

  /// Refinement threshold
  /**
   */
  virtual void refinementThreshold(Real a_refineThresh);

  /// Tag buffer size
  /**
   */
  void tagBufferSize(int a_tagBufferSize);

  /// GodunovPhysics class
  /**
   */
  void godunovPhysics(const GodunovPhysics* const a_gdnvPhysics);

  /// Initial dt multiplier
  /**
   */
  void initialDtMultiplier(Real a_initialDtMultiplier);

  /// Whether to compute timestep from cell-centered velocities (old way)
  /**
   */
  void dtFromCells(bool a_dtFromCells);

  /// Coordinate system factory
  /**
   */
  void coordinateSystemFactory(CoordSysFactory<FArrayBox, FluxBox>* a_coordSysFact);

  /**
   */
  void plotPrefix(const std::string& a_plotfile_prefix);

  /// Check that everything is defined
  /**
   */
  bool isDefined() const;

protected:
  // Some default values
  void setDefaultValues();

  int m_verbosity;

  // CFL number
  Real m_cfl;
  bool m_cflSet;
  //  PhysIBC* m_cons_ibc;

  /// spatial order of accuracy (default is 4)
  int m_spaceOrder;

  /// if true, limit face values in advection
  bool m_limitFaceValues;

  // if true, initial values are average
  bool m_initialAverage;

  // if true, use flattening
  bool m_useFlattening;

  // Avoid PPM?  Default false.
  bool m_noPPM;

  // Do deconvolution?  Default yes.
  bool m_doDeconvolution;

  // Do face deconvolution?  Default yes.
  bool m_doFaceDeconvolution;

  // use artificial viscosity (from divergence)?
  bool m_useArtificialViscosity;

  // coefficient of artificial viscosity (from divergence)
  Real m_artificialViscosity;

  // if true, use artificial viscosity
  bool m_useArtVisc;

  // if m_useArtVisc, artificial viscosity coefficient
  Real m_ratioArtVisc;

  // use forward Euler if true, or RK4 if false
  bool m_forwardEuler;

  /// if true, enforce min value
  bool m_enforceMinVal;

  /// min value to enforce
  Real m_minVal;

  // Physical dimension of the longest side of the domain
  Real m_domainLength;
  bool m_domainLengthSet;

  // Refinement threshold for gradient
  Real m_refineThresh;
  bool m_refineThreshSet;

  // Tag buffer size
  int  m_tagBufferSize;
  bool m_tagBufferSizeSet;

  // Initial dt multiplier
  Real m_initialDtMultiplier;
  bool m_initialDtMultiplierSet;

  // whether to compute timestep from cell data (the old way)
  bool m_dtFromCells;

  // Pointer to the object that supplies the physics needed by the integrator
  // (used as a factory here)
  const GodunovPhysics* m_gdnvPhysics;

  // Coordinate system factory class
  CoordSysFactory<FArrayBox,FluxBox>* m_coordSysFactPtr;
  bool m_coordSysFactSet;

  std::string m_plotfile_prefix;

private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelMappedConsFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelMappedConsFactory(const AMRLevelMappedConsFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
