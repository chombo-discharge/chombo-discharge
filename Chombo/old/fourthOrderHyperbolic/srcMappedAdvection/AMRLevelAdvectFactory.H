#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _AMRLEVELADVECTFACTORY_H_
#define _AMRLEVELADVECTFACTORY_H_

#include "AMRLevelFactory.H"
#include "BasicIBC.H"
#include "AMRLevelAdvect.H"
#include "CoordSys.H"
#include "UsingNamespace.H"

/// AMR Advect Equation factory
/**
 */
class AMRLevelAdvectFactory : public AMRLevelFactory
{
public:
  /// Null constructor
  /**
   */
  AMRLevelAdvectFactory();

  /// Virtual constructor
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  /**
   */
  virtual ~AMRLevelAdvectFactory();

  /// CFL number
  /**
   */
  virtual void CFL(Real a_cfl);

  /// spatial order of accuracy (default is 4)
  virtual void spaceOrder(int a_spaceOrder);

  /// sets whether to limit face values in advection (default is false)
  virtual void limitFaceValues(bool a_limitFaceValues);

  /// sets whether to enforce a min value in advection, along with valeu
  virtual void enforceMinVal(bool a_enforceMinVal, Real a_minVal);

  /// sets the low-order flux for min val enforcement [0:DCU,1:CTU]
  virtual void lowOrderFluxScheme(std::string a_lowOrderFluxScheme);

  /// sets whether to use redistribution to preserve positivity
  virtual void redistributeNegativeVal(bool a_redistributeNegativeVal,
                                       int a_maxRedistributionPasses);

  /// sets whether to use hyperviscous stabilization (default is false)
  virtual void useHyperviscosity(bool a_useHyperviscosity);

  /// sets hyperviscosity coefficient
  virtual void hyperviscosity(Real a_hyperviscosity);

  /// Physical dimension of the longest side of the domain
  /**
   */
  virtual void domainLength(Real a_domainLength);

  void x0(const Real& a_x0);

  void IBC(BasicIBC* a_w0)
  {
    m_advect_ibc = a_w0->new_basicIBC();
  }

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

  /// Initial dt multiplier
  /**
   */
  void initialDtMultiplier(Real a_initialDtMultiplier);

  /// Coordinate system factory
  /**
   */
  void coordinateSystemFactory(CoordSysFactory<FArrayBox,FluxBox>* a_coordSysFact);

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
  BasicIBC* m_advect_ibc;

  /// spatial order of accuracy (default is 4)
  int m_spaceOrder;

  /// if true, limit face values in advection
  bool m_limitFaceValues;

  /// if true, enforce min value
  bool m_enforceMinVal;

  /// min value to enforce
  Real m_minVal;

  // if enforcing minval, low order scheme to use
  std::string m_lowOrderFluxScheme;

  // if true, enforce postivity with redistribution
  bool m_redistributeNegativeVal;

  // if enforcing positivity with redistribution, number of passes over which to redistribute
  int m_maxRedistributionPasses;

  // if true, add hyperviscous fluxes to each flux for stabilization
  bool m_useHyperviscosity;
  Real m_hyperviscosity;

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

  // Coordinate system factory class
  CoordSysFactory<FArrayBox,FluxBox>* m_coordSysFactPtr;
  bool m_coordSysFactSet;

  std::string m_plotfile_prefix;

private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelAdvectFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelAdvectFactory(const AMRLevelAdvectFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
