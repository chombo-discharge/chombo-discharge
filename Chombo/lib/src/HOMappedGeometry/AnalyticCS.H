#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _ANALYTICCS_H_
#define _ANALYTICCS_H_

#include "NewFourthOrderCoordSys.H"
#include "RefCountedPtr.H"
#include "VectorFunction.H"
#include "TensorFunction.H"

#include "NamespaceHeader.H"

//! \class AnalyticCS
//! This class uses analytic vector/tensor functions to define a
//! coordinate mapping.
class AnalyticCS : public NewFourthOrderCoordSys
{
  public:

  //! Construct an analytic coordinate mapping.
  AnalyticCS(const RealVect& a_dX,
             RefCountedPtr<VectorFunction> a_X,
             RefCountedPtr<VectorFunction> a_Xi,
             RefCountedPtr<TensorFunction> a_J);

  //! Destructor.
  ~AnalyticCS();

  // Overridden methods

  RealVect realCoord(const RealVect& a_X) const;

  RealVect mappedCoord(const RealVect& a_x) const;

  Real dXdXi(const RealVect& a_Xi, int a_dirX, int a_dirXi) const;

  protected:

  RefCountedPtr<VectorFunction> m_X, m_Xi;
  RefCountedPtr<TensorFunction> m_J;

  // Cached Jacobian matrix and corresponding coordinate.
  mutable RealTensor m_J0;
  mutable RealVect m_Xi0;

  private:

  AnalyticCS();

};

//! \class AnalyticCSFactory
class AnalyticCSFactory: public NewCoordSysFactory
{
  public:

  /// constructor
  AnalyticCSFactory(RefCountedPtr<VectorFunction> a_X,
                    RefCountedPtr<VectorFunction> a_Xi,
                    RefCountedPtr<TensorFunction> a_J);

  ///
  virtual NewCoordSys* getCoordSys(const ProblemDomain& a_domain,
                                   const RealVect& a_dx) const;

  protected:

  RefCountedPtr<VectorFunction> m_X, m_Xi;
  RefCountedPtr<TensorFunction> m_J;
};


#include "NamespaceFooter.H"

#endif
