#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EXTRAP0THORDERBC_H_
#define _EXTRAP0THORDERBC_H_

#include "GhostBC.H"
#include "NamespaceHeader.H"

/// class to do 0th-order extrapolation BC's
class Extrap0thOrderBC : public BoxGhostBC
{
public:
  ///
  friend class DomainGhostBC;

  ///
  Extrap0thOrderBC();

  ///
  ~Extrap0thOrderBC();

  ///
  Extrap0thOrderBC(int dir, Side::LoHiSide sd,
             const Interval& a_comps);

  ///
  Extrap0thOrderBC(int dir, Side::LoHiSide sd);


protected:

  ///
  virtual void applyHomogeneousBCs(FArrayBox& a_state,
                                   const ProblemDomain& a_domain,
                                   Real a_dx) const;

  ///
  virtual void applyInhomogeneousBCs(FArrayBox& a_state,
                                     const ProblemDomain& a_domain,
                                     Real a_dx) const;

  ///
  virtual void applyExtrap0thOrderBCs(FArrayBox& a_state,
                              const ProblemDomain& a_domain,
                              Real a_dx) const;


  /// this doesn't do anything since we don't use the factors
  virtual void fillBCValues(FArrayBox& a_neumfac,
                            FArrayBox& a_dircfac,
                            FArrayBox& a_inhmval,
                            Real a_dx,
                            const ProblemDomain& a_domain) const;


  /// deprecated interface
  virtual void applyHomogeneousBCs(FArrayBox& a_state,
                                   const Box& a_domain,
                                   Real a_dx) const;

  /// deprecated interface
  virtual void applyInhomogeneousBCs(FArrayBox& a_state,
                                     const Box& a_domain,
                                     Real a_dx) const;

  /// deprecated interface
  virtual void applyExtrap0thOrderBCs(FArrayBox& a_state,
                              const Box& a_domain,
                              Real a_dx) const;


  /// this doesn't do anything since we don't use the factors
  virtual void fillBCValues(FArrayBox& a_neumfac,
                            FArrayBox& a_dircfac,
                            FArrayBox& a_inhmval,
                            Real a_dx,
                            const Box& a_domain) const;

  ///
  BoxGhostBC* new_boxghostbc() const;

private:

  Extrap0thOrderBC(const Extrap0thOrderBC& extrapin)
  {
  }

  void operator=(const Extrap0thOrderBC& extrapin)
  {
  }
};

#include "NamespaceFooter.H"
#endif


