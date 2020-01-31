#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EXTRAPBC_H_
#define _EXTRAPBC_H_

#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "LoHiSide.H"
#include "GhostBC.H"
#include "NamespaceHeader.H"
/// Class to implement linear extrapolation boundary conditions
/**
   BoxGhostBC class implementation for linear extrapolation
   boundary conditions
   Recall that BoxGhostBC is a class to encapsulate the
   operations  of ghost-cell boundary conditions at a face.
   These functions are only meant to be called by DomainGhostBC

 */
class ExtrapBC : public BoxGhostBC
{
public:
  friend class DomainGhostBC;

  ///
  ~ExtrapBC();

  ///
  ExtrapBC();

  ///
  ExtrapBC(int dir, Side::LoHiSide sd);

  ///
  ExtrapBC(int dir, Side::LoHiSide sd,
            const Interval& a_comps);

protected:

  ///
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real a_dx,
               const ProblemDomain& a_domain) const;

  ///
  virtual void applyHomogeneousBCs(FArrayBox& a_state,
                                   const ProblemDomain& a_domain,
                                   Real a_dx) const;

  ///
  virtual void applyInhomogeneousBCs(FArrayBox& a_state,
                                     const ProblemDomain& a_domain,
                                     Real a_dx) const;

  ///
  virtual void applyExtrapBCs(FArrayBox& a_state,
                              const ProblemDomain& a_domain,
                              Real a_dx) const;


  /// deprecated interface
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real a_dx,
               const Box& a_domain) const;

  /// deprecated interface
  virtual void applyHomogeneousBCs(FArrayBox& a_state,
                                   const Box& a_domain,
                                   Real a_dx) const;

  /// deprecated interface
  virtual void applyInhomogeneousBCs(FArrayBox& a_state,
                                     const Box& a_domain,
                                     Real a_dx) const;

  /// deprecated interface
  virtual void applyExtrapBCs(FArrayBox& a_state,
                              const Box& a_domain,
                              Real a_dx) const;


  ///
  BoxGhostBC* new_boxghostbc() const;

private:

  ExtrapBC(const ExtrapBC& neuin)
  {
  }

  void operator=(const ExtrapBC& neuin)
  {
  }
};



#include "NamespaceFooter.H"
#endif
