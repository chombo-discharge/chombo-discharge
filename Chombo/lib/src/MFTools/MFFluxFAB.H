#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFFLUXFAB_H_
#define _MFFLUXFAB_H_

#include "EBFluxFAB.H"
#include "BoxLayoutData.H"
#include "EBISLayout.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "MFIndexSpace.H"
#include "NamespaceHeader.H"

/// Container for face-centered fluxes for multifluid
class MFFluxFAB
{
public:
  ///
  MFFluxFAB()
  {
  }

  ///
  virtual ~MFFluxFAB();

  ///
  MFFluxFAB(const Vector<EBISBox>& a_phaseGraphs,
            const Box& a_region, const Vector<int>& a_nVar);

  ///
  void
  setVal(Vector<Real> a_value);

  ///
  void setVal(const Real& a_value);

  ///
  int nComp(int a_phase) const
  {
    return (*m_phase[a_phase]).nComp();
  }

  ///
  int numPhases() const
  {
    return m_phase.size();
  }

  ///
  const Box& box() const
  {
    return m_box;
  }

  /// return the single-phase EBFluxFAB
  EBFluxFAB& getPhase(int a_phase)
  {
    return (*m_phase[a_phase]);
  }

  /// return a const single-phase EBFluxFAB
  const EBFluxFAB& getPhase(int a_phase) const
  {
    return (*m_phase[a_phase]);
  }

  ///
  EBFluxFAB* getPhasePtr(int a_phase)
  {
    return m_phase[a_phase];
  }

  ///
  MFFluxFAB&  copy(const MFFluxFAB& a_src);

  ///
  void copy(const Box& RegionFrom,
            const Interval& destInt,
            const Box& RegionTo,
            const MFFluxFAB& source,
            const Interval& srcInt);

  static int preAllocatable()
  {
    return 1;
  }

  ///
  int size(const Box& R, const Interval& comps) const ;

  ///
  void linearOut(void* buf, const Box& R, const Interval& comps) const ;

  ///
  void linearIn(void* buf, const Box& R, const Interval& comps);

  ///invalid but necessary for leveldata to compile
  MFFluxFAB(const Box& a_region, int a_nVar)
  {
    MayDay::Error("invalid constructor called for baseebfluxfab");
  }


private:

  void operator=(const MFFluxFAB& rhs); //unimplemented

  Vector<EBFluxFAB*> m_phase;
  Box                m_box;

};

/// Factory class to produce MFFluxFABs.
/**
   Factory class to produce MFFluxFABs.
   This is needed for LevelData construction.
 */
class MFFluxFactory
  : public DataFactory<MFFluxFAB>
{
public:

  MFFluxFactory(Vector<EBISLayout>& a_ebis, const Vector<int>& a_ncomp);

  MFFluxFactory(const MFIndexSpace& a_mf, const DisjointBoxLayout& a_dbl,
                const Box& a_domain, const Vector<int>& a_ncomps,
                int ghost);

  void define(Vector<EBISLayout>& a_ebis, const Vector<int>& a_ncomp);

  virtual ~MFFluxFactory();
  /// factory function.
  /**
      Creates a new mffluxfab object
      and returns a pointer to it.  Responsibility
      for calling operator 'delete' on this pointer is passed to the user.
  */
  virtual MFFluxFAB* create(const Box& a_box, int a_ncompsIgnored,
                            const DataIndex& a_dit) const;

  const EBISLayout& getEBISLayout(int phase) const
  {
    return m_ebis[phase];
  }

private:

  Vector<EBISLayout> m_ebis;
  Vector<int>        m_ncomp;

};

#include "NamespaceFooter.H"
#endif
