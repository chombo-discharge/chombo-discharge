#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFCELLFAB_H_
#define _MFCELLFAB_H_

#include "EBCellFAB.H"
#include "BoxLayoutData.H"
#include "EBISLayout.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include "MFIndexSpace.H"
#include "NamespaceHeader.H"

class MFCellFAB
{
public:
  MFCellFAB()
  {
  }

  virtual ~MFCellFAB();

  MFCellFAB(const Vector<EBISBox>& a_phaseGraphs,
            const Box& a_region, const Vector<int>& a_nVar);

  void
  setVal(Vector<Real> a_value);

  void setVal(const Real& a_value);

  void setVal(Real a_val, const Box& a_region, int startcomp,int ncomp);

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

  EBCellFAB& getPhase(int a_phase)
  {
    return (*m_phase[a_phase]);
  }
  const EBCellFAB& getPhase(int a_phase) const
  {
    return (*m_phase[a_phase]);
  }
  EBCellFAB* getPhasePtr(int a_phase)
  {
    return m_phase[a_phase];
  }

  ///
  void copy(const Box& RegionFrom,
            const Interval& destInt,
            const Box& RegionTo,
            const MFCellFAB& source,
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
  MFCellFAB(const Box& a_region, int a_nVar)
  {
    MayDay::Error("invalid constructor called for baseebcellfab");
  }
  ///
  /**
     Negate the values in all locations
   */
  MFCellFAB& negate(void)
  {
    this->mult(-1);
    return *this;
  }

  ///
  /**
     Both fabs need the same ebisBox and region.
   */
  MFCellFAB& operator+=(const MFCellFAB& a_mffab);

  ///
  /**
     Both fabs need the same ebisBox and region.
  */
  MFCellFAB& plus(const MFCellFAB& a_ebfab,
                  int a_srccomp,
                  int a_destcomp,
                  int a_numcomp);


  MFCellFAB& plus(const MFCellFAB& a_src,
                  const Box&       a_srcbox,
                  const Box&       a_destbox,
                  const Real&      a_scale,
                  int              a_srccomp,
                  int              a_destcomp,
                  int              a_numcomp);
  ///
  /**
     Both fabs need the same ebisBox and region.
   */
  MFCellFAB& operator-=(const MFCellFAB& a_mffab);

  ///
  /**
     Both fabs need the same ebisBox and region.
  */
  MFCellFAB& minus(const MFCellFAB& a_ebfab,
                   int a_srccomp,
                   int a_destcomp,
                   int a_numcomp);


  ///
  /**
     Both fabs need the same ebisBox and region.
   */
  MFCellFAB& operator*=(const MFCellFAB& a_mffab);

  ///
  /**
     Both fabs need the same ebisBox and region.
  */
  MFCellFAB& mult(const MFCellFAB& a_ebfab,
                  int a_srccomp,
                  int a_destcomp,
                  int a_numcomp);


  ///
  /**
     Both fabs need the same ebisBox and region.
   */
  MFCellFAB& operator/=(const MFCellFAB& a_mffab);


  ///
  /**
     Both fabs need the same ebisBox and region.
  */
  MFCellFAB& divide(const MFCellFAB& a_ebfab,
                    int a_srccomp,
                    int a_destcomp,
                    int a_numcomp);


  ///
  /**
     Both fabs need the same ebisBox and region.
   */
  MFCellFAB& operator+=(const Real& a_scalar);

  ///
  /**
     Both fabs need the same ebisBox and region.
   */
  MFCellFAB& operator-=(const Real& a_scalar);

  ///
  /**
     Both fabs need the same ebisBox and region.
   */
  MFCellFAB& operator*=(const Real& a_scalar);

  ///
  /**
     Both fabs need the same ebisBox and region.
   */
  MFCellFAB& mult(Real a_scalar);

  ///
  /**
     Both fabs need the same ebisBox and region.
   */
  MFCellFAB& operator/=(const Real& a_scalar);

  ///
  /**
     Current FAB += a_src FAB * a_scalar.  Both fabs need the same ebisBox
     and region.
   */
  MFCellFAB& plus(const MFCellFAB& a_src,
                  Real             a_scalar);


  MFCellFAB& copy(const MFCellFAB& a_src);

  /// (Not implemented) Returns the Lp-norm of this MFCellFAB
  /**
     (Not implemented) Returns the Lp-norm of this MFCellFAB using components
     (a_comp : a_comp + a_numcomp - 1).  a_power < 0 -> ERROR.
     a_power = 0  -> infinity norm (max norm).
     a_power = 1  -> L1-norm
     a_power > 1  -> Lp-norm
  */
  virtual Real norm(int a_power,
                    int a_comp,
                    int a_numComp) const;

  /// (Not implemented) Returns the Lp-norm of this MFCellFAB within a region
  /**
     (Not implemented) Returns the Lp-norm of this MFCellFAB using components
     (a_comp : a_comp + a_numcomp - 1) and within the a_subbox.  a_power < 0
     -> ERROR.
     a_power = 0 -> infinity norm (max norm).
     a_power = 1 -> L1-norm
     a_power > 1 -> Lp-norm
  */
  virtual Real norm(const Box& a_subbox,
                    int        a_power,
                    int        a_comp,
                    int        a_numComp) const;

  /// (Not implemented) Returns a sum of powers of a subset of this MFCellFAB
  /**
     (Not implemented) Returns a sum of powers of a subset of this MFCellFAB,
     specifically components a_comp to a_comp+a_numcomp-1 within a_subbox.
     a_power >= 2 only.

  */
  virtual Real sumPow(const Box& a_subbox,
                      int        a_power,
                      int        a_comp,
                      int        a_numComp) const;

  /// (Not implemented) Return the dot product of this MFCellFAB with another
  /**
     (Not implemented) Return the dot product of this MFCellFAB and "ebfab2"
     over their overlap region and all components.
  */
  Real dotProduct(const MFCellFAB& ebfab2) const;

  Real dotProduct(const MFCellFAB& ebfab2, const Box& a_box) const;

private:

  void operator=(const MFCellFAB& rhs); //unimplemented

  Vector<EBCellFAB*> m_phase;
  Box                m_box;

};

/// Factory class to produce MFCellFABs.
/**
   Factory class to produce MFCellFABs.
   This is needed for LevelData construction.
 */
class MFCellFactory
  : public DataFactory<MFCellFAB>
{
public:

  MFCellFactory(Vector<EBISLayout>& a_ebis, const Vector<int>& a_ncomp);

  MFCellFactory(const MFIndexSpace& a_mf, const DisjointBoxLayout& a_dbl,
                const Box& a_domain, const Vector<int>& a_ncomps,
                int ghost);

  void define(Vector<EBISLayout>& a_ebis, const Vector<int>& a_ncomp);

  virtual ~MFCellFactory();
  /// factory function.
  /**
      Creates a new mfcellfab object
      and returns a pointer to it.  Responsibility
      for calling operator 'delete' on this pointer is passed to the user.
  */
  virtual MFCellFAB* create(const Box& a_box, int a_ncompsIgnored,
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
