#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _EBNORMALIZEBYVOLUMEFRACTION_H_
#define _EBNORMALIZEBYVOLUMEFRACTION_H_

#include "EBCellFAB.H"
#include "LevelData.H"
#include "EBLevelGrid.H"
#include "AggStencil.H"
#include "NamespaceHeader.H"

///
/** 
 *  This class represents an operation that takes a volume-fraction-weighted
 *  quantity kappa*Q (where kappa is the volume fraction) and computes Q
 *  from it by normalizing it w.r.t. the volume fractions of cells within a
 *  monotone path radius. The radius may either be specified explicitly or
 *  computed from a given EBLevelGrid.
 */
class EBNormalizeByVolumeFraction
{
public:

  ///
  /** Construct an EBNormalizeByVolumeFraction object associated with the
   *  associated with the given EBLevelGrid.
   *  \param EBLevelGrid The level grid for which data can be normalized
   *           by the resulting object.
   wacky interface is for backward compatibility
   if data is null, define is not called.
   */
  EBNormalizeByVolumeFraction(const EBLevelGrid          & a_eblg,
                              const LevelData<EBCellFAB> * a_data  = NULL,
                              const int                    a_radius = 1)
  {
    CH_assert(a_radius > 0);
    m_eblg      = a_eblg;
    m_radius    = a_radius;
    m_isDefined = false;
    if(a_data != NULL)
      {
        define(*a_data);
      }
  }

  
  ///define all the ste3ncils
  void define(const LevelData<EBCellFAB>& a_data);


  /** Destructor. */
  ~EBNormalizeByVolumeFraction()
    {
    }

  /** Compute the normalized quantity corresponding to the given
   *  volume-fraction-weighted quantity by volume-averaging that quantity
   *  over every monotone path radius for each cell in the level.
   *  \param a_Q A volume-fraction-weighted quantity to be normalized.
   *  \param a_compInterval An interval representing the range of components
   *                        to be normalized by this operation.
   */
  void operator()(LevelData<EBCellFAB>& a_Q,
                  const Interval& a_compInterval);

  /** Compute the normalized quantity corresponding to the given
   *  volume-fraction-weighted quantity by volume-averaging that quantity
   *  over every monotone path radius for each cell in the level. This
   *  version of the correction normalizes all components of \a a_Q.
   *  \param a_Q A volume-fraction-weighted quantity to be normalized.
   */
  void operator()(LevelData<EBCellFAB>& a_Q);
 
  void normalize(LevelData<EBCellFAB>& a_Q,
                 const Interval& a_compInterval) 
  {
    (*this)(a_Q, a_compInterval);
  }

  void normalize(LevelData<EBCellFAB>& a_Q) 
  {
    (*this)(a_Q);
  }
protected:

  void getLocalStencil(VoFStencil      & a_stencil, 
                       const VolIndex  & a_vof, 
                       const DataIndex & a_dit);

  EBNormalizeByVolumeFraction();
  EBNormalizeByVolumeFraction(const EBNormalizeByVolumeFraction&);
  EBNormalizeByVolumeFraction& operator=(const EBNormalizeByVolumeFraction&);

private:
  bool m_isDefined;
  int m_radius;
  LayoutData<RefCountedPtr<AggStencil<EBCellFAB, EBCellFAB> > > m_stencil;
  EBLevelGrid  m_eblg;
};

#include "NamespaceFooter.H"
#endif
