#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MBAGGSTENCIL_H_
#define _MBAGGSTENCIL_H_

#include "Vector.H"
#include "REAL.H"
#include "CH_Timer.H"
#include "IVSFAB.H"
#include "FArrayBox.H"
#include "MBStencil.H"
#include "NamespaceHeader.H"

/// 
/**
   Aggregated stencil for multiblock stuff.
   Takes IVSFAB<MBStencil> and turns that into
   arrays of offsets and weights so the stencil
   can be evaluated very quickly in comparison 
   with fab(iv, ivar) over and over.
 */
class MBAggStencil
{
public:
  ///
  /**
     IVS sent in separately because of the removeNoValid fiasco.
     sometimes the input stencil has extra cells.
   */
  MBAggStencil(const RefCountedPtr<IVSFAB<MBStencil> >  & a_stencil,
               const IntVectSet                         & a_region,
               const Box                                & a_srcBox,
               const Box                                & a_dstBox,
               const int                                & a_stencilDim);

  ///
  ~MBAggStencil()
  {
  }


  ///
  /**
     for applying a scalar operator to multiple operators.
     Sets the initial variable for destination to a_varDest
     Sets dataPtr for a_phi initial variable to a_varSrc
     Runs over ncomp components
   */
  void apply(FArrayBox                 & a_lph,
             const Vector<FArrayBox*>  & a_phi,
             const int                 & a_src,
             const int                 & a_dst,
             const int                 & a_nco,
             const bool                & a_incrementOnly) const;

  typedef struct
  {
    size_t m_offset;
    Real   m_weight;
    int    m_blockid;
  } stenelem_t;

  typedef Vector< stenelem_t > stencil_t;

protected:

  size_t getOffset(const IntVect& a_iv, const Box& a_grid);
  void computeSlabSizes();
  void computeOffsets(const RefCountedPtr<IVSFAB<MBStencil> >  & a_stencil,
                      const IntVectSet                         & a_region,
                      const Box                                & a_srcBox,
                      const Box                                & a_dstBox);

  //cached offsets and weights and such
  Vector<stencil_t>   m_mbstencil;
  Vector< size_t  >   m_dstAccess;

  // input data
  IntVect             m_srcBoxSize;
  IntVect             m_dstBoxSize;
  int                 m_stencilDim;
  //derived 
  bool                m_reducedDimStencil;
  size_t              m_srcSlabLinearSize;
  size_t              m_dstSlabLinearSize;
  int                 m_numSlabs;
private:
  ///weak construction is bad.
  MBAggStencil();

  ///deep copy for this object would kill performance
  MBAggStencil& operator=(const MBAggStencil& stenin);

  ///no reason to allow this one either.
  MBAggStencil(const MBAggStencil& stenin);

};

#include "NamespaceFooter.H"
#endif
