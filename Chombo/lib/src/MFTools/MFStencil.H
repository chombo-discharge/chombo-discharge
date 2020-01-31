#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _MFSTENCIL_H_
#define _MFSTENCIL_H_

#include "Stencils.H"
#include "Vector.H"
#include "VolIndex.H"
#include "REAL.H"
#include "Box.H"
#include "FArrayBox.H"
#include "EBIndexSpace.H"
#include "VolIndex.H"
#include "MFCellFAB.H"
#include "NamespaceHeader.H"

///
/**
   Simple stencil aggregation for multifluid problems.
 */
class MFStencil
{
public:
  /// structure to hold multifluid stencils
  struct
  {
    ///fluid where the answer will go
    int       destFluid;
    ////vof where answer will go
    VolIndex destVoF;
    ///one stencil for each fluid
    VoFStencil stenFluid[2];
  } typedef agg_t;

  ///
  /**
     Destructor
  */
  ~MFStencil()
  {;}



  ///
  /**
      pointSten is the stencil at each irregular cell
      gridBox is the box from the disjointboxlayout
      vectEBISBox[0] is the ebisbox from fluid 0
      vectEBISBox[1] is the ebisbox from fluid 1
      ghostVectLph is the ghost cells of the destination fab
      ghostVectPhi is the ghost cells of the input (phi) fab
      vardest is the variable location where data will be put
      lphi[ifluid](vof, vardest)  = blah blah
  */
  MFStencil(const Vector<agg_t>&    a_pointSten,
            const Box&              a_gridBox,
            const Vector<EBISBox>&  a_vectEBISBox,
            const IntVect&          a_ghostVectLph,
            const IntVect&          a_ghostVectPhi,
                  int               a_varDest);


  ///
  /**
     Applies stencil to each component of phi using the stencil weights and offsets to compute L.
     If incrementOnly = true, a_lofphi is incremented without any set to zero.
     If false, a_lofphi is set to zero and set equal to L(phi).
  */
  void apply(MFCellFAB& a_lofphi, const MFCellFAB& a_phi, bool incrementOnly = false) const;



  ///
  /**
     Cache lphi from inputs into internal storage.
  */
  void cache(const MFCellFAB& a_lphi) const;

  ///
  /**
     uncache lphi from internal storage into the input data holder.
  */
  void uncache(MFCellFAB& a_lphi) const;


  protected:

  void computeOffsets(const Vector<agg_t>& a_stencil);

  struct
  {
    int fluidid;
    int offset;
    bool multiValued;
  } typedef destTerm_t;

  struct
  {
    int fluidid;
    int offset;
    Real weight;
    //only used when debugging
    //VolIndex vof;
  } typedef stencilTerm;

  struct
  {
    Vector<stencilTerm> single;
    Vector<stencilTerm> multi;
  } typedef ebstencil_t;


  Box m_grid;
  EBISBox m_ebisBox[2];
  Box m_lphBox;
  Box m_phiBox;
  IntVect m_ghostVectPhi;
  IntVect m_ghostVectLph;
  int m_destVar;

  Vector<ebstencil_t> m_ebstencil;
  Vector<destTerm_t>  m_destTerms;
  mutable Vector<Real> m_cache;

  //only used when debugging
  //Vector<agg_t>  m_srcVoFs;
private:
  ///
  /*
     Default constructor.
  */
  MFStencil()
  {
    MayDay::Error("invalid operator");
  }

  ///
  /*
     Operator =.
  */
  void
  operator=(const MFStencil& stenin)
  {
    MayDay::Error("invalid operator");
  }

  ///
  /*
     Copy constructor.
  */
  MFStencil(const MFStencil& stenin)
  {
    MayDay::Error("invalid operator");
  }

};


#include "NamespaceFooter.H"
#endif
