#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Multipoles.H
// petermc, 25 Aug 2004

#ifndef _MULTIPOLES_H_
#define _MULTIPOLES_H_

#include <cstdlib>
#include <iostream>
#include "SPACE.H"
#include <cmath>
#include "BoxLayout.H"
#include "BoxLayoutData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "Tuple.H"
#include "REAL.H"
#include "InterpOnFace.H"

using std::vector;

#include "NamespaceHeader.H"

/// Computes and evaluates multipole expansions.
class Multipoles
/** Computes and evaluates multipole expansions.
*/
{

public:

  /**
     \name Constructors, destructor and defines
  */
  /*@{*/

  ///
  /** Default constructor leaves Multipoles undefined.
  */
  Multipoles();

  ///
  /**
     Full constructor.
     Calls full define() function with identical arguments.
  */
  Multipoles(/// NODE-centered domain with faces containing charges
             const Box&        a_srcBox,
             /// NODE-centered domain on faces of which to evaluate
             const Box&        a_dstBox,
             /// mesh spacing
             const RealVect&   a_dx,
             /// length of square patches
             int               a_patchSize,
             /// degree of multipole expansion
             int               a_multipoleOrder,
             /// coarsening ratio for destination faces
             int               a_dstFaceCoarsening,
             /// width of border for polynomial interpolation
             int               a_interpBorder,
             /// run in parallel?
             bool              a_parallel,
             /// solve for coarser-level solution also?
             bool              a_getOuterCoarse = false,
             /// coarsening ratio to coarser rectangular domain
             int               a_outerCoarsening = 1,
             /// width of layer of coarse nodes around a_dstBox coarsened by a_outerCoarsening
             int               a_outerCoarseBuffer = 0);

  ///
  /** Destructor.
   */
  virtual ~Multipoles();

  ///
  /**
     Full define function.
  */
  virtual void define(/// NODE-centered domain with faces containing charges
                      const Box&        a_srcBox,
                      /// NODE-centered domain on faces of which to evaluate
                      const Box&        a_dstBox,
                      /// mesh spacing
                      const RealVect&   a_dx,
                      /// length of square patches
                      int               a_patchSize,
                      /// degree of multipole expansion
                      int               a_multipoleOrder,
                      /// coarsening ratio for destination faces
                      int               a_dstFaceCoarsening,
                      /// width of border for polynomial interpolation
                      int               a_interpBorder,
                      /// run in parallel?
                      bool              a_parallel,
                      /// solve for coarser-level solution also?
                      bool              a_getOuterCoarse = false,
                      /// coarsening ratio to coarser rectangular domain
                      int               a_outerCoarsening = 1,
                      /// width of layer of coarse nodes around a_dstBox coarsened by a_outerCoarsening
                      int               a_outerCoarseBuffer = 0);

  /*@}*/

  /**
     \name Access functions
  */
  /*@{*/

  ///
  /** Returns <tt>true</tt> if this object was created with the defining
      constructor or if define() has been called.
  */
  bool isDefined() const;

  /*@}*/

  /**
     \name Parameter-setting functions
  */
  /*@{*/

  ///
  /** Set verbosity.  Default is 0.
   */
  void setVerbose(int a_verbose);

  /*@}*/

  /**
     \name Data modification functions
  */
  /*@{*/

  /** Delete internal memory and set pointers to NULL.
   */
  void clearMemory();

  /** Set default values.  Does not deal with memory.
   */
  void setDefaultValues();

  /*@}*/

  /**
     \name Functions that change the multipole coefficients
  */
  /*@{*/

  ///
  /** Compute coefficients of multipole expansions due to charges
      on all source faces.
  */
  void setCoeffs(const Vector<FArrayBox*>&   a_faceCharge);

  ///
  /** Compute coefficients of multipole expansions due to charges
      on one source face.
  */
  void setCoeffs(const FArrayBox&        a_faceCharge,
                 int                     a_idirSrc,
                 const Side::LoHiSide&   a_side);

  ///
  /** Clear multipole expansions for all source faces.
  */
  void clearCoeffs();

  /*@}*/

  /**
     \name Evaluation functions
  */
  /*@{*/

  ///
  /** Evaluate potential on fine points of all faces of destination box
      due to charges on all faces of source box.
      The function setCoeffs() must have been called already.
  */
  void eval(FArrayBox&   a_dst);

  ///
  /** Evaluate potential on fine points of one face of destination box
      due to charges on one face of source box.
      The function setCoeffs() must have been called already for
      the source face.
  */
  void eval(FArrayBox&              a_dstFace,
            int                     a_idirSrc,
            const Side::LoHiSide&   a_sideSrc,
            int                     a_idirDst,
            const Side::LoHiSide&   a_sideDst);

  ///
  /** Evaluate potential on outer coarse points
      due to charges on all source faces.
      The function setCoeffs() must have been called already,
      and m_getOuterCoarse must be true.
  */
  void evalOuterCoarse(FArrayBox&   a_dstOuter);

  ///
  /** Evaluate potential on outer coarse points
      due to charges on one source face.
      The function setCoeffs() must have been called already for
      the source face, and m_getOuterCoarse must be true.

  */
  void evalOuterCoarse(FArrayBox&              a_dstOuter,
                       int                     a_idirSrc,
                       const Side::LoHiSide&   a_sideSrc);

  /*@}*/

protected:

  enum DestFace
  {
    NEAR           = 0,
    FAR            = 1,
    FIRST_OBLIQUE  = 2,
    SECOND_OBLIQUE = 3
  };

  /** running in parallel?
   */
  bool m_parallel;

  /** has this Multipoles been defined yet?
   */
  bool m_isDefined;

  /** NODE-centered source box
   */
  Box m_srcBox;

  /** NODE-centered source box
   */
  Box m_dstBox;

  /** mesh spacing at this level
   */
  RealVect m_dx;

  /** patch size
   */
  int m_patchSize;

  /** length of source face
   */
  IntVect m_n1;

  /** length of destination face
   */
  IntVect m_n2;

  /** length of destination face coarsened by m_dstFaceCoarsening
   */
  IntVect m_n2c;

  /** length of destination face coarsened by m_outerCoarsening
   */
  IntVect m_n2cOuter;

  /** degree of multipole expansions
   */
  int m_multipoleOrder;

  /** number of coefficients in multipole expansions
   */
  int m_numCoeffs;

  /** coarsening ratio for destination faces
   */
  int m_dstFaceCoarsening;

  /** width of border for polynomial interpolation
   */
  int m_interpBorder;

  /** verbosity; default 0
   */
  int m_verbose;

  /** apply direct evaluation instead of multipole expansion?
   */
  bool m_direct;

  /** weights applied to charges on a patch to obtain multipole coefficients
   */
  Tuple<Vector<FArrayBox*>, SpaceDim> m_weights;

  BoxLayout m_coeffsLayout;

  /** coefficients of multipole expansions due to charges on faces of m_srcBox
   */
  BoxLayoutData<FArrayBox> m_coeffsPar;

  Vector<FArrayBox*> m_coeffs;

  /** (if m_direct) charges on faces of m_srcBox, saved in setCoeffs()
   */
  Vector<FArrayBox*> m_charges;

  /** (if m_direct) weights on charges on faces of m_srcBox
   */
  Vector<FArrayBox*> m_chargeWeights;

  /** CELL-centered boxes on which coefficients live
   */
  Tuple<Box, SpaceDim> m_bxPatchesFace;

  /** first and second other directions
      m_dir1[0] = 1, m_dir2[0] = 2
      m_dir1[1] = 0, m_dir2[1] = 2
      m_dir1[2] = 0, m_dir2[2] = 1
   */
  IntVect m_dir1, m_dir2;

  /** lower corner of each patch
   */
  Tuple<BaseFab<IntVect>*, SpaceDim> m_startPatch;

  /** weights for obtaining multipole coefficients from charges
   */
  Tuple<BaseFab<int>*, SpaceDim> m_whichWeights;

  /** extent of each patch in each dimension
   */
  Tuple<IntVect, SpaceDim> m_patchOffset;

  /** Have we computed coefficients due to charges on each face?
      Recall that Vector<bool> doesn't work, so use vector<bool>.
   */
  vector<bool> m_gotCoeffs;

  /** NODE-centered coarse boxes on which we evaluate multipole expansions
   */
  Vector<Box> m_bxDstCoarseFace;

  /** NODE-centered fine boxes on which we evaluate multipole expansions
   */
  Vector<Box> m_bxDstFace;

  /** For interpolation on faces
   */
  Tuple<InterpOnFace, SpaceDim> m_interp;

  /** For interpolation on faces
   */
  Tuple<InterpOnFace, SpaceDim> m_bxFaceCoarse;

  /** solve for coarser-level solution also?
   */
  bool m_getOuterCoarse;

  /** refinement ratio to coarser level on rectangular domain
   */
  int m_outerCoarsening;

  /** width of layer of coarse nodes outside m_dstBox
   */
  int m_outerCoarseBuffer;

  /** box of points coarsened by m_outerCoarsening, on which to
      evaluate multipole expansions outside m_dstBox
   */
  Box m_outerCoarseBox;

  /** number of patches in each dimension on bound(D1)
   */
  IntVect m_numPatches;

  /** for saving...
   */

  /** underlying layout for m_srcdstPar and m_sdcrsePar:
      a box for each source face
   */
  BoxLayout m_srcFaces;

  /** which face ID corresponds to which box
   */
  LayoutData<int> m_faceSrcIDs;

  /** Used only when m_parallel true:
      processors on which faces live; length 2*SpaceDim.
      m_faceSrcProc[faceSrcID] is index of processor on which face
      with ID faceSrcID lives.
   */
  Vector<int> m_faceSrcProc;

  /** Used only when m_parallel true:
      list of all distinct processors in the vector m_faceSrcProc.
   */
  Vector<int> m_faceSrcProcDistinct;

  // Vector<Real*> m_srcdst;
  Tuple<Vector<Real*>, SpaceDim> m_srcdst;

  Vector<Real*> m_srcdstDirect;

  LayoutData< Vector<Real*> > m_srcdstPar;

  // Real* m_sdcrse;
  // Vector<Real*> m_sdcrse;
  Tuple<Real*, SpaceDim> m_sdcrse;

  Vector<Real*> m_sdcrseDirect;

  LayoutData< Real* > m_sdcrsePar;

  int m_nterms;

  int m_ncpts;

  int* m_cpind;

  int* m_cpow;

  Real* m_cfac;

private:

  void plusReduce(FArrayBox&                a_sum,
                  const Vector< RefCountedPtr<FArrayBox> >& a_data);

  void findface(int&              a_idirSrc,
                Side::LoHiSide&   a_sideSrc,
                const Box&        a_bxSrcFace);

  void getCodeFlips(DestFace&               a_dstCode,
                    int&                    a_srcFlip,
                    int&                    a_dstFlip,
                    int                     a_idirSrc,
                    const Side::LoHiSide&   a_sideSrc,
                    int                     a_idirDst,
                    const Side::LoHiSide&   a_sideDst);

  unsigned long sizeMatrixPatches(int   a_idirSrc,
                                  int   a_idirDst);

  unsigned long sizeMatrixDirect(int   a_idirSrc,
                                 int   a_idirDst);

  unsigned long sizeOuterMatrixDirect(int   a_idirSrc);

  unsigned long sizeOuterMatrixPatches(int   a_idirSrc);

  void operator = (const Multipoles&)
  {
  }

  Multipoles(const Multipoles&)
  {
  }
};

#include "NamespaceFooter.H"

#endif
