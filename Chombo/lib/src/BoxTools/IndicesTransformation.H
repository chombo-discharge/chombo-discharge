#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#ifndef _INDICESTRANSFORMATION_H_
#define _INDICESTRANSFORMATION_H_

#include <iostream>
#include "IntVect.H"
#include "RealVect.H"
#include "BaseFab.H"
#include "Tuple.H"
#include "IndicesFunctions.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
/// Class to describe transformation of indices from one block to another.
class IndicesTransformation
{
public:

  /// null constructor leaves object in undefined state.
  IndicesTransformation();

  ///
  /**
     Constructor.

     IntVect pOld:  indices of a cell in old index space
     IntVect pNew:  indices of the same cell in new index space

     Then for each direction idir:
     pNew[idir] = a_sign[idir]*pOld[a_permutation[idir]] + a_translation[idir].
  */
  IndicesTransformation(const IntVect& a_permutation,
                        const IntVect& a_sign,
                        const IntVect& a_translation);

  ///
  /*
    Same as constructor.
   */
  void define(const IntVect& a_permutation,
              const IntVect& a_sign,
              const IntVect& a_translation);

  ///
  /*
    Define an IndicesTransformation that is a swap of two indices.
    Is Identity if the indices are equal.
   */
  void defineFromSwap(int a_ind1,
                      int a_ind2);

  ///
  /*
    Define an IndicesTransformation that is just a translation.
   */
  void defineFromTranslation(const IntVect& a_translation);

  ///
  /*
    Define an IndicesTransformation
    such that the NODE a_pivotOld maps to the NODE a_pivotNew
    and for each idir, [idir] maps to [a_permutation[idir]]
    with sign a_sign[idir].
   */
  void defineFromPivot(const IntVect& a_pivotOld,
                       const IntVect& a_pivotNew,
                       const IntVect& a_permutation,
                       const IntVect& a_sign);

  ///
  /*
    Define an IndicesTransformation
    such that a_srcBox's face in dimension a_srcDim side a_srcSide
    maps to a_dstBox's face in dimension a_dstDim side a_dstSide
    and for each idir, a_sign[idir] is -1 or +1 according to whether
    source direction idir is flipped or not in the opposite block.
   */
  void defineFromFaces(const Box& a_srcBox,
                       int a_srcDim,
                       Side::LoHiSide a_srcSide,
                       const Box& a_dstBox,
                       int a_dstDim,
                       Side::LoHiSide a_dstSide,
                       const IntVect& a_sign = IntVect::Unit);

  ///
  /*
    Is this transformation defined?
   */
  bool isDefined() const;

  bool operator==(const IndicesTransformation& a_itOther) const;

  bool operator!=(const IndicesTransformation& a_itOther) const;

  friend std::ostream& operator<< (std::ostream& a_os,
                                   const IndicesTransformation& a_it);

  ///
  /**
     IntVect pOld:  indices of a cell in old index space
     IntVect pNew:  indices of the same cell in new index space

     Then for each direction idir:
     pNew[idir] = m_sign[idir]*pOld[m_permutation[idir]] + m_translation[idir]
     and hence:
     pOld[m_permutation[idir]] = m_sign[idir]*pNew[idir]
                               - m_sign[idir]*m_translation[idir]
  */
  IntVect transformFwd(const IntVect& a_ivOld) const;

  IntVect transformBack(const IntVect& a_ivNew) const;

  IntVect transform(const IntVect& a_iv, bool a_forward = true) const;

  ///
  /**
     Transform NODE indices.
   */
  IntVect transformNode(const IntVect& a_iv) const;

  ///
  /**
     Transform coordinates in mapped space, not an index.
   */
  RealVect transformMapped(const RealVect& a_pointOld,
                           const RealVect& a_dxOld,
                           const RealVect& a_dxNew) const;

  ///
  /**
     Transform a vector, not an index.
     There may be permutation of indices and sign change, but no translation.
   */
  IntVect transformVectorFwd(const IntVect& a_vecOld) const;

  IntVect transformVectorBack(const IntVect& a_vecNew) const;

  IntVect transformVector(const IntVect& a_vec, bool a_forward = true) const;

  ///
  /**
     Transform indices of either cells or nodes or faces or edges.
   */
  IntVect transformWithType(const IntVect& a_iv,
                            const IntVect& a_tp,
                            bool a_forward = true) const;

  ///
  /**
     Transform the type of Box:
     no change if a_tp is the cell type or the node type,
     but if a face or edge type, then the transformed type may differ.
   */
  IntVect transformType(const IntVect& a_tp, bool a_forward = true) const;

  ///
  /**
     Transform a whole box.
   */
  Box transformFwd(const Box& a_bxOld) const;

  Box transformBack(const Box& a_bxNew) const;

  Box transform(const Box& a_bx, bool a_forward = true) const;

  ///
  /**
     Apply forward transformation on a Box within a BaseFab.
   */
  template<typename T> void transformFwd(BaseFab<T>& a_dstFab,
                                         const BaseFab<T>& a_srcFab,
                                         const Box& a_srcBox,
                                         const Interval& a_dstIntvl,
                                         const Interval& a_srcIntvl) const
  {
    transform(a_dstFab, a_srcFab, a_srcBox,
              a_dstIntvl, a_srcIntvl, true);
  }

  ///
  /**
     Apply inverse transformation on a Box within a BaseFab.
   */
  template<typename T> void transformBack(BaseFab<T>& a_dstFab,
                                          const BaseFab<T>& a_srcFab,
                                          const Box& a_srcBox,
                                          const Interval& a_dstIntvl,
                                          const Interval& a_srcIntvl) const
  {
    transform(a_dstFab, a_srcFab, a_srcBox,
              a_dstIntvl, a_srcIntvl, false);
  }

  ///
  /**
     Apply forward or inverse transformation on a Box within a BaseFab.
   */
  template<typename T> void transform(BaseFab<T>& a_dstFab,
                                      const BaseFab<T>& a_srcFab,
                                      const Box& a_srcBox,
                                      const Interval& a_dstIntvl,
                                      const Interval& a_srcIntvl,
                                      bool a_forward = true) const
  {
#if 0
    CH_assert(a_srcFab.box().contains(a_srcBox));
    Box transformedBox = transform(a_srcBox, a_forward);
    CH_assert(a_dstFab.box().contains(transformedBox));
    IntVect tp = a_srcBox.type();
    int dstComp = a_dstIntvl.begin();
    int srcComp = a_srcIntvl.begin();
    int srcCompHi = a_srcIntvl.end();
    for (; srcComp <= srcCompHi; ++srcComp, ++dstComp)
      {
        for (BoxIterator bit(a_srcBox); bit.ok(); ++bit)
          {
            IntVect ivSrc = bit();
            IntVect ivDst = transformWithType(ivSrc, tp, a_forward);
            a_dstFab(ivDst, dstComp) =
              a_srcFab(ivSrc, srcComp);
          }
      }
#else
    // const IndicesTransformation itThis = (a_forward) ? *this : this->inverse();

  if (a_dstIntvl.size() == 0)
    { // No components to set.
      return;
    }
  if (a_srcBox.isEmpty())
    { // Nothing to set.
      return;
    }
  CH_assert(a_dstIntvl.size() == a_srcIntvl.size());
  // dstIV[idir] = sign[idir]*srcIV[perm[idir]] + translation[idir]
  const Box& dstAllBox = a_dstFab.box();
  const Box& srcAllBox = a_srcFab.box();
  const IntVect& srcAllDims = srcAllBox.size();
  CH_assert(srcAllBox.contains(a_srcBox));
  const IntVect& srcType = a_srcBox.type();
  // Box dstSetBox = boxtransform(a_srcBox, perm, sign, translation);
  Box dstSetBox = transform(a_srcBox, a_forward);
  // Check that dstSetBox, the transformation of a_srcBox, is in dstAllBox.
  CH_assert(dstAllBox.contains(dstSetBox));
  const IntVect& srcSetLo = a_srcBox.smallEnd();
  const IntVect& srcSetDims = a_srcBox.size();
  IntVect dstSetLo = transformWithType(srcSetLo, srcType, a_forward);
  const size_t dstIndexLo =
    FortranArrayIndex(dstSetLo, dstAllBox);
  const size_t srcIndexLo =
    FortranArrayIndex(srcSetLo, srcAllBox);
  // Not unsigned, because increments can be negative.
  Tuple<long long int, SpaceDim> srcInc, dstInc, srcAllInc;
  for (int idir = 0; idir < SpaceDim; idir++)
    { // What happens if you move by BASISV(idir) from origin of set box?
      IntVect src1dir = srcSetLo + BASISV(idir);
      IntVect dst1dir = transformWithType(src1dir, srcType, a_forward);
      if ( ! (dstAllBox.contains(dst1dir)) )
        {
          // FortranArrayIndex(dst1dir, dstAllBox) will be negative.
          dstInc[idir] = 0;
        }
      else
        {
          size_t dstIndex1dirAll =
            FortranArrayIndex(dst1dir, dstAllBox);
          // Need to be careful when subtracting unsigned ints.
          if (dstIndex1dirAll >= dstIndexLo)
            {
              dstInc[idir] = dstIndex1dirAll - dstIndexLo;
            }
          else
            { // dstIndexLo > dstIndex1dirAll
              dstInc[idir] = -(dstIndexLo - dstIndex1dirAll);
            }
        }
      // same as FortranArrayIndex(src1dir, srcAllBox) - srcIndexLo;
      srcInc[idir] =
        FortranArrayIndex(BASISV(idir), srcAllDims);
      srcAllInc[idir] = srcInc[idir] * srcSetDims[idir];
    }

  Tuple<size_t, SpaceDim> srcBegin, dstBegin, srcEnd;
  srcBegin[CH_SPACEDIM-1] = srcIndexLo;
  dstBegin[CH_SPACEDIM-1] = dstIndexLo;

  size_t srcIndex, dstIndex;

  for (int srcComp = a_srcIntvl.begin(), dstComp = a_dstIntvl.begin();
       srcComp <= a_srcIntvl.end();
       ++srcComp, ++dstComp)
    {
      Interval src1Intvl(srcComp, srcComp);
      Interval dst1Intvl(dstComp, dstComp);
      BaseFab<T> dst1Fab(dst1Intvl, a_dstFab);
      const BaseFab<T> src1Fab(src1Intvl, const_cast<BaseFab<T>&>(a_srcFab));
      
      const T* src1Array = src1Fab.dataPtr();
      T* dst1Array = dst1Fab.dataPtr();

#if CH_SPACEDIM>5
      srcEnd[5] = srcBegin[5] + srcAllInc[5];
      for (srcBegin[4] = srcBegin[5], dstBegin[4] = dstBegin[5];
           srcBegin[4] < srcEnd[5];
           srcBegin[4] += srcInc[5], dstBegin[4] += dstInc[5])
        {
#endif
#if CH_SPACEDIM>4
          srcEnd[4] = srcBegin[4] + srcAllInc[4];
          for (srcBegin[3] = srcBegin[4], dstBegin[3] = dstBegin[4];
               srcBegin[3] < srcEnd[4];
               srcBegin[3] += srcInc[4], dstBegin[3] += dstInc[4])
            {
#endif
#if CH_SPACEDIM>3
              srcEnd[3] = srcBegin[3] + srcAllInc[3];
              for (srcBegin[2] = srcBegin[3], dstBegin[2] = dstBegin[3];
                   srcBegin[2] < srcEnd[3];
                   srcBegin[2] += srcInc[3], dstBegin[2] += dstInc[3])
                {
#endif
#if CH_SPACEDIM>2
                  srcEnd[2] = srcBegin[2] + srcAllInc[2];
                  for (srcBegin[1] = srcBegin[2], dstBegin[1] = dstBegin[2];
                       srcBegin[1] < srcEnd[2];
                       srcBegin[1] += srcInc[2], dstBegin[1] += dstInc[2])
                    {
#endif
#if CH_SPACEDIM>1
                      srcEnd[1] = srcBegin[1] + srcAllInc[1];
                      for (srcBegin[0] = srcBegin[1], dstBegin[0] = dstBegin[1];
                           srcBegin[0] < srcEnd[1];
                           srcBegin[0] += srcInc[1], dstBegin[0] += dstInc[1])
                        {
#endif
                          srcEnd[0] = srcBegin[0] + srcAllInc[0];
                          for (srcIndex = srcBegin[0], dstIndex = dstBegin[0];
                               srcIndex < srcEnd[0];
                               srcIndex += srcInc[0], dstIndex += dstInc[0])
                            {
#ifndef NDEBUG
                              // No need to check dstIndex >=0, srcIndex >= 0,
                              // because they are unsigned ints.
                              CH_assert(srcIndex < srcAllBox.numPts());
                              CH_assert(dstIndex < dstAllBox.numPts());
#endif
                              dst1Array[dstIndex] = src1Array[srcIndex];
                            }
#if CH_SPACEDIM>1
                        }
#endif
#if CH_SPACEDIM>2
                    }
#endif
#if CH_SPACEDIM>3
                }
#endif
#if CH_SPACEDIM>4
            }
#endif
#if CH_SPACEDIM>5
        }
#endif
    } // loop over components    
#endif
  }

  ///
  /**
     Return the inverse transformation.
   */
  IndicesTransformation inverse() const;

  ///
  /**
     Return the composite transformation:
     that is, this transformation followed by a_next.
   */
  IndicesTransformation compose(const IndicesTransformation& a_next) const;

  ///
  /**
     Return this transformation with index spaces coarsened by a_refRatio.
   */
  IndicesTransformation coarsen(int a_refRatio) const;

  ///
  /**
     Return this transformation with index spaces refined by a_refRatio.
   */
  IndicesTransformation refine(int a_refRatio) const;

  static const IntVect NoPermutation;

  ///
  /** The identity transformation, which has no permutation,
      no sign change, and no translation.
   */
  static const IndicesTransformation Identity;

  ///
  /** The undefined transformation.
   */
  static const IndicesTransformation Undefined;

  IntVect getPermutation() const
  {
    return m_permutation;
  }
  IntVect getSign() const
  {
    return m_sign;
  }
  IntVect getTranslation() const
  {
    return m_translation;
  }

  ///
  /** Initializes IndicesTransformation::Identity .
   */
  static int InitStatics();

protected:

  // Zero if undefined
  IntVect m_permutation;

  // Zero if undefined
  IntVect m_sign;

  // Zero if undefined
  IntVect m_translation;
};

//
//// Static initialization.  Gotta make sure there are no static object
//// definitions above here (except possibly stuff in headers).  Otherwise,
//// the danger is that some static object's constructor calls
//// IndicesTransformation::Identity, the very things the following
//// definition is supposed to initialize.
//// (I got this from IntVect.H)
static int s_dummyForIndicesTransformationH = IndicesTransformation::InitStatics();

#include "NamespaceFooter.H"

#endif // include guard
