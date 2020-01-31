#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// Multipoles.cpp
// petermc, 25 Aug 2004

#include "parstream.H" // for pout
#include "Multipoles.H"
#include "MultipoleCoeffsF_F.H"
#include "IntegrationWeightsF_F.H"
#include "FindFace.H"
#include "MayDay.H"
#include "Projections.H"
#include "InfiniteDomainF.H"
#include "BoxIterator.H"
#include "LayoutIterator.H"
#include "LoadBalance.H"
#include "EvalDirectF_F.H"
using std::cout;
using std::endl;
using std::cerr;

#include "NamespaceHeader.H"

// ---------------------------------------------------------
unsigned long Multipoles::sizeMatrixPatches(int   a_idirSrc,
                                            int   a_idirDst)
{
  unsigned long sz = m_numCoeffs *
    (unsigned long) (m_n2c[m_dir1[a_idirDst]] + 2*m_interpBorder + 1) *
    (unsigned long) (m_n2c[m_dir2[a_idirDst]] + 2*m_interpBorder + 1) *
    (unsigned long) (m_numPatches[m_dir1[a_idirSrc]] *
                     m_numPatches[m_dir2[a_idirSrc]]);
  return sz;
}


// ---------------------------------------------------------
unsigned long Multipoles::sizeMatrixDirect(int   a_idirSrc,
                                           int   a_idirDst)
{
  Box bxSrcFace = m_srcBox;
  bxSrcFace.setBig(a_idirSrc, m_srcBox.smallEnd(a_idirSrc));
  unsigned long sz =
    (unsigned long) (m_n2c[m_dir1[a_idirDst]] + 2*m_interpBorder + 1) *
    (unsigned long) (m_n2c[m_dir2[a_idirDst]] + 2*m_interpBorder + 1) *
    (unsigned long) (bxSrcFace.numPts());
  return sz;
}


// ---------------------------------------------------------
unsigned long Multipoles::sizeOuterMatrixDirect(int   a_idirSrc)
{
  Box bxSrcFace = m_chargeWeights[a_idirSrc]->box();
  unsigned long sz =
    (unsigned long) m_ncpts *
    (unsigned long) (bxSrcFace.numPts());
  return sz;
}


// ---------------------------------------------------------
unsigned long Multipoles::sizeOuterMatrixPatches(int   a_idirSrc)
{
  unsigned long sz = (unsigned long) (m_numCoeffs * m_ncpts) *
    (unsigned long) (m_numPatches[m_dir1[a_idirSrc]] *
                     m_numPatches[m_dir2[a_idirSrc]]);
  return sz;
}


// ---------------------------------------------------------
void Multipoles::getCodeFlips(DestFace&               a_dstCode,
                              int&                    a_srcFlip,
                              int&                    a_dstFlip,
                              int                     a_idirSrc,
                              const Side::LoHiSide&   a_sideSrc,
                              int                     a_idirDst,
                              const Side::LoHiSide&   a_sideDst)
{
  a_srcFlip = 0;
  a_dstFlip = 0;
  if (a_idirDst == a_idirSrc)
    {
      // no flip needed for source or dest
      if (a_sideDst == a_sideSrc)
        a_dstCode = NEAR;
      else
        a_dstCode = FAR;
    }
  else
    {
      if (a_idirDst == m_dir1[a_idirSrc])
        a_dstCode = FIRST_OBLIQUE;
      else if (a_idirDst == m_dir2[a_idirSrc])
        a_dstCode = SECOND_OBLIQUE;
      else
        {
          cerr << "Multipoles::getCodeFlips has illegal side"
               << endl;
          MayDay::Error("returning");
          return;
        }
      if (! m_parallel)
        { // In non-parallel case, we do not have srcdst stored for
          // high source face.
          // For low source face, no dest flip needed.
          // For high source face, flip dest dimension
          // of source from dest's point of view.
          if (a_sideSrc == Side::Hi)
            {
              if (a_idirSrc == m_dir1[a_idirDst])
                a_dstFlip = 1;
              else if (a_idirSrc == m_dir2[a_idirDst])
                a_dstFlip = 2;
            }
        }
      // For low dest face, no source flip needed.
      // For high dest face, flip source dimension
      // of dest from source's point of view.
      if (a_sideDst == Side::Hi)
        {
          if (a_idirDst == m_dir1[a_idirSrc])
            a_srcFlip = 1;
          else if (a_idirDst == m_dir2[a_idirSrc])
            a_srcFlip = 2;
        }
    }
}


// ---------------------------------------------------------
void Multipoles::plusReduce(FArrayBox&                a_sum,
                            const Vector< RefCountedPtr<FArrayBox> >& a_data)
{
  int ncomp = a_sum.nComp();
  // CH_assert(a_data.nComp() == ncomp);
  // const Box& sumFabBox = a_sum.box();
  a_sum.setVal(0.);
  for (int ind = 0; ind < a_data.size(); ind++)
    {
      const FArrayBox& dataFab = *a_data[ind];
      CH_assert(dataFab.nComp() == ncomp);
      // Add over intersection of domains of a_sum and dataFab.
      a_sum.plus(dataFab, 0, 0, ncomp);
    }
}


// ---------------------------------------------------------
IntVect ceilIntVectDivide(const IntVect& iv,
                          int ref)
{
  IntVect quotient = iv / ref;
  for (int idir = 0; idir < SpaceDim; idir++)
    if (ref * quotient[idir] < iv[idir])
      quotient.setVal(idir, quotient[idir]+1);

  return quotient;
}

// ---------------------------------------------------------
// default constructor
Multipoles::Multipoles()
{
  setDefaultValues();
}


// ---------------------------------------------------------
// complete constructor
Multipoles::Multipoles(const Box&        a_srcBox,
                       const Box&        a_dstBox,
                       const RealVect&   a_dx,
                       int               a_patchSize,
                       int               a_multipoleOrder,
                       int               a_dstFaceCoarsening,
                       int               a_interpBorder,
                       bool              a_parallel,
                       bool              a_getOuterCoarse,
                       int               a_outerCoarsening,
                       int               a_outerCoarseBuffer)
{
  setDefaultValues();
  define(a_srcBox, a_dstBox, a_dx,
         a_patchSize, a_multipoleOrder,
         a_dstFaceCoarsening, a_interpBorder, a_parallel,
         a_getOuterCoarse, a_outerCoarsening, a_outerCoarseBuffer);
}


// ---------------------------------------------------------
void
Multipoles::define(const Box&        a_srcBox,
                   const Box&        a_dstBox,
                   const RealVect&   a_dx,
                   int               a_patchSize,
                   int               a_multipoleOrder,
                   int               a_dstFaceCoarsening,
                   int               a_interpBorder,
                   bool              a_parallel,
                   bool              a_getOuterCoarse,
                   int               a_outerCoarsening,
                   int               a_outerCoarseBuffer)
{
  // pout() << "Proc " << procID() << endl;
  clearMemory();
  m_srcBox = a_srcBox;
  m_dstBox = a_dstBox;
  m_dx = a_dx;
  if (a_multipoleOrder < 0)
    {
      // direct evaluation from source faces to destination faces
      m_direct = true;
    }
  else
    {
      // multipole expansions from source faces to destination faces
      m_direct = false;
      m_patchSize = a_patchSize;
      m_multipoleOrder = a_multipoleOrder;
    }
  m_dstFaceCoarsening = a_dstFaceCoarsening;
  m_interpBorder = a_interpBorder;
  m_parallel = a_parallel;
  // Force parallelism off if m_direct == true.  FIX this!
  if (m_direct) m_parallel = false;
  m_getOuterCoarse = a_getOuterCoarse;
  if (m_getOuterCoarse)
    {
      m_outerCoarsening = a_outerCoarsening;
      m_outerCoarseBuffer = a_outerCoarseBuffer;
    }
  else
    {
      m_outerCoarsening = -1;
      m_outerCoarseBuffer = 0;
    }
  // Force parallelism off if getting an outer coarse layer.  FIX this!
  if (m_getOuterCoarse && m_outerCoarseBuffer > 0) m_parallel = false;

  m_n1 = m_srcBox.bigEnd() - m_srcBox.smallEnd();
  m_n2 = m_dstBox.bigEnd() - m_dstBox.smallEnd();

  // IntVect m_s2 = (m_n2 - m_n1) / 2;

  // m_n2c:  length of coarsen(D2, dstFaceCoarsening)
  // m_dstFaceCoarsening must divide n2....
  // m_n2c = ceil(m_n2 / m_dstFaceCoarsening)
  m_n2c = ceilIntVectDivide(m_n2, m_dstFaceCoarsening);

  m_bxDstFace.resize(2*SpaceDim);
  m_bxDstCoarseFace.resize(2*SpaceDim);

  m_gotCoeffs.resize(2*SpaceDim);

  m_dir1 = IntVect(D_DECL(1, 0, 0));
  m_dir2 = IntVect(D_DECL(2, 2, 1));

  int faceDstID = 0;
  for (SideIterator sit; sit.ok(); sit.next())
    {
      Side::LoHiSide side = sit();
      for (int idirDst = 0; idirDst < SpaceDim; idirDst++)
        {
          m_bxDstFace[faceDstID] = m_dstBox;
          if (side == Side::Lo)
            m_bxDstFace[faceDstID].setBig(idirDst, m_dstBox.smallEnd(idirDst));
          else if (side == Side::Hi)
            m_bxDstFace[faceDstID].setSmall(idirDst, m_dstBox.bigEnd(idirDst));
          else
            MayDay::Error("bad side, returning");

          m_bxDstCoarseFace[faceDstID] =
            coarsen(m_bxDstFace[faceDstID] - m_bxDstFace[faceDstID].smallEnd(),
                    m_dstFaceCoarsening);
          IntVect ivBorder =
            m_interpBorder * (IntVect::Unit - BASISV(idirDst));
          m_bxDstCoarseFace[faceDstID].grow(ivBorder);
          m_gotCoeffs[faceDstID] = false;

          faceDstID++;
        }
    }

  if (m_parallel)
    { // Get m_srcFaces, the layout of source faces.
      Vector<Box> srcFacesVec(2*SpaceDim);
      int faceSrcID = 0;
      for (SideIterator sitSrc; sitSrc.ok(); sitSrc.next())
        {
          Side::LoHiSide sideSrc = sitSrc();
          for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
            {
              Box bxSrcFace = m_srcBox;
              if (sideSrc == Side::Lo)
                bxSrcFace.setBig(idirSrc, m_srcBox.smallEnd(idirSrc));
              else if (sideSrc == Side::Hi)
                bxSrcFace.setSmall(idirSrc, m_srcBox.bigEnd(idirSrc));
              else
                MayDay::Error("bad side, returning");

              srcFacesVec[faceSrcID] = bxSrcFace;
              faceSrcID++;
            }
        }

      // Vector<int> m_faceSrcProc;
      // m_faceSrcProc[faceSrcID] is the processor holding face with
      // ID faceSrcID.
      m_faceSrcProc.resize(2*SpaceDim);
      int eekflag = LoadBalance(m_faceSrcProc, srcFacesVec);
      if (eekflag != 0)
        {
          MayDay::Error("Error in LoadBalance");
        }
      m_srcFaces.define(srcFacesVec, m_faceSrcProc);

      // m_faceSrcProcDistinct is Vector listing IDs of all distinct
      // processors in the Vector m_faceSrcProc.
      m_faceSrcProcDistinct.clear();
      for (int iproc = 0; iproc < 2*SpaceDim; iproc++)
        {
          int thisProc = m_faceSrcProc[iproc];
          bool alreadyHave = false;
          for (int iprocOld = 0; iprocOld < iproc; iprocOld++)
            if (m_faceSrcProc[iprocOld] == thisProc)
              alreadyHave = true;

          if (! alreadyHave) m_faceSrcProcDistinct.push_back(thisProc);
        }

      // Get m_faceSrcIDs, which tells you which face m_srcFaces is on.
      // This is a very kludgy way to do it.
      m_faceSrcIDs.define(m_srcFaces);
      for (DataIterator srcDit = m_srcFaces.dataIterator();
           srcDit.ok(); ++srcDit)
        {
          const Box& thisFace = m_srcFaces.get(srcDit());
          int& thisFaceID = m_faceSrcIDs[srcDit()];
          int faceSrcID = 0;
          for (SideIterator sitSrc; sitSrc.ok(); sitSrc.next())
            {
              Side::LoHiSide sideSrc = sitSrc();
              for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
                {
                  Box bxSrcFace = m_srcBox;
                  if (sideSrc == Side::Lo)
                    bxSrcFace.setBig(idirSrc, m_srcBox.smallEnd(idirSrc));
                  else if (sideSrc == Side::Hi)
                    bxSrcFace.setSmall(idirSrc, m_srcBox.bigEnd(idirSrc));
                  else
                    MayDay::Error("bad side, returning");

                  if (bxSrcFace == thisFace)
                    {
                      thisFaceID = faceSrcID;
                      break;
                    }
                  faceSrcID++;
                }
            }
        }
    }

  if (m_direct)
    {
      m_charges.resize(2*SpaceDim);
      m_chargeWeights.resize(2*SpaceDim);

      int faceSrcID = 0;
      for (SideIterator sitSrc; sitSrc.ok(); sitSrc.next())
        {
          Side::LoHiSide sideSrc = sitSrc();
          for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
            {
              Box bxSrcFace = m_srcBox;
              if (sideSrc == Side::Lo)
                bxSrcFace.setBig(idirSrc, m_srcBox.smallEnd(idirSrc));
              else if (sideSrc == Side::Hi)
                bxSrcFace.setSmall(idirSrc, m_srcBox.bigEnd(idirSrc));
              else
                MayDay::Error("bad side, returning");

              m_chargeWeights[faceSrcID] = new FArrayBox(bxSrcFace, 1);

              IntVect lengthSrcFace =
                bxSrcFace.bigEnd() - bxSrcFace.smallEnd();

              // Decide which integration rule to use.
              // Depends on the length over which you are integrating.
              int length1 = lengthSrcFace[m_dir1[idirSrc]];
              int intRule1;
              FORT_BESTINTEGRATIONRULE(CHF_INT(intRule1),
                                       CHF_CONST_INT(length1));

              int length2 = lengthSrcFace[m_dir2[idirSrc]];
              int intRule2;
              FORT_BESTINTEGRATIONRULE(CHF_INT(intRule2),
                                       CHF_CONST_INT(length2));
              IntVect intRules =
                intRule1 * BASISV(m_dir1[idirSrc]) +
                intRule2 * BASISV(m_dir2[idirSrc]);

              FArrayBox& integWeights = *m_chargeWeights[faceSrcID];
              FORT_INTEGRATIONWEIGHTS(CHF_FRA1(integWeights, 0),
                                      CHF_BOX(bxSrcFace),
                                      CHF_CONST_REALVECT(m_dx),
                                      CHF_CONST_INTVECT(intRules));

              // Real dA = 1.;
              // for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
              // if (idirOther != idirSrc)
              // dA *= m_dx[idirOther];
              // m_chargeWeights[faceSrcID]->setVal(dA);

              faceSrcID++;
            }
        }

      // The only purpose of all this mess here is to calculate szsrcdstTotal.
      unsigned long szsrcdstTotal = 0;
      for (int idirDst = 0; idirDst < SpaceDim; idirDst++)
        for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
          szsrcdstTotal += 4 * sizeMatrixDirect(idirSrc, idirDst);

      // hm, check that this isn't more than 2^32
      pout() << "m_srcdstDirect using " << (szsrcdstTotal) << " scalars, "
           << ((Real(szsrcdstTotal) * Real(sizeof(Real)))/ 1048576.) << " MB"
           << endl;
      // m_srcdstDirect = new Real[szsrcdst];

      RealVect dxCoarse = m_dstFaceCoarsening * m_dx;

      m_srcdstDirect.resize(6*6);
      int ori = 0;
      for (SideIterator sitDst; sitDst.ok(); sitDst.next())
        {
          Side::LoHiSide sideDst = sitDst();
          int isideDst = sign(sideDst);
          for (int idirDst = 0; idirDst < SpaceDim; idirDst++)
            {
              int ddirDst = isideDst * (idirDst + 1);

              int faceDstID = (sideDst == Side::Lo) ?
                idirDst : (idirDst + SpaceDim);

              RealVect baseDstCoarse =
                m_bxDstFace[faceDstID].smallEnd() * m_dx;

              for (SideIterator sitSrc; sitSrc.ok(); sitSrc.next())
                {
                  Side::LoHiSide sideSrc = sitSrc();
                  int isideSrc = sign(sideSrc);
                  for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
                    {
                      int ddirSrc = isideSrc * (idirSrc + 1);

                      int faceSrcID = (sideSrc == Side::Lo) ?
                        idirSrc : (idirSrc + SpaceDim);

                      unsigned long szsrcdst =
                        sizeMatrixDirect(idirSrc, idirDst);

                      Box bxSrcFace = m_srcBox;
                      if (sideSrc == Side::Lo)
                        bxSrcFace.setBig(idirSrc, m_srcBox.smallEnd(idirSrc));
                      else if (sideSrc == Side::Hi)
                        bxSrcFace.setSmall(idirSrc, m_srcBox.bigEnd(idirSrc));
                      else
                        MayDay::Error("bad side, returning");

                      RealVect baseSrc =
                        bxSrcFace.smallEnd() * m_dx;

                      IntVect lengthSrcFace =
                        m_srcBox.bigEnd() - m_srcBox.smallEnd();
                      lengthSrcFace.setVal(idirSrc, 0);

                      m_srcdstDirect[ori] =
                        (Real *) malloc(szsrcdst * sizeof(Real));
                      int errcode = 0;
                      FORT_GETDIRECTRECTMATRIX(m_srcdstDirect[ori],
                                               m_chargeWeights[faceSrcID]->dataPtr(),
                                               &ddirDst, &ddirSrc,
                                               lengthSrcFace.dataPtr(),
                                               &m_n2c[m_dir1[idirDst]],
                                               &m_n2c[m_dir2[idirDst]],
                                               &m_interpBorder,
                                               m_dx.dataPtr(),
                                               dxCoarse.dataPtr(),
                                               baseSrc.dataPtr(),
                                               baseDstCoarse.dataPtr(),
                                               &m_verbose, &errcode);
                      if (errcode != 0)
                        {
                          cerr << "getallsrcdst returned error code "
                               << errcode << endl;
                          MayDay::Error("returning");
                          return;
                        }

                      ori++;

                    }
                }
            }
        }

    }
  else // not direct:  using multipole expansions
    {
      // get the srcdst matrix

      int ntmax = ( (m_multipoleOrder+1)*(m_multipoleOrder+2)*
                    (m_multipoleOrder+3)*(m_multipoleOrder+4) )/24;
      // m_cpind = new int[ntmax];
      m_cpind = (int *) malloc(ntmax * sizeof(int));
      // m_cpow = new int[3*ntmax];
      m_cpow = (int *) malloc(3*ntmax * sizeof(int));
      // m_cfac = new Real[ntmax];
      m_cfac = (Real *) (malloc(ntmax * sizeof(Real)));
      FORT_DERIV1OVERR(m_cpind, m_cfac, m_cpow, &m_nterms, &m_multipoleOrder);

      // m_n2c = m_n2 / m_dstFaceCoarsening;
      // for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
      // if (m_n2c[idirSrc] * m_dstFaceCoarsening < m_n2[idirSrc])
      // m_n2c.setVal(idirSrc, m_n2c[idirSrc]+1);

      // CH_assert(m_n2c * m_dstFaceCoarsening == m_n2);

      // m_numPatches:  number of patches in each dimension on bound(D1).
      // It's OK if m_patchSize does not divide m_n1.
      m_numPatches = ceilIntVectDivide(m_n1, m_patchSize);
      // m_numPatches = m_n1 / m_patchSize;
      // for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
      // if (m_numPatches[idirSrc] * m_patchSize < m_n1[idirSrc])
      // m_numPatches.setVal(idirSrc, m_numPatches[idirSrc]+1);

      m_numCoeffs = ((m_multipoleOrder+1)*(m_multipoleOrder+2))/2;

      // Tuple<Vector<FArrayBox*>, SpaceDim> m_weights;
      // weights from charges to multipole coefficients
      for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
        {
          IntVect ivPatchesHigh =
            (m_numPatches-1) * (BASISV(m_dir1[idirSrc]) +
                                BASISV(m_dir2[idirSrc]));
          m_bxPatchesFace[idirSrc] = Box(IntVect::Zero, ivPatchesHigh);

          // Tuple m_startPatch<BaseFab<IntVect>*, SpaceDim>;
          m_startPatch[idirSrc] =
            new BaseFab<IntVect>(m_bxPatchesFace[idirSrc], 1);

          BaseFab<IntVect>& startPatchDir = *m_startPatch[idirSrc];
          for (BoxIterator bit(m_bxPatchesFace[idirSrc]); bit.ok(); ++bit)
            {
              IntVect patch = bit();
              int start1, start2;

              if (patch[m_dir1[idirSrc]] == 0)
                start1 = 0;
              else
                {
                  int endSize1 = (m_numPatches[m_dir1[idirSrc]] == 1) ?
                    m_n1[m_dir1[idirSrc]] :
                    m_patchSize -
                    (m_numPatches[m_dir1[idirSrc]] * m_patchSize -
                     m_n1[m_dir1[idirSrc]]) / 2;
                  start1 = endSize1 +
                    (patch[m_dir1[idirSrc]] - 1) * m_patchSize;
                }

              if (patch[m_dir2[idirSrc]] == 0)
                start2 = 0;
              else
                {
                  int endSize2 = (m_numPatches[m_dir2[idirSrc]] == 1) ?
                    m_n1[m_dir2[idirSrc]] :
                    m_patchSize -
                    (m_numPatches[m_dir2[idirSrc]] * m_patchSize -
                     m_n1[m_dir2[idirSrc]]) / 2;
                  start2 = endSize2 +
                    (patch[m_dir2[idirSrc]] - 1) * m_patchSize;
                }

              startPatchDir(patch, 0) =
                start1 * BASISV(m_dir1[idirSrc]) +
                start2 * BASISV(m_dir2[idirSrc]);
            }

          m_whichWeights[idirSrc] =
            new BaseFab<int>(m_bxPatchesFace[idirSrc], 1);

          // extent of a patch:  will be used later
          m_patchOffset[idirSrc] =
            m_patchSize * (BASISV(m_dir1[idirSrc]) +
                           BASISV(m_dir2[idirSrc]));

          int* patchLength1;
          int* patchLength2;
          if (m_numPatches * m_patchSize == m_n1)
            {
              m_weights[idirSrc].resize(1);

              patchLength1 = new int[1];
              patchLength2 = new int[1];
              patchLength1[0] = m_patchSize; patchLength2[0] = m_patchSize;

              m_whichWeights[idirSrc]->setVal(0);
            }
          else
            {
              // more complicated:  some patches cover fewer points.
              m_weights[idirSrc].resize(4);

              int endSize1 = (m_numPatches[m_dir1[idirSrc]] == 1) ?
                m_n1[m_dir1[idirSrc]] :
                m_patchSize -
                (m_numPatches[m_dir1[idirSrc]] * m_patchSize -
                 m_n1[m_dir1[idirSrc]]) / 2;
              int endSize2 = (m_numPatches[m_dir2[idirSrc]] == 1) ?
                m_n1[m_dir2[idirSrc]] :
                m_patchSize -
                (m_numPatches[m_dir2[idirSrc]] * m_patchSize -
                 m_n1[m_dir2[idirSrc]]) / 2;
              patchLength1 = new int[4];
              patchLength2 = new int[4];
              patchLength1[0] = m_patchSize; patchLength2[0] = m_patchSize;
              patchLength1[1] = endSize1;    patchLength2[1] = m_patchSize;
              patchLength1[2] = m_patchSize; patchLength2[2] = endSize2;
              patchLength1[3] = endSize1;    patchLength2[3] = endSize2;

              // m_bxPatchesFace[idirSrc] contains indices of the patches
              // on face in direction idirSrc.
              // Set edge{1|2}{hi|lo} to contain the indices of the
              // boxes on the {high|low} edges in dimension isrc{1|2}.

              Box edge1lo = m_bxPatchesFace[idirSrc];
              edge1lo.setBig(m_dir1[idirSrc],
                             m_bxPatchesFace[idirSrc].smallEnd(m_dir1[idirSrc]));
              Box edge1hi = m_bxPatchesFace[idirSrc];
              edge1hi.setSmall(m_dir1[idirSrc],
                               m_bxPatchesFace[idirSrc].bigEnd(m_dir1[idirSrc]));

              Box edge2lo = m_bxPatchesFace[idirSrc];
              edge2lo.setBig(m_dir2[idirSrc],
                             m_bxPatchesFace[idirSrc].smallEnd(m_dir2[idirSrc]));
              Box edge2hi = m_bxPatchesFace[idirSrc];
              edge2hi.setSmall(m_dir2[idirSrc],
                               m_bxPatchesFace[idirSrc].bigEnd(m_dir2[idirSrc]));

              m_whichWeights[idirSrc]->setVal(0);
              m_whichWeights[idirSrc]->setVal(1, edge1lo, 0);
              m_whichWeights[idirSrc]->setVal(1, edge1hi, 0);
              m_whichWeights[idirSrc]->setVal(2, edge2lo, 0);
              m_whichWeights[idirSrc]->setVal(2, edge2hi, 0);
              m_whichWeights[idirSrc]->setVal(3, edge1lo & edge2lo, 0);
              m_whichWeights[idirSrc]->setVal(3, edge1lo & edge2hi, 0);
              m_whichWeights[idirSrc]->setVal(3, edge1hi & edge2lo, 0);
              m_whichWeights[idirSrc]->setVal(3, edge1hi & edge2hi, 0);
            }

          //   In this illustration, npchs = 3.
          //   Patch centers are indicated by X.
          //
          //          :<-- ptchsz  -->:               :<-- ptchsz  -->:
          //              :<- endsz ->:<-- ptchsz  -->:<- endsz ->:   :
          //          +---+-----------+---------------+-----------+---+
          //          |   :           :               :           :   |
          //    - - - + - +-----------+---------------+-----------+ - | - - -
          //      ^   |   |           |               |           |   |   ^
          //      |   |   |           |               |           |   |   |
          //    endsz |   |     X     |       X       |     X     |   |   |
          //      |   |   |  shape 2  |    shape 3    |  shape 2  |   |   |
          //      v   |   |           |               |           |   |   |
          //    - - - + - +-----------+---------------+-----------+ - +   |
          //      ^   |   |           |               |           |   |   |
          //      |   |   |           |               |           |   |   |
          //      |   |   |           |               |           |   |   |
          //   ptchsz |   |     X     |       X       |     X     |   |   n1
          //      |   |   |  shape 1  |    shape 0    |  shape 1  |   |   |
          //      |   |   |           |               |           |   |   |
          //      v   |   |           |               |           |   |   |
          //    - - - + - +-----------+---------------+-----------+ - +   |
          //      ^   |   |           |               |           |   |   |
          //      |   |   |           |               |           |   |   |
          //    endsz |   |     X     |       X       |     X     |   |   |
          //      |   |   |  shape 2  |    shape 3    |  shape 2  |   |   |
          //      v   |   |           |               |           |   |   v
          //    - - - + - +-----------+---------------+-----------+ - + - - -
          //          |   :           :               :           :   |
          //          +---+-----------+---------------+-----------+---+
          //              :                                       :
          //              :<--------------- n1 ------------------>:
          //              :                                       :

          for (int shape = 0; shape < m_weights[idirSrc].size(); shape++)
            {
              // Decide which integration rule to use.
              // Depends on the length over which you are integrating.

              // Note that if n1 = 20 and r = 8
              // then we split 20 into patches as 6 + 8 + 6.
              // But 6 is not divisible by 4, so we must use Simpson's rule.
              // Does this have an effect?
              int intRule1;
              FORT_BESTINTEGRATIONRULE(CHF_INT(intRule1),
                                       CHF_CONST_INT(patchLength1[shape]));
              int intRule2;
              FORT_BESTINTEGRATIONRULE(CHF_INT(intRule2),
                                       CHF_CONST_INT(patchLength2[shape]));
              IntVect intRules =
                intRule1 * BASISV(m_dir1[idirSrc]) +
                intRule2 * BASISV(m_dir2[idirSrc]);

              IntVect patchMax =
                patchLength1[shape] * BASISV(m_dir1[idirSrc]) +
                patchLength2[shape] * BASISV(m_dir2[idirSrc]);
              // underlying box of patch for this particular shape
              Box bxPatchFace(IntVect::Zero, patchMax,
                              IndexType::TheNodeType());
              FArrayBox* weightsFacePtr =
                new FArrayBox(bxPatchFace, m_numCoeffs);
              m_weights[idirSrc][shape] = weightsFacePtr;
              RealVect center =
                (bxPatchFace.bigEnd() - bxPatchFace.smallEnd()) * m_dx / 2.;
              FArrayBox weightsIntegration(bxPatchFace, 1);
              FORT_INTEGRATIONWEIGHTS(CHF_FRA1(weightsIntegration, 0),
                                      CHF_BOX(bxPatchFace),
                                      CHF_CONST_REALVECT(m_dx),
                                      CHF_CONST_INTVECT(intRules));
              FORT_WEIGHTCOEFFS(CHF_FRA((*weightsFacePtr)),
                                CHF_CONST_INT(m_multipoleOrder),
                                CHF_CONST_REALVECT(m_dx),
                                CHF_CONST_REALVECT(center),
                                CHF_BOX(bxPatchFace),
                                CHF_CONST_FRA1(weightsIntegration, 0));
            }
        }

      if (m_parallel)
        {
          m_srcdstPar.define(m_srcFaces);
          for (DataIterator srcDit = m_srcFaces.dataIterator();
               srcDit.ok(); ++srcDit)
            {
              // const Box& bxSrcFace = m_srcFaces.get(srcDit());
              Vector<Real*>& srcdstVec = m_srcdstPar[srcDit()];
              // srcdstVec.resize(6);
              srcdstVec.resize(4);
              // m_srcdstPar.resize(6*6);
              // int ori = 0;

              int faceSrcID = m_faceSrcIDs[srcDit()];
              Side::LoHiSide sideSrc = (faceSrcID < SpaceDim) ?
                Side::Lo : Side::Hi;
              int idirSrc = faceSrcID % SpaceDim;

              int isideSrc = sign(sideSrc);
              int ddirSrc = isideSrc * (idirSrc + 1);

              // Calculate and write out szsrcdstTotal.
              unsigned long szsrcdstTotal =
                // for NEAR and FAR
                2 * sizeMatrixPatches(idirSrc, idirSrc) +
                // for FIRST_OBLIQUE
                sizeMatrixPatches(idirSrc, m_dir1[idirSrc]) +
                // for SECOND_OBLIQUE
                sizeMatrixPatches(idirSrc, m_dir2[idirSrc]);

              // hm, check that this isn't more than 2^32
              pout() << "m_srcdstPar using " << (szsrcdstTotal)
                     << " scalars, "
                     << ((Real(szsrcdstTotal) * Real(sizeof(Real)))/ 1048576.) << " MB"
                     << endl;

              // m_srcdstPar = new Real[szsrcdst];

              // For the parallel faces:  NEAR and FAR.
              int idirDst = idirSrc;

              unsigned long szsrcdstParallel =
                sizeMatrixPatches(idirSrc, idirDst);

              int errcode = 0;

              // NEAR face
              srcdstVec[NEAR] =
                (Real *) malloc(szsrcdstParallel * sizeof(Real));
              int ddirDst = isideSrc * (idirDst + 1);
              FORT_GETRECTMATRIX(srcdstVec[NEAR],
                                 &ddirDst, &ddirSrc,
                                 m_n1.dataPtr(), m_n2.dataPtr(),
                                 &m_dstFaceCoarsening, &m_patchSize,
                                 &m_numPatches[m_dir1[idirSrc]],
                                 &m_numPatches[m_dir2[idirSrc]],
                                 &m_n2c[m_dir1[idirDst]],
                                 &m_n2c[m_dir2[idirDst]],
                                 &m_multipoleOrder, &m_interpBorder,
                                 m_dx.dataPtr(),
                                 m_cpind, m_cfac, m_cpow, &m_nterms,
                                 &m_verbose, &errcode);
              if (errcode != 0)
                {
                  cerr << "getrectmatrix returned error code "
                       << errcode << endl;
                  MayDay::Error("returning");
                  return;
                }

              // FAR face
              srcdstVec[FAR] =
                (Real *) malloc(szsrcdstParallel * sizeof(Real));
              ddirDst = -isideSrc * (idirDst + 1);
              FORT_GETRECTMATRIX(srcdstVec[FAR],
                                 &ddirDst, &ddirSrc,
                                 m_n1.dataPtr(), m_n2.dataPtr(),
                                 &m_dstFaceCoarsening, &m_patchSize,
                                 &m_numPatches[m_dir1[idirSrc]],
                                 &m_numPatches[m_dir2[idirSrc]],
                                 &m_n2c[m_dir1[idirDst]],
                                 &m_n2c[m_dir2[idirDst]],
                                 &m_multipoleOrder, &m_interpBorder,
                                 m_dx.dataPtr(),
                                 m_cpind, m_cfac, m_cpow, &m_nterms,
                                 &m_verbose, &errcode);
              if (errcode != 0)
                {
                  cerr << "getrectmatrix returned error code "
                       << errcode << endl;
                  MayDay::Error("returning");
                  return;
                }

              // FIRST_OBLIQUE face
              idirDst = m_dir1[idirSrc];
              unsigned long szsrcdstFirstOblique =
                sizeMatrixPatches(idirSrc, idirDst);
              srcdstVec[FIRST_OBLIQUE] =
                (Real *) malloc(szsrcdstFirstOblique * sizeof(Real));
              ddirDst = -(idirDst + 1);
              FORT_GETRECTMATRIX(srcdstVec[FIRST_OBLIQUE],
                                 &ddirDst, &ddirSrc,
                                 m_n1.dataPtr(), m_n2.dataPtr(),
                                 &m_dstFaceCoarsening, &m_patchSize,
                                 &m_numPatches[m_dir1[idirSrc]],
                                 &m_numPatches[m_dir2[idirSrc]],
                                 &m_n2c[m_dir1[idirDst]],
                                 &m_n2c[m_dir2[idirDst]],
                                 &m_multipoleOrder, &m_interpBorder,
                                 m_dx.dataPtr(),
                                 m_cpind, m_cfac, m_cpow, &m_nterms,
                                 &m_verbose, &errcode);
              if (errcode != 0)
                {
                  cerr << "getrectmatrix returned error code "
                       << errcode << endl;
                  MayDay::Error("returning");
                  return;
                }

              // SECOND_OBLIQUE face
              idirDst = m_dir2[idirSrc];
              unsigned long szsrcdstSecondOblique =
                sizeMatrixPatches(idirSrc, idirDst);
              srcdstVec[SECOND_OBLIQUE] =
                (Real *) malloc(szsrcdstSecondOblique * sizeof(Real));
              ddirDst = -(idirDst + 1);
              FORT_GETRECTMATRIX(srcdstVec[SECOND_OBLIQUE],
                                 &ddirDst, &ddirSrc,
                                 m_n1.dataPtr(), m_n2.dataPtr(),
                                 &m_dstFaceCoarsening, &m_patchSize,
                                 &m_numPatches[m_dir1[idirSrc]],
                                 &m_numPatches[m_dir2[idirSrc]],
                                 &m_n2c[m_dir1[idirDst]],
                                 &m_n2c[m_dir2[idirDst]],
                                 &m_multipoleOrder, &m_interpBorder,
                                 m_dx.dataPtr(),
                                 m_cpind, m_cfac, m_cpow, &m_nterms,
                                 &m_verbose, &errcode);
              if (errcode != 0)
                {
                  cerr << "getrectmatrix returned error code "
                       << errcode << endl;
                  MayDay::Error("returning");
                  return;
                }

            }
          {
            // m_coeffsPar should be distributed, too.
            // m_coeffsPar should live on both high and low faces.
            // m_coeffsPar.resize(2*SpaceDim);
            FindFace tfmFace(m_srcBox, m_bxPatchesFace);
            m_coeffsLayout.deepCopy(m_srcFaces);
            m_coeffsLayout.transform(tfmFace);
            m_coeffsLayout.closeNoSort();

            m_coeffsPar.define(m_coeffsLayout, m_numCoeffs);
          }
        }
      else // ! m_parallel
        {
          // The only purpose of this mess here is to calculate szsrcdstTotal.
          unsigned long szsrcdstTotal = 0;
          for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
            {
              szsrcdstTotal +=
                // for NEAR and FAR
                2 * sizeMatrixPatches(idirSrc, idirSrc) +
                // for FIRST_OBLIQUE
                sizeMatrixPatches(idirSrc, m_dir1[idirSrc]) +
                // for SECOND_OBLIQUE
                sizeMatrixPatches(idirSrc, m_dir2[idirSrc]);
            }

          // hm, check that this isn't more than 2^32
          pout() << "m_srcdst using "
                 << (szsrcdstTotal) << " scalars, "
                 << ((Real(szsrcdstTotal) * Real(sizeof(Real)))/ 1048576.) << " MB"
                 << endl;
          // m_srcdst = new Real[szsrcdst];

          // m_srcdst.resize(6*6);
          Side::LoHiSide sideSrc = Side::Lo;
          for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
            {
              Vector<Real*>& srcdstSrc = m_srcdst[idirSrc];
              srcdstSrc.resize(4);

              int isideSrc = sign(sideSrc);
              int ddirSrc = isideSrc * (idirSrc + 1);

              // For the parallel faces:  NEAR and FAR.

              int idirDst = idirSrc;
              unsigned long szsrcdstParallel =
                sizeMatrixPatches(idirSrc, idirDst);

              int errcode = 0;

              // NEAR face
              srcdstSrc[NEAR] =
                (Real *) malloc(szsrcdstParallel * sizeof(Real));
              int ddirDst = -(idirDst + 1);
              FORT_GETRECTMATRIX(srcdstSrc[NEAR],
                                 &ddirDst, &ddirSrc,
                                 m_n1.dataPtr(), m_n2.dataPtr(),
                                 &m_dstFaceCoarsening, &m_patchSize,
                                 &m_numPatches[m_dir1[idirSrc]],
                                 &m_numPatches[m_dir2[idirSrc]],
                                 &m_n2c[m_dir1[idirDst]],
                                 &m_n2c[m_dir2[idirDst]],
                                 &m_multipoleOrder, &m_interpBorder,
                                 m_dx.dataPtr(),
                                 m_cpind, m_cfac, m_cpow, &m_nterms,
                                 &m_verbose, &errcode);
              if (errcode != 0)
                {
                  cerr << "getrectmatrix returned error code "
                       << errcode << endl;
                  MayDay::Error("returning");
                  return;
                }

              // FAR face
              srcdstSrc[FAR] =
                (Real *) malloc(szsrcdstParallel * sizeof(Real));
              ddirDst = (idirDst + 1);
              FORT_GETRECTMATRIX(srcdstSrc[FAR],
                                 &ddirDst, &ddirSrc,
                                 m_n1.dataPtr(), m_n2.dataPtr(),
                                 &m_dstFaceCoarsening, &m_patchSize,
                                 &m_numPatches[m_dir1[idirSrc]],
                                 &m_numPatches[m_dir2[idirSrc]],
                                 &m_n2c[m_dir1[idirDst]],
                                 &m_n2c[m_dir2[idirDst]],
                                 &m_multipoleOrder, &m_interpBorder,
                                 m_dx.dataPtr(),
                                 m_cpind, m_cfac, m_cpow, &m_nterms,
                                 &m_verbose, &errcode);
              if (errcode != 0)
                {
                  cerr << "getrectmatrix returned error code "
                       << errcode << endl;
                  MayDay::Error("returning");
                  return;
                }

              // FIRST_OBLIQUE face
              idirDst = m_dir1[idirSrc];
              unsigned long szsrcdstFirstOblique =
                sizeMatrixPatches(idirSrc, idirDst);
              srcdstSrc[FIRST_OBLIQUE] =
                (Real *) malloc(szsrcdstFirstOblique * sizeof(Real));
              ddirDst = -(idirDst + 1);
              FORT_GETRECTMATRIX(srcdstSrc[FIRST_OBLIQUE],
                                 &ddirDst, &ddirSrc,
                                 m_n1.dataPtr(), m_n2.dataPtr(),
                                 &m_dstFaceCoarsening, &m_patchSize,
                                 &m_numPatches[m_dir1[idirSrc]],
                                 &m_numPatches[m_dir2[idirSrc]],
                                 &m_n2c[m_dir1[idirDst]],
                                 &m_n2c[m_dir2[idirDst]],
                                 &m_multipoleOrder, &m_interpBorder,
                                 m_dx.dataPtr(),
                                 m_cpind, m_cfac, m_cpow, &m_nterms,
                                 &m_verbose, &errcode);
              if (errcode != 0)
                {
                  cerr << "getrectmatrix returned error code "
                       << errcode << endl;
                  MayDay::Error("returning");
                  return;
                }

              // SECOND_OBLIQUE face
              idirDst = m_dir2[idirSrc];
              unsigned long szsrcdstSecondOblique =
                sizeMatrixPatches(idirSrc, idirDst);
              srcdstSrc[SECOND_OBLIQUE] =
                (Real *) malloc(szsrcdstSecondOblique * sizeof(Real));
              ddirDst = -(idirDst + 1);
              FORT_GETRECTMATRIX(srcdstSrc[SECOND_OBLIQUE],
                                 &ddirDst, &ddirSrc,
                                 m_n1.dataPtr(), m_n2.dataPtr(),
                                 &m_dstFaceCoarsening, &m_patchSize,
                                 &m_numPatches[m_dir1[idirSrc]],
                                 &m_numPatches[m_dir2[idirSrc]],
                                 &m_n2c[m_dir1[idirDst]],
                                 &m_n2c[m_dir2[idirDst]],
                                 &m_multipoleOrder, &m_interpBorder,
                                 m_dx.dataPtr(),
                                 m_cpind, m_cfac, m_cpow, &m_nterms,
                                 &m_verbose, &errcode);
              if (errcode != 0)
                {
                  cerr << "getrectmatrix returned error code "
                       << errcode << endl;
                  MayDay::Error("returning");
                  return;
                }
            }
          m_coeffs.resize(2*SpaceDim);
        }
    }

  /*
    m_interp for interpolation of multipole expansions on faces
  */
  for (int idirDst = 0; idirDst < SpaceDim; idirDst++)
    {
      // faceBox has coordinate of normal direction fixed to zero.
      Box faceBox(m_dstBox);
      faceBox.setRange(idirDst, 0, 1);
      m_interp[idirDst].define(faceBox, idirDst,
                               m_dstFaceCoarsening, m_interpBorder);
    }

  if (m_getOuterCoarse)
    {
      // assumptions:  isotropic, cubic
      Box coarseDstBox = coarsenInner(m_dstBox, m_outerCoarsening);
      m_outerCoarseBox = grow(coarseDstBox, m_outerCoarseBuffer);
      m_ncpts = m_outerCoarseBox.numPts() - coarseDstBox.numPts();
      m_n2cOuter = coarseDstBox.bigEnd() - coarseDstBox.smallEnd();

      // FORT_GETNUMCOARSEPOINTS(&m_ncpts, m_n2.dataPtr(),
      // &m_outerCoarsening, &m_outerCoarseBuffer);

      if (m_outerCoarseBuffer > 0)
        {
        if (m_direct)
          {
            m_sdcrseDirect.resize(6);
            // hm, check that this isn't more than 2^32
            unsigned long szsdcrseTotal = 0;
            for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
              {
                szsdcrseTotal +=
                  2 * sizeOuterMatrixDirect(idirSrc);
              }
            pout() << "m_sdcrseDirect using " << szsdcrseTotal << " scalars, "
                   << ((Real(szsdcrseTotal) * sizeof(Real)) / 1048576.) << " MB"
                   << endl;

            int faceSrcID = 0;
            for (SideIterator sitSrc; sitSrc.ok(); sitSrc.next())
              {
                // Side::LoHiSide sideSrc = sitSrc();
                for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
                  {
                    unsigned long szsdcrse =
                      sizeOuterMatrixDirect(idirSrc);

                    m_sdcrseDirect[faceSrcID] =
                      (Real *) malloc(szsdcrse * sizeof(Real));

                    Box srcFaceBox = m_chargeWeights[faceSrcID]->box();

                    RealVect baseSrc = srcFaceBox.smallEnd() * m_dx;
                    RealVect dxCoarse = m_outerCoarsening * m_dx;
                    IntVect offset; // same as offset in getrectcoarsematrix
                    for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
                      offset[idirOther] =
                        ((m_n2[idirOther] / 2) % m_outerCoarsening) /
                        m_outerCoarsening;
                    RealVect baseDst = (m_outerCoarseBox.smallEnd() +
                                        m_outerCoarseBuffer * IntVect::Unit)
                      * dxCoarse +
                      offset * m_dx;

                    IntVect srcFaceExtent =
                      srcFaceBox.bigEnd() - srcFaceBox.smallEnd();

                    int errcode = 0;
                    FORT_GETRECTDIRECTCOARSEMATRIX(m_sdcrseDirect[faceSrcID],
                                                   m_chargeWeights[faceSrcID]->dataPtr(),
                                                   &m_ncpts,
                                                   m_n2cOuter.dataPtr(),
                                                   &m_outerCoarseBuffer,
                                                   srcFaceExtent.dataPtr(),
                                                   m_dx.dataPtr(),
                                                   dxCoarse.dataPtr(),
                                                   baseSrc.dataPtr(),
                                                   baseDst.dataPtr(),
                                                   &m_verbose, &errcode);
                    if (errcode != 0)
                      {
                        cerr << "getrectdirectcoarsematrix returned error code "
                             << errcode
                             << " on face with ID " << faceSrcID
                             << endl;
                        MayDay::Error("returning");
                        return;
                      }

                    faceSrcID++;
                  }
              }
          }
        else if (m_parallel) // not direct:  using multipole expansions
          {
            m_sdcrsePar.define(m_srcFaces);
            for (DataIterator srcDit = m_srcFaces.dataIterator();
                 srcDit.ok(); ++srcDit)
              {
                // const Box& bxSrcFace = m_srcFaces.get(srcDit());
                Real*& sdcrseFace = m_sdcrsePar[srcDit()];

                int faceSrcID = m_faceSrcIDs[srcDit()];
                Side::LoHiSide sideSrc = (faceSrcID < SpaceDim) ?
                  Side::Lo : Side::Hi;
                int idirSrc = faceSrcID % SpaceDim;

                // m_sdcrsePar.resize(6);
                // hm, check that this isn't more than 2^32
                unsigned long szsdcrse =
                  sizeOuterMatrixPatches(idirSrc);

                pout() << "m_sdcrsePar using " << szsdcrse<< " scalars, "
                       << ((Real(szsdcrse) * sizeof(Real)) / 1048576.) << " MB"
                       << endl;
                // m_sdcrsePar = new Real[szsdcrse];
                // m_sdcrsePar = (Real *) malloc(szsdcrse * sizeof(Real));

                sdcrseFace =
                  (Real *) malloc(szsdcrse * sizeof(Real));

                int icp = (idirSrc + 1) * sign(sideSrc);

                int errcode = 0;
                FORT_GETRECTCOARSEMATRIX(sdcrseFace, &icp,
                                         m_n1.dataPtr(), m_n2.dataPtr(),
                                         &m_outerCoarsening,
                                         &m_patchSize,
                                         &m_numPatches[m_dir1[idirSrc]],
                                         &m_numPatches[m_dir2[idirSrc]],
                                         &m_multipoleOrder, m_dx.dataPtr(),
                                         &m_outerCoarseBuffer,
                                         m_cpind, m_cfac, m_cpow, &m_nterms,
                                         &m_ncpts, &m_verbose, &errcode);

                if (errcode != 0)
                  {
                    cerr << "getrectcoarsematrix returned error code "
                         << errcode
                         << " on face with ID " << faceSrcID
                         << endl;
                    MayDay::Error("returning");
                    return;
                  }
              }
          }
        else // ! m_parallel
          {
            // m_sdcrse.resize(6);
            // hm, check that this isn't more than 2^32
            unsigned long szsdcrseTotal = 0;
            for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
              {
                szsdcrseTotal +=
                  sizeOuterMatrixPatches(idirSrc);
              }

            pout() << "m_sdcrse using " << szsdcrseTotal << " scalars, "
                   << ((Real(szsdcrseTotal) * sizeof(Real)) / 1048576.) << " MB"
                   << endl;
            // m_sdcrse = new Real[szsdcrse];
            // m_sdcrse = (Real *) malloc(szsdcrse * sizeof(Real));

            Side::LoHiSide sideSrc = Side::Lo;
            for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
              {
                int icp = (idirSrc + 1) * sign(sideSrc);

                unsigned long szsdcrse =
                  sizeOuterMatrixPatches(idirSrc);
                m_sdcrse[idirSrc] =
                  (Real *) malloc(szsdcrse * sizeof(Real));
                int errcode = 0;
                FORT_GETRECTCOARSEMATRIX(m_sdcrse[idirSrc], &icp,
                                         m_n1.dataPtr(), m_n2.dataPtr(),
                                         &m_outerCoarsening,
                                         &m_patchSize,
                                         &m_numPatches[m_dir1[idirSrc]],
                                         &m_numPatches[m_dir2[idirSrc]],
                                         &m_multipoleOrder, m_dx.dataPtr(),
                                         &m_outerCoarseBuffer,
                                         m_cpind, m_cfac, m_cpow, &m_nterms,
                                         &m_ncpts, &m_verbose, &errcode);

                if (errcode != 0)
                  {
                    cerr << "getrectcoarsematrix returned error code "
                         << errcode
                         << " on face " << idirSrc
                         << endl;
                    MayDay::Error("returning");
                    return;
                  }
              }
          }
        } // if (m_outerCoarseBuffer > 0)
    }

  m_isDefined = true;
}


// ---------------------------------------------------------
// destructor
Multipoles::~Multipoles()
{
  clearMemory();
}


// ---------------------------------------------------------
bool
Multipoles::isDefined() const
{
  return m_isDefined;
}


// ---------------------------------------------------------
void
Multipoles::setDefaultValues()
{
  m_isDefined = false;
  m_verbose = 0;
}


// ---------------------------------------------------------
void
Multipoles::setVerbose(int a_verbose)
{
  m_verbose = a_verbose;
}


// ---------------------------------------------------------
// returns Multipoles to a basically undefined state
void
Multipoles::clearMemory()
{
  // if (m_levelOpPtr != NULL)
  // delete m_levelOpPtr;
  // m_levelOpPtr = NULL;
  if (m_isDefined)
    {
      clearCoeffs();
      if (m_direct)
        {
          for (int face = 0; face < 6; face++)
            delete m_chargeWeights[face];
          for (int ori = 0; ori < 6*6; ori++)
            free(m_srcdstDirect[ori]);
        }
      else
        {
          free(m_cpind);
          free(m_cpow);
          free(m_cfac);
          if (m_parallel)
            {
              for (DataIterator srcDit = m_srcFaces.dataIterator();
                   srcDit.ok(); ++srcDit)
                {
                  Vector<Real*>& srcdstVec = m_srcdstPar[srcDit()];
                  for (int i = 0; i < srcdstVec.size(); i++)
                    free(srcdstVec[i]);
                }

              if (m_getOuterCoarse && m_outerCoarseBuffer > 0)
                for (DataIterator srcDit = m_srcFaces.dataIterator();
                     srcDit.ok(); ++srcDit)
                  {
                    Real*& sdcrseFace = m_sdcrsePar[srcDit()];
                    free(sdcrseFace);
                  }
            }
          else
            {
              for (int idir = 0; idir < SpaceDim; idir++)
                for (int i = 0; i < m_srcdst[idir].size(); i++)
                  free(m_srcdst[idir][i]);
              if (m_getOuterCoarse && m_outerCoarseBuffer > 0)
                for (int idir = 0; idir < SpaceDim; idir++)
                  free(m_sdcrse[idir]);
            }
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              delete m_startPatch[idir];
              delete m_whichWeights[idir];
              Vector<FArrayBox*> weightsFace = m_weights[idir];
              for (int i = 0; i < weightsFace.size(); i++)
                delete weightsFace[i];
            }
        }
    }
}


// ---------------------------------------------------------
void
Multipoles::setCoeffs(const Vector<FArrayBox*>&   a_faceCharge)
{
  CH_assert(isDefined());
  if (m_parallel)
    // for (DataIterator srcDit = m_srcFaces.dataIterator();
    for (DataIterator srcDit = m_coeffsLayout.dataIterator();
         srcDit.ok(); ++srcDit)
      {
        int faceSrcID = m_faceSrcIDs[srcDit()];
        //      Side::LoHiSide sideSrc = (faceSrcID < SpaceDim) ?
        //        Side::Lo : Side::Hi;
        int idirSrc = faceSrcID % SpaceDim;

        FArrayBox& thisFaceCharge = *(a_faceCharge[faceSrcID]);
        FArrayBox& coeffsFace = m_coeffsPar[srcDit()];

        BaseFab<IntVect>& startPatchDir = *m_startPatch[idirSrc];
        IntVect startFace = a_faceCharge[faceSrcID]->smallEnd();
        Vector<FArrayBox*> m_weightsDir = m_weights[idirSrc];
        if (m_weightsDir.size() == 1)
          {
            FArrayBox* weightsPatchPtr = m_weightsDir[0];
            for (BoxIterator bit(m_bxPatchesFace[idirSrc]); bit.ok(); ++bit)
              {
                IntVect patch = bit();
                // You cannot shift a_faceCharge (const),
                // so you must shift weights
                // to match the appropriate patch of a_faceCharge.
                IntVect shiftNodes = startFace + patch * m_patchOffset[idirSrc];
                Box nodesOnPatch = weightsPatchPtr->box() + shiftNodes;
                weightsPatchPtr->shift(shiftNodes);
                CH_assert(coeffsFace.box().contains(patch));
                FORT_MULTIPOLECOEFFS(CHF_FRA(coeffsFace),
                                     CHF_CONST_INTVECT(patch),
                                     CHF_BOX(nodesOnPatch),
                                     CHF_CONST_FRA1(thisFaceCharge, 0),
                                     CHF_CONST_FRA((*weightsPatchPtr)));
                weightsPatchPtr->shift(-shiftNodes);
              }
          }
        else
          {
            const BaseFab<int>& whichWeights(*(m_whichWeights[idirSrc]));
            for (BoxIterator bit(m_bxPatchesFace[idirSrc]); bit.ok(); ++bit)
              {
                IntVect patch = bit();
                // You cannot shift a_faceCharge (const),
                // so you must shift weights
                // to match the appropriate patch of a_faceCharge.
                IntVect shiftNodes =
                  startFace + startPatchDir(patch, 0);

                // whichWeights(patch, 0) is one of 0, 1, 2, 3
                FArrayBox* weightsPatchPtr =
                  m_weightsDir[whichWeights(patch, 0)];
                Box nodesOnPatch = weightsPatchPtr->box() + shiftNodes;

                weightsPatchPtr->shift(shiftNodes);
                CH_assert(coeffsFace.box().contains(patch));
                FORT_MULTIPOLECOEFFS(CHF_FRA(coeffsFace),
                                     CHF_CONST_INTVECT(patch),
                                     CHF_BOX(nodesOnPatch),
                                     CHF_CONST_FRA1(thisFaceCharge, 0),
                                     CHF_CONST_FRA((*weightsPatchPtr)));
                weightsPatchPtr->shift(-shiftNodes);
              }
          }
        m_gotCoeffs[faceSrcID] = true;
      }
  else // ! m_parallel
    {
      int faceSrcID = 0;
      for (SideIterator sit; sit.ok(); sit.next())
        for (int idir = 0; idir < SpaceDim; idir++)
          {
            setCoeffs(*a_faceCharge[faceSrcID], idir, sit());
            faceSrcID++;
          }
    }
}


// ---------------------------------------------------------
void
Multipoles::setCoeffs(const FArrayBox&        a_faceCharge,
                      int                     a_idirSrc,
                      const Side::LoHiSide&   a_side)
{
  CH_assert(isDefined());
  int faceSrcID = (a_side == Side::Lo) ? a_idirSrc : (a_idirSrc + SpaceDim);
  if (m_direct)
    {
      FArrayBox* faceChargePtr =
        new FArrayBox(a_faceCharge.box(), 1);
      faceChargePtr->copy(a_faceCharge);
      // *faceChargePtr *= *m_chargeWeights[faceSrcID];
      m_charges[faceSrcID] = faceChargePtr;
    }
  else if (! m_parallel)
    {
      // Recall m_bxPatchesFace[a_idirSrc] is CELL-centered:
      // 0:m_numPatches-1 in the two parallel directions,
      // 0 in normal direction.
      FArrayBox* coeffsFacePtr =
        new FArrayBox(m_bxPatchesFace[a_idirSrc], m_numCoeffs);
      m_coeffs[faceSrcID] = coeffsFacePtr;

      IntVect startFace = a_faceCharge.smallEnd();
      BaseFab<IntVect>& startPatchDir = *m_startPatch[a_idirSrc];
      Vector<FArrayBox*> m_weightsDir = m_weights[a_idirSrc];
      if (m_weightsDir.size() == 1)
        {
          FArrayBox* weightsPatchPtr = m_weightsDir[0];
          for (BoxIterator bit(m_bxPatchesFace[a_idirSrc]); bit.ok(); ++bit)
            {
              IntVect patch = bit();
              // You cannot shift a_faceCharge (const),
              // so you must shift weights
              // to match the appropriate patch of a_faceCharge.
              IntVect shiftNodes =
                startFace + startPatchDir(patch, 0);
              Box nodesOnPatch = weightsPatchPtr->box() + shiftNodes;

              weightsPatchPtr->shift(shiftNodes);
              CH_assert(coeffsFacePtr->box().contains(patch));
              FORT_MULTIPOLECOEFFS(CHF_FRA((*coeffsFacePtr)),
                                   CHF_CONST_INTVECT(patch),
                                   CHF_BOX(nodesOnPatch),
                                   CHF_CONST_FRA1(a_faceCharge, 0),
                                   CHF_CONST_FRA((*weightsPatchPtr)));
              weightsPatchPtr->shift(-shiftNodes);
            }
        }
      else
        {
          const BaseFab<int>& whichWeights(*(m_whichWeights[a_idirSrc]));
          for (BoxIterator bit(m_bxPatchesFace[a_idirSrc]); bit.ok(); ++bit)
            {
              IntVect patch = bit();
              // You cannot shift a_faceCharge (const),
              // so you must shift weights
              // to match the appropriate patch of a_faceCharge.
              IntVect shiftNodes =
                startFace + startPatchDir(patch, 0);
              // whichWeights(patch, 0) is one of 0, 1, 2, 3
              FArrayBox* weightsPatchPtr =
                m_weightsDir[whichWeights(patch, 0)];
              Box nodesOnPatch = weightsPatchPtr->box() + shiftNodes;

              weightsPatchPtr->shift(shiftNodes);
              CH_assert(coeffsFacePtr->box().contains(patch));
              FORT_MULTIPOLECOEFFS(CHF_FRA((*coeffsFacePtr)),
                                   CHF_CONST_INTVECT(patch),
                                   CHF_BOX(nodesOnPatch),
                                   CHF_CONST_FRA1(a_faceCharge, 0),
                                   CHF_CONST_FRA((*weightsPatchPtr)));
              weightsPatchPtr->shift(-shiftNodes);
            }
        }
    }
  else
    {
      cerr << "called Multipoles::setCoeffs in parallel; not done yet"
           << endl;
      MayDay::Error("returning");
      return;
    }
  m_gotCoeffs[faceSrcID] = true;
}


// ---------------------------------------------------------
void
Multipoles::clearCoeffs()
{
  CH_assert(isDefined());
  for (int faceSrcID = 0; faceSrcID < 2*SpaceDim; faceSrcID++)
    if (m_gotCoeffs[faceSrcID])
      {
        if (m_direct)
          {
            delete m_charges[faceSrcID];
          }
        else if (! m_parallel)
          {
            delete m_coeffs[faceSrcID];
          }
        m_gotCoeffs[faceSrcID] = false;
      }
}


// ---------------------------------------------------------
void
Multipoles::eval(FArrayBox&   a_dst)
{
  CH_assert(isDefined());
  const Interval intvl0(0, 0);
  // Set a_dst to zero on faces.
  int faceDstID = 0;
  for (SideIterator sitDst; sitDst.ok(); sitDst.next())
    for (int idirDst = 0; idirDst < SpaceDim; idirDst++)
      {
        Box dstBox = m_bxDstFace[faceDstID];
        a_dst.setVal(0., dstBox, 0);
        faceDstID++;
      }
  // Add up interpolations onto each face.
  faceDstID = 0;
  for (SideIterator sitDst; sitDst.ok(); sitDst.next())
    {
      Side::LoHiSide sideDst = sitDst();
      for (int idirDst = 0; idirDst < SpaceDim; idirDst++)
        {
          FArrayBox dstFaceCoarseAll(m_bxDstCoarseFace[faceDstID], 1);
          // dstFaceCoarseAll will contain the sum of potentials due
          // to all 2*SpaceDim source faces.
          dstFaceCoarseAll.setVal(0.);
          if (! m_parallel)
            {
              int faceSrcID = 0;
              for (SideIterator sitSrc; sitSrc.ok(); sitSrc.next())
                {
                  Side::LoHiSide sideSrc = sitSrc();
                  for (int idirSrc = 0; idirSrc < SpaceDim; idirSrc++)
                    {
                      int ori = (2*SpaceDim)*faceDstID + faceSrcID;
                      if (m_direct)
                        {
                          const FArrayBox& srcFaceCharges = *m_charges[faceSrcID];
                          Box srcFaceBox = srcFaceCharges.box();

                          IntVect lengthSrcFace =
                            srcFaceBox.bigEnd() - srcFaceBox.smallEnd();

                          int errcode = 0;
                          FORT_ADDEVAL1DIRECTFACECOARSE(dstFaceCoarseAll.dataPtr(),
                                                        srcFaceCharges.dataPtr(),
                                                        lengthSrcFace.dataPtr(),
                                                        &m_n2c[m_dir1[idirDst]],
                                                        &m_n2c[m_dir2[idirDst]],
                                                        &m_interpBorder,
                                                        m_srcdstDirect[ori],
                                                        &m_verbose, &errcode);

                          if (errcode != 0)
                            {
                              cerr << "addeval1directfacecoarse returned error code "
                                   << errcode << endl;
                              MayDay::Error("returning");
                              return;
                            }
                        }
                      else // ! m_direct
                        {
                          DestFace dstCode;
                          int srcFlip, dstFlip;
                          getCodeFlips(dstCode, srcFlip, dstFlip,
                                       idirSrc, sideSrc, idirDst, sideDst);

                          // int ddirSrc = sign(sideSrc) * (idirSrc + 1);
                          // int ddirDst = sign(sideDst) * (idirDst + 1);

                          int errcode = 0;
                          FORT_ADDEVALSYM1FACECOARSE(dstFaceCoarseAll.dataPtr(),
                                                     m_coeffs[faceSrcID]->dataPtr(),
                                                     &dstFlip, &srcFlip,
                                                     &m_numPatches[m_dir1[idirSrc]],
                                                     &m_numPatches[m_dir2[idirSrc]],
                                                     &m_n2c[m_dir1[idirDst]],
                                                     &m_n2c[m_dir2[idirDst]],
                                                     &m_multipoleOrder,
                                                     &m_numCoeffs, &m_interpBorder,
                                                     m_srcdst[idirSrc][dstCode],
                                                     &m_verbose, &errcode);
                          if (errcode != 0)
                            {
                              cerr << "addevalsym1facecoarse returned error code "
                                   << errcode << endl;
                              MayDay::Error("returning");
                              return;
                            }
                        }
                      faceSrcID++;
                    }
                }
            }
          else // m_parallel
            { // not direct:  evaluate multipole expansions
              // FArrayBox dstFaceCoarseSum(m_bxDstCoarseFace[faceDstID], 1);
              // There is a separate dstFaceCoarseSum for each processor
              // that holds any source faces:
              // these are in dstFaceCoarseSumLayoutData.
              // Later we'll add them all up.

              // dstFaceCoarseLayout consists of a copy of the same box
              // m_bxDstCoarseFace[faceDstID] on
              // every processor on which the program is running.
              int nproc = numProc();
              Vector<int> procs(nproc);
              for (int iproc = 0; iproc < nproc; iproc++) procs[iproc] = iproc;

              Vector<Box> dstFaceCoarseBoxes(nproc);
              dstFaceCoarseBoxes.assign(m_bxDstCoarseFace[faceDstID]);
              BoxLayout dstFaceCoarseLayout(dstFaceCoarseBoxes, procs);

              // dstFaceCoarseSumLayout consists of a copy of the same box
              // m_bxDstCoarseFace[faceDstID] on
              // every processor that stores a source face.
              // (There are at most 2*SpaceDim such processors.)
              int ndistinct = m_faceSrcProcDistinct.size();
              Vector<Box> dstFaceCoarseSumBoxes(ndistinct);
              Vector<int> procsdistinct(ndistinct);
              dstFaceCoarseSumBoxes.assign(m_bxDstCoarseFace[faceDstID]);
              for (int iproc = 0; iproc < ndistinct; iproc++)
                procsdistinct[iproc] = m_faceSrcProcDistinct[iproc];
              BoxLayout dstFaceCoarseSumLayout(dstFaceCoarseSumBoxes,
                                               procsdistinct);

              // Each processor holds at most one FArrayBox
              // dstFaceCoarseSumLayoutData[dit()], which will contain the
              // sum of the potentials due to all source faces stored
              // by that processor.
              BoxLayoutData<FArrayBox>
                dstFaceCoarseSumLayoutData(dstFaceCoarseSumLayout, 1);

              // dstDomain is used in generalCopyTo().
              ProblemDomain dstDomain(m_bxDstCoarseFace[faceDstID]);

              // loop over the processors that store source faces
              for (DataIterator dstDit = dstFaceCoarseSumLayout.dataIterator();
                   dstDit.ok(); ++dstDit)
                {
                  // dstFaceCoarseSum is sum of potentials due to
                  // source faces stored on this processor.
                  FArrayBox& dstFaceCoarseSum =
                    dstFaceCoarseSumLayoutData[dstDit()];
                  dstFaceCoarseSum.setVal(0.);
                  for (DataIterator srcDit = m_srcFaces.dataIterator();
                       srcDit.ok(); ++srcDit)
                    {
                      // FArrayBox dstFaceCoarse(m_bxDstCoarseFace[faceDstID], 1);
                      // const Box& bxSrcFace = m_srcFaces.get(srcDit());
                      const Vector<Real*>& srcdstVec = m_srcdstPar[srcDit()];

                      int faceSrcID = m_faceSrcIDs[srcDit()];
                      Side::LoHiSide sideSrc = (faceSrcID < SpaceDim) ?
                        Side::Lo : Side::Hi;
                      int idirSrc = faceSrcID % SpaceDim;

                      // int isideSrc = sign(sideSrc);
                      // int ddirSrc = isideSrc * (idirSrc + 1);

                      DestFace dstCode;
                      int srcFlip, dstFlip;
                      getCodeFlips(dstCode, srcFlip, dstFlip,
                                   idirSrc, sideSrc, idirDst, sideDst);
                      int errcode = 0;
                      FORT_ADDEVALSYM1FACECOARSE(dstFaceCoarseSum.dataPtr(),
                                                 m_coeffsPar[srcDit()].dataPtr(),
                                                 &dstFlip, &srcFlip,
                                                 &m_numPatches[m_dir1[idirSrc]],
                                                 &m_numPatches[m_dir2[idirSrc]],
                                                 &m_n2c[m_dir1[idirDst]],
                                                 &m_n2c[m_dir2[idirDst]],
                                                 &m_multipoleOrder,
                                                 &m_numCoeffs, &m_interpBorder,
                                                 srcdstVec[dstCode],
                                                 &m_verbose, &errcode);
                      if (errcode != 0)
                        {
                          cerr << "addevalsym1facecoarse returned error code "
                               << errcode << endl;
                          MayDay::Error("returning");
                          return;
                        }
                    } // end iteration over m_srcFaces
                } // end iteration over dest procs
              // pout() << "dstFaceCoarseLayout == " << dstFaceCoarseLayout
              // << endl;
              // pout() << "dstFaceCoarseSumLayout == " << dstFaceCoarseSumLayout
              // << endl;
              LayoutData< Vector<RefCountedPtr<FArrayBox> > >
                dstFaceCoarseInt;
              dstFaceCoarseSumLayoutData.generalCopyTo(// dest layout
                                                       dstFaceCoarseLayout,
                                                       // dest
                                                       dstFaceCoarseInt,
                                                       // copy component 0
                                                       intvl0,
                                                       // domain of source, dest
                                                       dstDomain);

              // Do the communication BEFORE interpolation.
              // This reduces computation (not as many interpolations)
              // as well as communication (send coarse data only).
              // Recall that on each processor, there's only one element in
              // dstFaceCoarseLayout.
              // FArrayBox dstFaceCoarseAll(m_bxDstCoarseFace[faceDstID], 1);
              for (DataIterator dstDit = dstFaceCoarseLayout.dataIterator();
                   dstDit.ok(); ++dstDit)
                {
                  plusReduce(dstFaceCoarseAll, dstFaceCoarseInt[dstDit()]);
                }

              // FArrayBox dstFaceAll(dstBox, 1);
            }
          Box dstBox = m_bxDstFace[faceDstID];
          FArrayBox dstFace(dstBox, 1);
          // dstFace.setVal(0.);
          int shiftFine = m_bxDstFace[faceDstID].smallEnd(idirDst);
          dstFace.shift(idirDst, -shiftFine);
          // Box of dstFaceCoarse has small end at (0, 0, 0).
          m_interp[idirDst].interpolate(dstFace, dstFaceCoarseAll);
          dstFace.shift(idirDst, +shiftFine);
          // Multiply dstFace by 1/3 on corners and by 1/2 on edges:
          // this is averaging on intersections of destination faces.
          // int m_dir1[idirDst] = (idirDst == 0) ? 1 : 0;
          // int m_dir2[idirDst] = (idirDst == 2) ? 1 : 2;
          IntVect e1 = BASISV(m_dir1[idirDst]);
          IntVect e2 = BASISV(m_dir2[idirDst]);
          // multiply corners by 1/3
          IntVect ivLoLo = dstBox.smallEnd();
          IntVect ivHiHi = dstBox.bigEnd();
          IntVect ivLoHi(ivLoLo);
          ivLoHi.setVal(m_dir2[idirDst], ivHiHi[m_dir2[idirDst]]);
          IntVect ivHiLo(ivLoLo);
          ivHiLo.setVal(m_dir1[idirDst], ivHiHi[m_dir1[idirDst]]);
          dstFace(ivLoLo, 0) /= 3.;
          dstFace(ivLoHi, 0) /= 3.;
          dstFace(ivHiLo, 0) /= 3.;
          dstFace(ivHiHi, 0) /= 3.;
          // multiply edges (except corners) by 1/2
          Box bxEdge1Lo(ivLoLo + e2, ivLoHi - e2, IndexType::TheNodeType());
          Box bxEdge1Hi(ivHiLo + e2, ivHiHi - e2, IndexType::TheNodeType());
          Box bxEdge2Lo(ivLoLo + e1, ivHiLo - e1, IndexType::TheNodeType());
          Box bxEdge2Hi(ivLoHi + e1, ivHiHi - e1, IndexType::TheNodeType());
          dstFace.divide(2., bxEdge1Lo);
          dstFace.divide(2., bxEdge1Hi);
          dstFace.divide(2., bxEdge2Lo);
          dstFace.divide(2., bxEdge2Hi);
          // this WAS a_dst.plus()
          a_dst.minus(dstFace, dstBox, 0, 0);
          faceDstID++;
        }
    }
}


// ---------------------------------------------------------
void
Multipoles::eval(FArrayBox&              a_dstFace,
                 int                     a_idirSrc,
                 const Side::LoHiSide&   a_sideSrc,
                 int                     a_idirDst,
                 const Side::LoHiSide&   a_sideDst)
{
  CH_assert(isDefined());
  int faceSrcID = (a_sideSrc == Side::Lo) ? a_idirSrc : (a_idirSrc + SpaceDim);
  CH_assert(m_gotCoeffs[faceSrcID]);
  int faceDstID = (a_sideDst == Side::Lo) ? a_idirDst : (a_idirDst + SpaceDim);
  // FArrayBox dstFaceCoarse(m_bxDstFace[faceDstID], 1);
  FArrayBox dstFaceCoarse(m_bxDstCoarseFace[faceDstID], 1);
  dstFaceCoarse.setVal(0.);
  int ori = (2*SpaceDim)*faceDstID + faceSrcID;
  // cout << "from " << ((a_sideSrc == Side::Lo) ? -1 : +1) * (a_idirSrc+1)
  // << " to " << ((a_sideDst == Side::Lo) ? -1 : +1) * (a_idirDst+1)
  // << " ori = " << ori << endl;
  if (m_direct)
    {
      const FArrayBox& srcFaceCharges = *m_charges[faceSrcID];
      Box srcFaceBox = srcFaceCharges.box();

      IntVect lengthSrcFace =
        srcFaceBox.bigEnd() - srcFaceBox.smallEnd();

      int errcode = 0;
      FORT_ADDEVAL1DIRECTFACECOARSE(dstFaceCoarse.dataPtr(),
                                    srcFaceCharges.dataPtr(),
                                    lengthSrcFace.dataPtr(),
                                    &m_n2c[m_dir1[a_idirDst]],
                                    &m_n2c[m_dir2[a_idirDst]],
                                    &m_interpBorder,
                                    m_srcdstDirect[ori], &m_verbose, &errcode);

      if (errcode != 0)
        {
          cerr << "addeval1directfacecoarse returned error code "
               << errcode << endl;
          MayDay::Error("returning");
          return;
        }
    }
  else if (! m_parallel)
    {
      DestFace dstCode;
      int srcFlip, dstFlip;
      getCodeFlips(dstCode, srcFlip, dstFlip,
                   a_idirSrc, a_sideSrc, a_idirDst, a_sideDst);

      int errcode = 0;
      FORT_ADDEVALSYM1FACECOARSE(dstFaceCoarse.dataPtr(),
                                 m_coeffs[faceSrcID]->dataPtr(),
                                 &dstFlip, &srcFlip,
                                 &m_numPatches[m_dir1[a_idirSrc]],
                                 &m_numPatches[m_dir2[a_idirSrc]],
                                 &m_n2c[m_dir1[a_idirDst]],
                                 &m_n2c[m_dir2[a_idirDst]],
                                 &m_multipoleOrder,
                                 &m_numCoeffs, &m_interpBorder,
                                 m_srcdst[a_idirSrc][dstCode],
                                 &m_verbose, &errcode);
      if (errcode != 0)
        {
          cerr << "addevalsym1facecoarse returned error code "
               << errcode << endl;
          MayDay::Error("returning");
          return;
        }
    }
  else // m_parallel
    {
      cerr << "called Multipoles::eval in parallel; not done yet"
           << endl;
      MayDay::Error("returning");
      return;
    }
  int shiftFine = m_bxDstFace[faceDstID].smallEnd(a_idirDst);
  a_dstFace.shift(a_idirDst, -shiftFine);
  // Box of dstFaceCoarse has small end at (0, 0, 0).
  m_interp[a_idirDst].interpolate(a_dstFace, dstFaceCoarse);
  a_dstFace.shift(a_idirDst, +shiftFine);
}


// ---------------------------------------------------------
void
Multipoles::evalOuterCoarse(FArrayBox&   a_dstOuter)
{
  if (m_outerCoarseBuffer > 0) // otherwise nothing to do
    {
      CH_assert(a_dstOuter.box() == m_outerCoarseBox);
      a_dstOuter.setVal(0.);
      for (SideIterator sit; sit.ok(); sit.next())
        {
          Side::LoHiSide side = sit();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              FArrayBox dstFromFace(a_dstOuter.box(), 1);
              evalOuterCoarse(dstFromFace, idir, side);
              a_dstOuter.plus(dstFromFace);
            }
        }
    }
}


// ---------------------------------------------------------
void
Multipoles::evalOuterCoarse(FArrayBox&              a_dstOuter,
                            int                     a_idirSrc,
                            const Side::LoHiSide&   a_sideSrc)
{
  if (m_outerCoarseBuffer > 0)
    {
      int faceSrcID = (a_sideSrc == Side::Lo) ?
        a_idirSrc : (a_idirSrc + SpaceDim);

      if (m_direct)
        {
          const FArrayBox& faceCharges = *m_charges[faceSrcID];
          const Box& srcFaceBox = faceCharges.box();
          RealVect baseSrc = srcFaceBox.smallEnd() * m_dx;

          RealVect dxCoarse = m_outerCoarsening * m_dx;
          IntVect offset; // same as offset in getrectcoarsematrix
          for (int idirOther = 0; idirOther < SpaceDim; idirOther++)
            offset[idirOther] = ((m_n2[idirOther] / 2) % m_outerCoarsening) /
              m_outerCoarsening;
          RealVect baseDst = (m_outerCoarseBox.smallEnd() +
                              m_outerCoarseBuffer * IntVect::Unit) * dxCoarse +
            offset * m_dx;

          IntVect srcFaceExtent = srcFaceBox.bigEnd() - srcFaceBox.smallEnd();

          int errcode;
          FORT_EVALDIRECTRECTOUTERCOARSE(a_dstOuter.dataPtr(),
                                         faceCharges.dataPtr(),
                                         &m_ncpts,
                                         m_n2cOuter.dataPtr(),
                                         &m_outerCoarseBuffer,
                                         srcFaceExtent.dataPtr(),
                                         m_sdcrseDirect[faceSrcID],
                                         &m_verbose, &errcode);
          if (errcode != 0)
            {
              cerr << "evaldirectrectoutercoarse returned error code "
                   << errcode << endl;
              MayDay::Error("returning");
              return;
            }
        }
      else if (m_parallel)
        {
          // We have to go through the whole ugly procedure again of
          // doing a reduce, etc.
          // FIX this -- it's NOT correct now.
          // You have to be a little more clever than with
          // dstFaceCoarseSum above, because you just want a layer
          // -- don't send around the whole FArrayBox.
          // Probably copy to an array of size m_ncpts, and send that.
          // Perhaps that might also be helpful in the serial case?
          // Or, split up the FArrayBox into rectangular pieces:
          // need 6 of them.
          for (DataIterator srcDit = m_srcFaces.dataIterator();
               srcDit.ok(); ++srcDit)
            {
              int errcode;
              FORT_EVALRECTOUTERCOARSE(a_dstOuter.dataPtr(),
                                       m_coeffsPar[srcDit()].dataPtr(),
                                       &m_ncpts,
                                       &m_numPatches[m_dir1[a_idirSrc]],
                                       &m_numPatches[m_dir2[a_idirSrc]],
                                       m_n2cOuter.dataPtr(),
                                       &m_outerCoarseBuffer,
                                       &m_numCoeffs,
                                       m_sdcrsePar[srcDit()],
                                       &m_verbose, &errcode);
              if (errcode != 0)
                {
                  cerr << "evalrectoutercoarse returned error code "
                       << errcode << endl;
                  MayDay::Error("returning");
                  return;
                }
            }
        }
      else // ! m_parallel
        {
          // If source face is on high side, then set dstFlip to
          // 1 or 2 or 3, depending on the dimension.
          int dstFlip = (a_sideSrc == Side::Lo) ? 0 : (a_idirSrc + 1);
          int errcode;
          FORT_EVALSYMRECTOUTERCOARSE(a_dstOuter.dataPtr(),
                                      m_coeffs[faceSrcID]->dataPtr(),
                                      &dstFlip,
                                      &m_ncpts,
                                      &m_numPatches[m_dir1[a_idirSrc]],
                                      &m_numPatches[m_dir2[a_idirSrc]],
                                      m_n2cOuter.dataPtr(),
                                      &m_outerCoarseBuffer,
                                      &m_numCoeffs,
                                      m_sdcrse[a_idirSrc],
                                      &m_verbose, &errcode);

          if (errcode != 0)
            {
              cerr << "evalsymrectoutercoarse returned error code "
                   << errcode << endl;
              MayDay::Error("returning");
              return;
            }
        }
    }
}

#include "NamespaceFooter.H"
