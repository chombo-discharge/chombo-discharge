#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// #include <cstdio>

#include "FourthOrderMappedCFInterp.H"
#include "BoxIterator.H"
#include "GradientF_F.H"
#include "MappedCFInterpF_F.H"

#include "NamespaceHeader.H"

//////////////////////////////////////////////////////////////////////////////
// Constructor - set up some defaults
FourthOrderMappedCFInterp::FourthOrderMappedCFInterp()
{
  m_defined = false;
  m_haveJ = false;
  m_haveCoeffs = false;
}

//////////////////////////////////////////////////////////////////////////////
// Destructor - free up storage
FourthOrderMappedCFInterp::~FourthOrderMappedCFInterp()
{
  clearMemory();
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderMappedCFInterp::clearMemory()
{
  if (m_defined)
    {
      for (BoxIterator bit(m_stencilBox); bit.ok(); ++bit)
        {
          IntVect offset = bit();
          delete m_stencils(offset, 0);
        }
      freeJ();
    }
}


//////////////////////////////////////////////////////////////////////////////
// Define the object so that time stepping can begin
void FourthOrderMappedCFInterp::define(/// layout at finer level
                                       const DisjointBoxLayout&  a_layout,
                                       /// number of variables
                                       int                       a_numStates,
                                       /// refinement ratio between finer level and coarser level
                                       int                       a_refineCoarse,
                                       /// problem domain on finer level
                                       const ProblemDomain&      a_domain,
                                       /// width of coarse ghost cell layer with fine cells to be filled in
                                       int                       a_coarseGhostsFill,
                                       /// dimensions that are fixed, not interpolated
                                       Interval                  a_fixedDims)
{
  CH_TIME("FourthOrderMappedCFInterp::define");
  // Free up existing storage, if already defined.
  clearMemory();

  // Cache data
  m_domain = a_domain;
  m_refineCoarse = a_refineCoarse;
  m_numStates = a_numStates;
  m_coarseGhostsFill = a_coarseGhostsFill;
  m_fixedDims = a_fixedDims;
  m_nInterpDims = SpaceDim - m_fixedDims.size();
  // m_ghostVect = a_ghostVect;

  m_interpUnit = IntVect::Unit;
  m_refineVect = m_refineCoarse * IntVect::Unit;
  for (int dirf = m_fixedDims.begin(); dirf <= m_fixedDims.end(); dirf++)
    {
      m_interpUnit[dirf] = 0;
      m_refineVect[dirf] = 1;
    }

  /*
    For coarse level.
   */

  m_coarseDomain = coarsen(m_domain, m_refineVect);

  m_layout = a_layout;
  // petermc, 19 Dec 2008:  prevents crash in case of calling define again
  m_layoutCoarsened = DisjointBoxLayout();
  coarsen(m_layoutCoarsened, m_layout, m_refineVect);

  /// maximum distance of stencil from domain boundary
  m_maxStencilDist = 2;

  // For 4th-order interpolation, take degree 3.
  m_degree = 3;

  m_stencilBox = Box(-m_maxStencilDist*m_interpUnit,
                     m_maxStencilDist*m_interpUnit);
  m_stencils.define(m_stencilBox, 1);
  for (BoxIterator bit(m_stencilBox); bit.ok(); ++bit)
    {
      IntVect offset = bit();
      m_stencils(offset, 0) =
        new FourthOrderMappedCFInterpStencil(offset, m_degree, m_fixedDims);
    }

  // Now find m_stencilHere:
  // For each coarse cell whose fine cells will be filled in,
  // figure out which stencil to use.
  m_coarseGhostVect = m_coarseGhostsFill * m_interpUnit;
  // LevelData< BaseFab<IntVect> > m_stencilHere;
  m_stencilHere.define(m_layoutCoarsened, 1, m_coarseGhostVect);

  DataIterator dit = m_layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bxCoarse = m_layoutCoarsened[dit];
      // petermc, 9 Jan 2009:  grow by m_coarseGhostVect for filling in ghosts
      Box coarseGhostedBox = grow(bxCoarse, m_coarseGhostVect);

      BaseFab<IntVect>& stencilHereFab = m_stencilHere[dit];
      setStencil(stencilHereFab, coarseGhostedBox);
    }

  // layers of ghost cells needed to interpolate to fine data,
  // which have m_coarseGhostsFill layers of ghost cells
  m_ghostVect = (m_coarseGhostsFill + m_maxStencilDist) * m_interpUnit;
  // used only in interpToFine function
  m_coarseData.define(m_layoutCoarsened, m_numStates, m_ghostVect);

  /*
    For fine level.
   */

  Box degreeBox(IntVect::Zero, m_degree * m_interpUnit);
  m_numNonzeroPowers = 0;
  for (BoxIterator bit(degreeBox); bit.ok(); ++bit)
    {
      IntVect iv = bit();
      int sumPowers = iv.sum();
      if ( (sumPowers > 0) && (sumPowers <= m_degree) )
        { // IntVect iv passes the test, so count it.
          m_numNonzeroPowers++;
        }
    }

  // Now find <(z - z0)^p>k for each fine cell k.
  m_fineBox = Box(IntVect::Zero, m_refineVect - IntVect::Unit);
  Box fine1Box = grow(m_fineBox, m_interpUnit);
  m_fineAvgPowers.define(fine1Box, m_numNonzeroPowers);
  int numFine = fine1Box.numPts();
  FORT_FINEAVGPOWERS(CHF_FRA(m_fineAvgPowers),
                     CHF_CONST_INT(m_refineCoarse),
                     CHF_CONST_INT(numFine),
                     CHF_CONST_INT(m_numNonzeroPowers),
                     CHF_CONST_INT(m_degree),
                     CHF_BOX(fine1Box),
                     CHF_BOX(degreeBox));

  // undivided differences
  m_gradFineAvgPowers.define(m_fineBox, m_nInterpDims * m_numNonzeroPowers);
  getGradient(m_gradFineAvgPowers, m_fineAvgPowers, m_fineBox);

  // Everything is defined now, except J and U, JU.
  m_defined = true;
  m_haveJ = false;
  m_haveCoeffs = false;
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderMappedCFInterp::setStencil(BaseFab<IntVect>&  a_stencil,
                                           const Box&         a_bxCoarse)
{
  CH_TIME("FourthOrderMappedCFInterp::setStencil");
  // This is copied from FourthOrderPatchInterp::setStencil().
  const Box& coarseDomainBox = m_coarseDomain.domainBox();
  const IntVect& coarseDomainLo = coarseDomainBox.smallEnd();
  const IntVect& coarseDomainHi = coarseDomainBox.bigEnd();
  for (BoxIterator bit(a_bxCoarse); bit.ok(); ++bit)
    {
      IntVect ivc = bit();
      // Set IntVect stencilHereFab(ivc, 0).
      // Find distance to coarseDomain boundary
      IntVect dist = IntVect::Zero;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (! m_fixedDims.contains(idir))
            {
              if (! m_coarseDomain.isPeriodic(idir))
                {
                  int offLo = coarseDomainLo[idir] - ivc[idir] - 1;
                  int offHi = coarseDomainHi[idir] - ivc[idir] + 1;
                  if (offLo < 0 && offHi >0) // condition means ivc in coarseDomain
                    {
                      if ((offLo >= -m_maxStencilDist) &&
                          (offHi <= m_maxStencilDist))
                        { // both of these:  very narrow domain, you are in trouble
                          MayDay::Error("FourthOrderMappedCFInterp::setStencil bad boxes");
                        }
                      if (offLo >= -m_maxStencilDist) // -1 or -2
                        dist[idir] = offLo;
                      if (offHi <= m_maxStencilDist) // 1 or 2
                        dist[idir] = offHi;
                      // Otherwise, dist[idir] = 0.
                    }
                }
            }
        }
      a_stencil(ivc, 0) = dist;
    }
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderMappedCFInterp::getGradient(FArrayBox&        a_gradFab,
                                            const FArrayBox&  a_srcFab,
                                            const Box&        a_bx)
{
  CH_TIME("FourthOrderMappedCFInterp::getGradient");
  int nSrcComp = a_srcFab.nComp();
  int gradComp = 0;
  Real bogusdx = 1.;
  for (int srcComp = 0; srcComp < nSrcComp; srcComp++)
    {
      for (int idir = 0; idir < SpaceDim; idir++)
        if ( !m_fixedDims.contains(idir) )
          {
            FORT_GRADCC(CHF_FRA1(a_gradFab, gradComp),
                        CHF_CONST_FRA1(a_srcFab, srcComp),
                        CHF_BOX(a_bx),
                        CHF_CONST_REAL(bogusdx),
                        CHF_INT(idir));
            gradComp++;
          }
    }
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderMappedCFInterp::freeJ()
{
  CH_TIME("FourthOrderMappedCFInterp::freeJ");
  if (m_haveJ)
    {
      DataIterator dit = m_layout.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& bxCoarse = m_layoutCoarsened[dit];
          Box coarseGhostedBox = grow(bxCoarse, m_coarseGhostVect);
          BaseFab< Vector<Real>* >& A1Fab = *(m_A1[dit]);
          BaseFab< Vector<Real>* >& QRmatrixFab = *(m_QRmatrix[dit]);
          for (BoxIterator bit(coarseGhostedBox); bit.ok(); ++bit)
            {
              IntVect ivc = bit();
              delete A1Fab(ivc, 0);
              delete QRmatrixFab(ivc, 0);
            }
          delete m_A1[dit];
          delete m_QRmatrix[dit];
        }
    }
  m_haveJ = false;
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderMappedCFInterp::setJ(
                                     const LevelData<FArrayBox>&  a_FnLevJ,
                                     const LevelData<FArrayBox>&  a_CrFnLevJ)
{
  CH_TIME("FourthOrderMappedCFInterp::setJ");
  // Free up existing arrays with J, if any.
  freeJ();
  // IntVect fineGhostVect = m_refineVect * m_coarseGhostVect;
  // FIXME: a_FnLevJ has 4 ghosts, but I need m_interpUnit more for gradient.
  IntVect fineGhostVect = m_refineVect * m_coarseGhostVect;
  IntVect fineGhostVect1 = fineGhostVect + m_interpUnit;
  m_JfineAvgPowers.define(m_layout, m_numNonzeroPowers, fineGhostVect1);
  m_FnLevJ.define(m_layout, 1, fineGhostVect1);
  DataIterator dit = m_layout.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& FnLevJFab = m_FnLevJ[dit];
      FnLevJFab.copy(a_FnLevJ[dit]);
      const Box& bx = m_layout[dit];
      const Box& bxCoarse = m_layoutCoarsened[dit];
      Box coarseGhostedBox = grow(bxCoarse, m_coarseGhostVect);
      Box fineGhostedBox = grow(bx, fineGhostVect);
      Box fineGhosted1Box = grow(bx, fineGhostVect1);
      FArrayBox& JfineAvgPowersFab = m_JfineAvgPowers[dit];
      FORT_REPEATPRODFAB(CHF_FRA(JfineAvgPowersFab),
                         CHF_CONST_FRA(FnLevJFab),
                         CHF_CONST_FRA(m_fineAvgPowers),
                         CHF_CONST_INTVECT(m_refineVect),
                         CHF_BOX(fineGhosted1Box));
      FArrayBox gradJFineFab(fineGhostedBox, m_nInterpDims);
      getGradient(gradJFineFab, FnLevJFab, fineGhostedBox);
      FArrayBox prodFab(fineGhostedBox, 1);
      Real scaling = 1./12.;
      int compDir = 0;
      for (int comp = 0; comp < m_numStates; comp++)
        {
          int interpDir = 0;
          for (int idir = 0; idir < SpaceDim; idir++)
            if ( !m_fixedDims.contains(idir) )
              {
                Interval interpIntvl(interpDir, interpDir);
                FArrayBox gradJFineFabDir(interpIntvl, gradJFineFab);
                Interval compIntvl(compDir, compDir);
                FArrayBox gradFineAvgPowersDir(compIntvl, m_gradFineAvgPowers);
                FORT_REPEATPRODFAB(CHF_FRA(prodFab),
                                   CHF_CONST_FRA(gradJFineFabDir),
                                   CHF_CONST_FRA(gradFineAvgPowersDir),
                                   CHF_CONST_INTVECT(m_refineVect),
                                   CHF_BOX(fineGhostedBox));
                JfineAvgPowersFab.plus(prodFab, scaling, 0, comp);
                interpDir++;
                compDir++;
              }
        }
    }

  m_alpha.define(m_layoutCoarsened, 1, m_coarseGhostVect);
  m_beta.define(m_layoutCoarsened, 1, m_coarseGhostVect);
  m_B.define(m_layoutCoarsened, m_numNonzeroPowers, m_coarseGhostVect);
  m_normB.define(m_layoutCoarsened, 1, m_coarseGhostVect);
  m_A1.define(m_layoutCoarsened);
  m_QRtau.define(m_layoutCoarsened, m_numNonzeroPowers, m_coarseGhostVect);
  m_QRmatrix.define(m_layoutCoarsened);
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bxCoarse = m_layoutCoarsened[dit];
      Box coarseGhostedBox = grow(bxCoarse, m_coarseGhostVect);
      const FArrayBox& JfineAvgPowersFab = m_JfineAvgPowers[dit];
      const FArrayBox& CrFnLevJFab = a_CrFnLevJ[dit];
      // Set BFab to sum of JfineAvgPowersFab over the fine cells within.
      FArrayBox& BFab = m_B[dit];
      FORT_AVGWITHIN(CHF_FRA(BFab),
                     CHF_CONST_FRA(JfineAvgPowersFab),
                     CHF_CONST_INTVECT(m_refineVect),
                     CHF_BOX(coarseGhostedBox),
                     CHF_BOX(m_fineBox));
      // At each coarse cell, set normBFab to
      // sqrt(<Jcoarse>^2 + sum_i of BFab[i]**2)
      FArrayBox& normBFab = m_normB[dit];
      FArrayBox& alphaFab = m_alpha[dit];
      FArrayBox& betaFab = m_beta[dit];
      FORT_GETMAPPEDCFINTERPARRAYS(CHF_FRA1(normBFab, 0),
                                   CHF_FRA1(alphaFab, 0),
                                   CHF_FRA1(betaFab, 0),
                                   CHF_CONST_FRA1(CrFnLevJFab, 0),
                                   CHF_CONST_FRA(BFab),
                                   CHF_BOX(coarseGhostedBox));

      const BaseFab<IntVect>& stencilHereFab = m_stencilHere[dit];
      m_A1[dit] = new BaseFab< Vector<Real>* >(coarseGhostedBox, 1);
      BaseFab< Vector<Real>* >& A1Fab = *(m_A1[dit]);
      m_QRmatrix[dit] = new BaseFab< Vector<Real>* >(coarseGhostedBox, 1);
      BaseFab< Vector<Real>* >& QRmatrixFab = *(m_QRmatrix[dit]);
      FArrayBox& QRtauFab = m_QRtau[dit];
      for (BoxIterator bit(coarseGhostedBox); bit.ok(); ++bit)
        {
          IntVect ivc = bit();
          IntVect stencilIndex = stencilHereFab(ivc, 0);
          const FourthOrderMappedCFInterpStencil& stencilAtCell =
            *(m_stencils(stencilIndex, 0));
          const Vector<Real>& avgPowers = stencilAtCell.avgPowers();
          int stencilSize = avgPowers.size() / m_numNonzeroPowers;

          Vector<Real> A1(stencilSize);
          Vector<Real> A2(stencilSize * m_numNonzeroPowers);
          // Now A2 has length stencilSize * m_numNonzeroPowers,
          // grouped by rows of length m_numNonzeroPowers.
          // Find QR decomposition of A2, and mathematically get:
          // QQTU matrix (P-1)xN, upper part of tranpose(QQ), which we save;
          // RU matrix with (P*(P-1))/2 nonzeroes, which we save.
          // Numerically we have QRmatrixFab(ivc, 0) and QRtauFab(ivc, :).
          // We overwrite A2 with QR decomposition information.
          FORT_GETMAPPEDCFINTERPQR(CHF_VR(A1),
                                   CHF_VR(A2),
                                   CHF_FRA(QRtauFab),
                                   CHF_CONST_VR(avgPowers),
                                   CHF_CONST_FRA1(alphaFab, 0),
                                   CHF_CONST_FRA1(betaFab, 0),
                                   CHF_CONST_FRA(BFab),
                                   CHF_CONST_INTVECT(ivc),
                                   CHF_CONST_INT(stencilSize),
                                   CHF_CONST_INT(m_numNonzeroPowers));
          A1Fab(ivc, 0) = new Vector<Real>(A1);
          QRmatrixFab(ivc, 0) = new Vector<Real>(A2);
        }
    }
  m_haveJ = true;
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderMappedCFInterp::setCrFnLevUandJU(
                                                 /// coarse <U>
                                                 const LevelData<FArrayBox>&   a_CrFnLevU,
                                                 /// coarse <JU>
                                                 const LevelData<FArrayBox>&   a_CrFnLevJU)
{
  CH_TIME("FourthOrderMappedCFInterp::setCrFnLevUandJU");
  CH_assert(m_defined);
  CH_assert(m_haveJ);

  m_coeff0.define(m_layoutCoarsened, m_numStates, m_coarseGhostVect);
  m_coeffs.define(m_layoutCoarsened, m_numStates * m_numNonzeroPowers, m_coarseGhostVect);

  // Simple copyTo also fills in ghost cells of m_coarseData.
  // a_CrFnLevU.copyTo(m_coarseData);
  DataIterator dit = m_layoutCoarsened.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      const Box& bxCoarse = m_layoutCoarsened[dit];
      Box coarseGhostedBox = grow(bxCoarse, m_coarseGhostVect);
      const FArrayBox& CrFnLevUFab = a_CrFnLevU[dit];
      const FArrayBox& CrFnLevJUFab = a_CrFnLevJU[dit];
      const FArrayBox& QRtauFab = m_QRtau[dit];
      const FArrayBox& BFab = m_B[dit];
      const FArrayBox& normBFab = m_normB[dit];
      const FArrayBox& alphaFab = m_alpha[dit];
      const FArrayBox& betaFab = m_beta[dit];
      FArrayBox& coeff0Fab = m_coeff0[dit];
      FArrayBox& coeffsFab = m_coeffs[dit];
      const BaseFab< Vector<Real>* >& QRmatrixFab = *(m_QRmatrix[dit]);
      const BaseFab< Vector<Real>* >& A1Fab = *(m_A1[dit]);
      for (BoxIterator bit(coarseGhostedBox); bit.ok(); ++bit)
        {
          IntVect ivc = bit();
          // We need <U> on every stencil cell of ivc.
          const BaseFab<IntVect>& stencilHereFab = m_stencilHere[dit];
          IntVect stencilIndex = stencilHereFab(ivc, 0);
          const FourthOrderMappedCFInterpStencil& stencilAtCell =
            *(m_stencils(stencilIndex, 0));
          const Vector<int>& stencilIndices = stencilAtCell.indicesFull();
          int stencilSize = stencilIndices.size() / SpaceDim;
          const Vector<Real>& QRmatrix = *(QRmatrixFab(ivc, 0));
          const Vector<Real>& A1 = *(A1Fab(ivc, 0));
          FORT_GETMAPPEDCFINTERPCOEFFS(CHF_FRA(coeff0Fab),
                                       CHF_FRA(coeffsFab),
                                       CHF_CONST_FRA(CrFnLevUFab),
                                       CHF_CONST_FRA(CrFnLevJUFab),
                                       CHF_CONST_FRA(BFab),
                                       CHF_CONST_FRA1(normBFab, 0),
                                       CHF_CONST_FRA1(alphaFab, 0),
                                       CHF_CONST_FRA1(betaFab, 0),
                                       CHF_CONST_INTVECT(ivc),
                                       CHF_CONST_VI(stencilIndices),
                                       CHF_CONST_INT(stencilSize),
                                       CHF_CONST_INT(m_numNonzeroPowers),
                                       CHF_CONST_INT(m_numStates),
                                       CHF_CONST_VR(QRmatrix),
                                       CHF_CONST_FRA(QRtauFab),
                                       CHF_CONST_VR(A1));
        }
    }
  m_haveCoeffs = true;
}


//////////////////////////////////////////////////////////////////////////////
void FourthOrderMappedCFInterp::interpToFineJU(
                                               /// interpolated fine <JU>
                                               LevelData<FArrayBox>&         a_FnLevJU)
{
  CH_TIME("FourthOrderMappedCFInterp::interpToFineJU");
  CH_assert(m_defined);
  CH_assert(m_haveJ);
  CH_assert(m_haveCoeffs);

  // Simple copyTo also fills in ghost cells of m_coarseData.
  // a_CrFnLevU.copyTo(m_coarseData);
  DataIterator dit = m_layoutCoarsened.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& FnLevJUFab = a_FnLevJU[dit];
      // const Box& bx = m_layout[dit];
      const Box& bx = FnLevJUFab.box();
      const FArrayBox& coeff0Fab = m_coeff0[dit];
      const FArrayBox& coeffsFab = m_coeffs[dit];
      const FArrayBox& FnLevJFab = m_FnLevJ[dit];
      const FArrayBox& JfineAvgPowersFab = m_JfineAvgPowers[dit];
      FORT_SOLVEMAPPEDCFINTERP(CHF_FRA(FnLevJUFab),
                               CHF_CONST_FRA(coeff0Fab),
                               CHF_CONST_FRA(coeffsFab),
                               CHF_CONST_FRA1(FnLevJFab, 0),
                               CHF_CONST_FRA(JfineAvgPowersFab),
                               CHF_CONST_INTVECT(m_refineVect),
                               CHF_CONST_INT(m_numStates),
                               CHF_CONST_INT(m_numNonzeroPowers),
                               CHF_BOX(bx));
      // dummy statement in order to get around gdb bug
      int dummy_unused = 0; dummy_unused = 0;
    }
}

#include "NamespaceFooter.H"
