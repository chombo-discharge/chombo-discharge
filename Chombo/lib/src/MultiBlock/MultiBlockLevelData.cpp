#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

// petermc, 10 Apr 2008

#include "MultiBlockLevelData.H"
#include "LoHiSide.H"
#include "BoxIterator.H"
#include "LoadBalance.H"
#include "MultiBlockCopier.H"
#include "IntegrationF_F.H"
#include "AMRNodeOpF_F.H"  // FORT_NODEGRAD, in AMRElliptic
#include "AMRPoissonOpF_F.H" // FORT_OPERATORLAP, in AMRElliptic
#include "CellCentersF_F.H"
#include "MomentsF_F.H"
#include "LeastSquaresF_F.H"
#include "LeastSquaresCallF.H"

#include "NamespaceHeader.H"

// ---------------------------------------------------------
MultiBlockLevelData::MultiBlockLevelData()
{
  m_isDefined = false;
  m_gotIntegralStuff = false;
  m_gotGradJacobianCellCenters = false;
  m_order = 2;
  m_ghostFactor = 1;
  m_useAverage = false;
  m_averageGhost = false;
  m_verbose = 0;
}


// ---------------------------------------------------------
MultiBlockLevelData::~MultiBlockLevelData()
{
  if (m_isDefined)
    {
      delete [] m_workArray;
    }
  if (m_gotIntegralStuff)
    {
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          delete m_neighborhoodSize[dit];
          delete m_moments[dit];
          delete m_ghostMoments[dit];
          delete m_neighborhood[dit];
        }
    }
}


// ---------------------------------------------------------
MultiBlockLevelData::MultiBlockLevelData(/// multiblocked mapped domain
                                         const MappedDomain&        a_mappedDomain,
                                         /// CELL-centered grids
                                         const DisjointBoxLayout&   a_grids,
                                         /// ghost width
                                         const IntVect&             a_ghost,
                                         /// radius of neighborhood from which to interpolate
                                         int                        a_radius,
                                         /// max degree of moments of displacement
                                         int                        a_degree)
{
  define(a_mappedDomain, a_grids, a_ghost, a_radius, a_degree);
}


// ---------------------------------------------------------
void MultiBlockLevelData::define(/// multiblocked mapped domain
                                 const MappedDomain&        a_mappedDomain,
                                 /// CELL-centered grids
                                 const DisjointBoxLayout&   a_grids,
                                 /// ghost width
                                 const IntVect&             a_ghost,
                                 /// radius of neighborhood from which to interpolate
                                 int                        a_radius,
                                 /// max degree of moments of displacement
                                 int                        a_degree)
{
  CH_TIME("MultiBlockLevelData::define");
  m_mappedDomain = a_mappedDomain;
  m_grids = a_grids;
  m_ghost = a_ghost;
  m_radius = a_radius;
  m_degree = a_degree;
  m_degreeBox = Box(IntVect::Zero, m_degree * IntVect::Unit);
  m_ghostMore = m_ghostFactor*m_ghost;

  m_neqns = 1;
  m_nvars = 1;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_neqns *= 2*m_radius + 1;
      m_nvars *= m_degree + idir + 1;
    }
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      m_nvars /= idir + 1;
    }
  if (m_verbose >= 2)
    {
      pout() << "m_neqns = " << m_neqns
             << ", m_nvars = " << m_nvars << endl;
    }
  m_lwork = FORT_LSWORKSIZE(&m_neqns, &m_nvars);
  m_workArray = new Real[m_lwork];

  m_block.define(m_grids);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const Box& bx = m_grids[dit];
      m_block[dit] = m_mappedDomain.findBlock(bx);
    }

  m_cellCentersValid.define(m_grids, SpaceDim, m_ghostMore);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& cellCentersValidFab = m_cellCentersValid[dit];
      const Box& bxGhosted = cellCentersValidFab.box();
      // Maybe I should move all this to MappedDomain.
      const MappedBlock& mb = m_mappedDomain.block(m_block[dit]);
      const BlockMap& bm = mb.map();
      bm.getCellCenterCoordinates(cellCentersValidFab, bxGhosted);
      // dummy statement in order to get around gdb bug
      { int dummy_unused = 0; dummy_unused = 0; }
    }

  // We haven't yet done exchange on m_cellCentersValid,
  // so ghost cells of m_cellCentersValid still contain coordinates
  // of centers of ghost cells.  We'll overwrite them with valid data
  // (coordinates of centers of valid cells of neighboring boxes) only
  // after copying them to m_cellCentersGhosted.

  // m_gridsGhostedInternal.deepCopy(m_grids);
  if (! (m_ghost == IntVect::Zero))
    { // Construct   LayoutData< IntVectSet > m_ghostBetweenBlockCells.
      m_cellCentersGhosted.define(m_grids, SpaceDim, m_ghost);
      m_ghostBetweenBlockCells.define(m_grids);
      const Box& wholeDomain = m_mappedDomain.domain();
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          m_cellCentersGhosted[dit].copy(m_cellCentersValid[dit]);

          const MappedBlock& mblock = m_mappedDomain.block(m_block[dit]);
          const Box& blockDomain = mblock.box();

          const Box& bxBase = m_grids[dit];
          Box bxGhostedInternal(bxBase);
          bxGhostedInternal.grow(m_ghost);
          // Remove ghost cells that are external.
          // You could, instead, loop over all the faces of mblock,
          // removing ghost cells on the BOUNDARY faces.
          bxGhostedInternal &= wholeDomain;
          // ivs is ghosted internal block with the block domain removed.
          IntVectSet ivs(bxGhostedInternal);
          if (m_useAverage && (m_order == 4)) // added by petermc, 14 May 2008
            ivs -= grow(blockDomain, -1);
          else
            ivs -= blockDomain;
          m_ghostBetweenBlockCells[dit] = ivs;
          // CH_assert(bxGhostedInternal.contains(ivs.minBox()));
        }
    }
  exchangeCopyOnly(m_cellCentersValid);
  //  m_gridsGhostedInternal.close();

  m_isDefined = true;
}


// ---------------------------------------------------------
void MultiBlockLevelData::setVerbose(int   a_verbose)
{
  m_verbose = a_verbose;
}


// ---------------------------------------------------------
void MultiBlockLevelData::setOrder(int   a_order)
{
  m_order = a_order;
  if (! ((m_order == 2) || (m_order == 4)) )
    {
      MayDay::Error("MultiBlockLevelData::setOrder:  a_order must be 2 or 4");
    }

}


// ---------------------------------------------------------
void MultiBlockLevelData::setGhostFactor(int   a_ghostFactor)
{
  m_ghostFactor = a_ghostFactor;
}


// ---------------------------------------------------------
void MultiBlockLevelData::setUseAverage(bool   a_useAverage)
{
  m_useAverage = a_useAverage;
}


// ---------------------------------------------------------
void MultiBlockLevelData::setAverageGhost(bool   a_averageGhost)
{
  m_averageGhost = a_averageGhost;
}


// ---------------------------------------------------------
void MultiBlockLevelData::setGaussQuadrature(int   a_integralLength)
{
  CH_TIME("MultiBlockLevelData::setGaussQuadrature");
  CH_assert(m_isDefined);
  m_integralLength = a_integralLength;
  m_mappedDomainRefinedForQuadrature = refine(m_mappedDomain, m_integralLength);
  m_quadratureMethod = GAUSS;

  m_gridsRefinedForQuadrature.deepCopy(m_grids);
  m_gridsRefinedForQuadrature.refine(m_integralLength);
  m_gridsRefinedForQuadrature.close();

  m_integralCellBox = Box(IntVect::Zero,
                          (m_integralLength - 1) * IntVect::Unit,
                          IndexType::TheCellType());
  m_physPowFab.define(m_integralCellBox, SpaceDim*m_degree);
  FArrayBox integWeightsFab(m_integralCellBox, 1);
  FArrayBox cartesianOriginCoordsFab(m_integralCellBox, SpaceDim);
  FORT_GAUSSQUADRATURE(CHF_FRA(cartesianOriginCoordsFab),
                       CHF_FRA1(integWeightsFab, 0),
                       CHF_BOX(m_integralCellBox),
                       CHF_CONST_REALVECT(RealVect::Zero), // low corner
                       CHF_CONST_REALVECT(RealVect::Unit)); // high corner
  // Set m_physIntegrationPointsValidAll and m_weightedJacobianValidAll.
  setIntegralPointsAndWeights(cartesianOriginCoordsFab, integWeightsFab);
}


// ---------------------------------------------------------
void MultiBlockLevelData::setNewtonCotesQuadrature(int   a_integralLength)
{
  CH_TIME("MultiBlockLevelData::setNewtonCotesQuadrature");
  CH_assert(m_isDefined);
  m_integralLength = a_integralLength;
  m_quadratureMethod = NEWTONCOTES;

  m_gridsIntegrationAll.deepCopy(m_grids);
  m_gridsIntegrationAll.grow(m_ghostMore);
  m_gridsIntegrationAll.refine(m_integralLength);
  m_gridsIntegrationAll.surroundingNodes();
  m_gridsIntegrationAll.close();

  m_integralCellBox = Box(IntVect::Zero,
                          m_integralLength * IntVect::Unit,
                          IndexType::TheNodeType());
  m_physPowFab.define(m_integralCellBox, SpaceDim*m_degree);
  FArrayBox integWeightsFab(m_integralCellBox, 1);
  IntVect integralRules = m_integralLength * IntVect::Unit;
  RealVect integralInterval = (1. / (m_integralLength*1.)) * RealVect::Unit;
  FORT_INTEGRATIONWEIGHTS(CHF_FRA1(integWeightsFab, 0),
                          CHF_BOX(m_integralCellBox),
                          CHF_CONST_REALVECT(integralInterval),
                          CHF_CONST_INTVECT(integralRules));
  FArrayBox cartesianOriginCoordsFab(m_integralCellBox, SpaceDim);
  FORT_INTEGRATIONPOINTS(CHF_FRA(cartesianOriginCoordsFab),
                         CHF_BOX(m_integralCellBox),
                         CHF_CONST_REALVECT(RealVect::Zero),
                         CHF_CONST_REALVECT(integralInterval));

  // Set m_physIntegrationPointsValidAll and m_weightedJacobianAll.
  setIntegralPointsAndWeights(cartesianOriginCoordsFab, integWeightsFab);
}


// ---------------------------------------------------------
void MultiBlockLevelData::setIntegralPointsAndWeights(/// cartesian coordinates within unit cell
                                                      const FArrayBox&   a_cartesianOriginCoordsFab,
                                                      /// weights of points
                                                      const FArrayBox&   a_integWeightsFab)
{
  CH_TIME("MultiBlockLevelData::setIntegralPointsAndWeights");
  // m_physIntegrationPointsValidAll.define(m_gridsIntegrationAll, SpaceDim);
  m_physIntegrationPointsValidAll.define(m_gridsRefinedForQuadrature, SpaceDim,
                                         m_integralLength * m_ghostMore);

  m_physIntegrationPointsGhostedAll.define(m_gridsRefinedForQuadrature, SpaceDim,
                                           m_integralLength * m_ghostMore);

  // m_weightedJacobianValidAll.define(m_gridsIntegrationAll, 1);
  m_weightedJacobianValidAll.define(m_gridsRefinedForQuadrature, 1,
                                    m_integralLength * m_ghostMore);

  FArrayBox cartesianCoordsFab(m_integralCellBox, SpaceDim);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& physIntegrationPointsValidFab =
        m_physIntegrationPointsValidAll[dit];
      FArrayBox& physIntegrationPointsGhostedFab =
        m_physIntegrationPointsGhostedAll[dit];
      FArrayBox& weightedJacobianValidFab = m_weightedJacobianValidAll[dit];
      const MappedBlock& mb = m_mappedDomain.block(m_block[dit]);
      const BlockMap& bm = mb.map();

      // changed by petermc, 14 May 2008, m_ghost to m_ghostMore
      const Box& bxCells = grow(m_grids[dit], m_ghostMore);
      // Do NOT do const Box& bxCells = physIntegrationPointsFab.box();
      // because bxCells consist of coarse cells, not fine ones.
      for (BoxIterator bitCell(bxCells); bitCell.ok(); ++bitCell)
        {
          IntVect ivCell = bitCell();
          cartesianCoordsFab.copy(a_cartesianOriginCoordsFab);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              cartesianCoordsFab.plus(ivCell[idir]*1., idir);
            }
          IntVect shiftCell = m_integralLength * ivCell;
          physIntegrationPointsValidFab.shift(-shiftCell);
          weightedJacobianValidFab.shift(-shiftCell);
          // Set m_physIntegrationPointsFab and m_weightedJacobianFab on
          // m_integralCellBox only.  (Will shift back later.)
          bm.setJacobian(weightedJacobianValidFab,
                         m_integralCellBox,
                         cartesianCoordsFab);
          weightedJacobianValidFab.mult(a_integWeightsFab,
                                        m_integralCellBox, 0, 0);
          bm.setPhysicalFromMap(physIntegrationPointsValidFab,
                                m_integralCellBox,
                                cartesianCoordsFab);
          physIntegrationPointsValidFab.shift(+shiftCell);
          weightedJacobianValidFab.shift(+shiftCell);
          //          physIntegrationPointsGhostedFab.shift(-shiftCell);
          //          bm.setPhysicalFromMap(physIntegrationPointsGhostedFab,
          //                                m_integralCellBox,
          //                                cartesianCoordsFab);
          //          physIntegrationPointsGhostedFab.shift(+shiftCell);
        }
      physIntegrationPointsGhostedFab.copy(physIntegrationPointsValidFab);
    }
  if (m_verbose >= 1) pout() << "MultiBlockLevelData got integral stuff" << endl;
  m_gotIntegralStuff = true;
  // Still haven't done exchange on m_weightedJacobianValidAll
  // or on m_physIntegrationPointsValidAll.
  setCellVolumes();
  if (m_averageGhost)
    {
      m_weightedJacobianGhostedAll.define(m_gridsRefinedForQuadrature, 1,
                                          m_integralLength * m_ghostMore);
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          m_weightedJacobianGhostedAll[dit].copy(m_weightedJacobianValidAll[dit]);
        }
    }
  exchangeCopyOnly(m_weightedJacobianValidAll, true);
  exchangeCopyOnly(m_physIntegrationPointsValidAll, true);
  // NO exchange on m_physIntegrationPointsGhostedAll.
  if (m_verbose >= 1) pout() << "MultiBlockLevelData got Jacobians and volumes" << endl;

  setAllMoments();
  if (m_verbose >= 1) pout() << "MultiBlockLevelData got all moments" << endl;
}


// ---------------------------------------------------------
void MultiBlockLevelData::setAllMoments()
{
  CH_TIME("MultiBlockLevelData::setAllMoments");
  CH_assert(m_gotIntegralStuff); // need volumes too
  m_neighborhoodSize.define(m_grids);
  m_moments.define(m_grids);
  m_ghostMoments.define(m_grids);
  m_neighborhood.define(m_grids);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      int iblock = m_block[dit];
      // const Box& bxGhostedInternal = m_gridsGhostedInternal[dit];
      // This should be const, but we'll be shifting it around.
      FArrayBox& physIntegrationPointsFab =
        m_physIntegrationPointsValidAll[dit];
      const FArrayBox& weightedJacobianValidFab =
        m_weightedJacobianValidAll[dit];
      const FArrayBox* weightedJacobianGhostedFabPtr = (m_averageGhost) ?
        &m_weightedJacobianValidAll[dit] : NULL;
      const FArrayBox* cellVolumesGhostedFabPtr = (m_averageGhost) ?
        &m_cellVolumesGhosted[dit] : NULL;
      // we should have called exchangeCopyOnly() on this.
      const FArrayBox& cellVolumesValidFab =
        m_cellVolumesValid[dit];

      // m_ghost is NOT zero, so we do have m_cellCentersGhosted.
      // On ghost cells, it contains centers of ghost cells.
      const FArrayBox& cellCentersGhostedFab = m_cellCentersGhosted[dit];
      // On ghost cells, this contains centers of corresponding valid cells.
      const FArrayBox& cellCentersValidFab = m_cellCentersValid[dit];
      // const Box& cellCentersValidBox = cellCentersValidFab.box();
      const IntVectSet& ivsGhostBetweenBlock = m_ghostBetweenBlockCells[dit];

      /*
        We're going to fill these in by calling setMoments.
        The destructor should delete them.
      */
      {
        CH_TIME("MultiBlockLevelData::setAllMoments_allocate");

        m_neighborhoodSize[dit] =
          new IVSFAB<int>(ivsGhostBetweenBlock, 1);

        m_moments[dit] =
          new IVSFAB<Real>(ivsGhostBetweenBlock, m_neqns*m_nvars);

        m_ghostMoments[dit] =
          new IVSFAB<Real>(ivsGhostBetweenBlock, m_nvars);

        m_neighborhood[dit] =
          new IVSFAB<IntVect>(ivsGhostBetweenBlock, m_neqns);
      }

      IVSFAB<int>& neighborhoodSize = *m_neighborhoodSize[dit];
      IVSFAB<Real>& moments = *m_moments[dit];
      IVSFAB<Real>& ghostMoments = *m_ghostMoments[dit];
      IVSFAB<IntVect>& neighborhood = *m_neighborhood[dit];

      for (IVSIterator ivsit(ivsGhostBetweenBlock); ivsit.ok(); ++ivsit)
        {
          IntVect ghostCell = ivsit();

          // CH_assert(bxGhostedInternal.contains(ghostCell));
          /*
            We will fill in calcCenterFab(ghostCell, 0), and optionally
            condFab(ghostCell, 0) and calcAvgFab(ghostCell, 0).
          */

          // physical coordinates of center of ghost cell
          RealVect ghostCellCenter;
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              ghostCellCenter[idir] = cellCentersGhostedFab(ghostCell, idir);
            }

          Vector<Real> momentsVec;
          Vector<Real> ghostMomentsVec;
          Vector<IntVect> neighborhoodVec;
          setMoments(momentsVec, ghostMomentsVec, neighborhoodVec,
                     iblock, ghostCell, ghostCellCenter, cellCentersValidFab,
                     cellVolumesValidFab, cellVolumesGhostedFabPtr,
                     weightedJacobianValidFab, weightedJacobianGhostedFabPtr,
                     physIntegrationPointsFab);
          {
            CH_TIME("MultiBlockLevelData::setAllMoments_shuffle");
            int npoints = neighborhoodVec.size();
            neighborhoodSize(ghostCell, 0) = npoints;
            for (int comp = 0; comp < npoints * m_nvars; comp++)
              {
                moments(ghostCell, comp) = momentsVec[comp];
              }
            if (m_averageGhost)
              for (int comp = 0; comp < m_nvars; comp++)
                {
                  ghostMoments(ghostCell, comp) = ghostMomentsVec[comp];
                }
            for (int comp = 0; comp < npoints; comp++)
              {
                neighborhood(ghostCell, comp) = neighborhoodVec[comp];
              }
          }
          // int npoints = moments.size() / m_nvars;
        } // end iteration over IntVectSet ivsGhostBetweenBlock
    } // end loop over grids
}


// ---------------------------------------------------------
void MultiBlockLevelData::integrateOnCells(BoxLayoutData<FArrayBox>&         a_integral,
                                           const BoxLayoutData<FArrayBox>&   a_function)
{
  CH_TIME("MultiBlockLevelData::integrateOnCells");
  CH_assert(m_gotIntegralStuff);
  FArrayBox funOnCellFab(m_integralCellBox, 1);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& integralFab = a_integral[dit];
      const FArrayBox& functionFab = a_function[dit];
      const Box& bxCells = integralFab.box();
      const FArrayBox& weightedJacobianFab = (m_averageGhost) ?
        m_weightedJacobianGhostedAll[dit] : m_weightedJacobianValidAll[dit];
      for (BoxIterator bitCell(bxCells); bitCell.ok(); ++bitCell)
        {
          IntVect ivCell = bitCell();
          IntVect shiftCell = m_integralLength * ivCell;
          Box shiftedIntegrationPointsBox = m_integralCellBox + shiftCell;

          funOnCellFab.copy(weightedJacobianFab,
                            shiftedIntegrationPointsBox, 0, // source
                            m_integralCellBox, 0, // dest
                            1); // 1 component

          funOnCellFab.mult(functionFab,
                            shiftedIntegrationPointsBox, // source subbox
                            m_integralCellBox, // dest subbox
                            0, 0); // source and dest components

          integralFab(ivCell, 0) = funOnCellFab.sum(0);
        }
    }
}


// ---------------------------------------------------------
void MultiBlockLevelData::setGradJacobianCellCenters()
{
  CH_TIME("MultiBlockLevelData::setGradJacobianCellCenters");
  m_gradJacobianCellCenters.define(m_grids, SpaceDim);
  // Keep a ghost layer for exchange.
  // Why not just set up an FArrayBox with an extra layer where we recompute?
  // Because if ghost layer is 0, it's a pain to get the cell centers.
  LevelData<FArrayBox> jacobianCellCenters(m_grids, 1, IntVect::Unit);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      const MappedBlock& mb = m_mappedDomain.block(m_block[dit]);
      const BlockMap& bm = mb.map();

      const Box& validBox = m_grids[dit];
      FArrayBox& jacobianCellCentersFab = jacobianCellCenters[dit];

      FArrayBox cartesianCentersFab(validBox, SpaceDim);
      FORT_SETCELLCENTERS(CHF_FRA(cartesianCentersFab),
                          CHF_BOX(validBox));
      bm.setJacobian(jacobianCellCentersFab,
                     validBox,
                     cartesianCentersFab);
    }
  // This does NOT do the right thing across block boundaries, but that's OK,
  // because you're not going to use the results at block boundaries.
  jacobianCellCenters.exchange();

  Real one = 1.;
  // Now fill in m_gradJacobianCellCenters.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& gradJacobianCellCentersFab = m_gradJacobianCellCenters[dit];
      const FArrayBox& jacobianCellCentersFab = jacobianCellCenters[dit];
      const Box& validBox = m_grids[dit];

      const Box& blockDomain = m_mappedDomain.block(m_block[dit]).box();
      Box interiorBlockDomain = grow(blockDomain, -1);
      Box interiorBox = validBox & interiorBlockDomain;

      // Check that box containing gradients contains what you need of jacobianCellCentersFab.
      CH_assert(jacobianCellCentersFab.box().contains(grow(gradJacobianCellCentersFab.box(), 1)));
      FORT_NODEGRAD(CHF_FRA(gradJacobianCellCentersFab),
                    CHF_CONST_FRA(jacobianCellCentersFab),
                    CHF_BOX(interiorBox),
                    CHF_CONST_REAL(one));
    }
  m_gotGradJacobianCellCenters = true;
}


// ---------------------------------------------------------
void MultiBlockLevelData::interpCenterFromAvg(/// function values on CELL centers, predefined and filled in by this function
                                              LevelData<FArrayBox>&         a_funCenter,
                                              /// averaged function values over CELLs
                                              const LevelData<FArrayBox>&   a_funAvg)
{
  CH_TIME("MultiBlockLevelData::interpCenterFromAvg");
  // base boxes of valid cells (excluding ghost cells)
  const DisjointBoxLayout& validLayout = a_funAvg.disjointBoxLayout();
  if (!m_useAverage || m_order == 2)
    {
      for (DataIterator dit = a_funCenter.dataIterator(); dit.ok(); ++dit)
        { // Fill in valid cells only.
          a_funCenter[dit].copy(a_funAvg[dit], validLayout[dit]);
        }
      // Fill in ghost cells:  We DO want the block interiors filled.
      // Ghost cells between boxes in the same block will be filled.
      // May not be necessary because we'll call exchange() after this anyway,
      // within exchangeAvg().
      // a_funCenter.exchange();
    }
  else if (m_order == 4)
    {
      if (!m_gotGradJacobianCellCenters)
        {
          setGradJacobianCellCenters();
        }
      // These are used for FORT_OPERATORLAP.
      Real alpha = 0.;
      Real beta = 1.;
      Real one = 1.;
      bool ghostSufficient = true;
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          if (m_ghost[idir] == 0) ghostSufficient = false;
        }
      LevelData<FArrayBox>* funAvgPtr;
      if (ghostSufficient)
        {
          funAvgPtr = (LevelData<FArrayBox>*) &a_funAvg;
        }
      else
        {
          funAvgPtr = new LevelData<FArrayBox>(m_grids, 1, IntVect::Unit);
          a_funAvg.copyTo(*funAvgPtr);
        }
      for (DataIterator dit = a_funCenter.dataIterator(); dit.ok(); ++dit)
        {
          FArrayBox& funCenterFab = a_funCenter[dit];
          const FArrayBox& funAvgFab = funAvgPtr->operator[](dit);
          const FArrayBox& cellVolumesValidFab = m_cellVolumesValid[dit];

          const Box& validBox = m_grids[dit];
          const Box& blockDomain = m_mappedDomain.block(m_block[dit]).box();
          Box interiorBlockDomain = grow(blockDomain, -1);
          Box interiorBox = validBox & interiorBlockDomain;

          funCenterFab.copy(funAvgFab, interiorBox);

          FArrayBox gradFunFab(interiorBox, SpaceDim);
          // Check that box containing gradients contains what you need of gradFunFab.
          CH_assert(funAvgFab.box().contains(grow(gradFunFab.box(), 1)));
          FORT_NODEGRAD(CHF_FRA(gradFunFab),
                        CHF_CONST_FRA(funAvgFab),
                        CHF_BOX(interiorBox),
                        CHF_CONST_REAL(one));

          const FArrayBox& gradJacobianCellCentersFab =
            m_gradJacobianCellCenters[dit];
          // Set termFab = gradJacobianCellCentersFab dot gradFunFab / 12.
          FArrayBox termFab(interiorBox, 1);
          FArrayBox prodFab(interiorBox, 1);
          termFab.setVal(0.);
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              prodFab.copy(gradJacobianCellCentersFab, idir, 0);
              prodFab.mult(gradFunFab, idir, 0);
              termFab += prodFab;
            }
          termFab *= (-1./12.);
          // termFab /= m_cellVolumes[dit];
          termFab.divide(cellVolumesValidFab, interiorBox, 0, 0); // src, dest
          // funCenterFab.plus(termFab);
          funCenterFab.plus(termFab, interiorBox, 0, 0); // src, dest

          FORT_OPERATORLAP(CHF_FRA(termFab),
                           CHF_CONST_FRA(funAvgFab),
                           CHF_BOX(interiorBox),
                           CHF_CONST_REAL(one),
                           CHF_CONST_REAL(alpha),
                           CHF_CONST_REAL(beta));
          termFab *= -(1./24.);
          funCenterFab.plus(termFab, interiorBox, 0, 0); // src, dest
        }
      if (!ghostSufficient)
        {
          delete funAvgPtr;
        }
      // Fill in ghost cells:  We DO want the block interiors filled.
      // Ghost cells between boxes in the same block will be filled.
      // May not be necessary because we'll call exchange() after this anyway,
      // within exchangeAvg().
      // a_funCenter.exchange();
    }
  else
    {
      MayDay::Error("MultiBlockLevelData::interpCenterFromAvg:  m_order must be 2 or 4");
    }
}


// ---------------------------------------------------------
void MultiBlockLevelData::setCellVolumes()
{
  CH_TIME("MultiBlockLevelData::setCellVolumes");
  CH_assert(m_gotIntegralStuff);
  LevelData<FArrayBox> ones(m_gridsRefinedForQuadrature, 1,
                            m_integralLength * m_ghostMore);
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
    {
      ones[dit].setVal(1.);
    }
  // petermc, 14 May 2008:  changed from m_ghost to m_ghostMore
  m_cellVolumesValid.define(m_grids, 1, m_ghostMore);
  integrateOnCells(m_cellVolumesValid, ones);
  // Note that we haven't yet done exchange on m_cellVolumesValid.
  // We've filled in ghost cells and want to retain them iff m_averageGhost.
  if (m_averageGhost)
    {
      m_cellVolumesGhosted.define(m_grids, 1, m_ghostMore);
      for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit)
        {
          m_cellVolumesGhosted[dit].copy(m_cellVolumesValid[dit]);
        }
    }
  // Now that we've saved the ghost cells in m_cellVolumesGhosted,
  // overwrite them with valid cells.
  exchangeCopyOnly(m_cellVolumesValid);
}


// ---------------------------------------------------------
const LevelData<FArrayBox>& MultiBlockLevelData::cellCenters()
{
  CH_assert(m_isDefined);
  return
    ((m_ghost == IntVect::Zero) ? m_cellCentersValid : m_cellCentersGhosted);
}


// ---------------------------------------------------------
const LevelData<FArrayBox>& MultiBlockLevelData::physIntegrationPointsValidAll()
{
  CH_assert(m_gotIntegralStuff);
  return m_physIntegrationPointsValidAll;
}


// ---------------------------------------------------------
const LevelData<FArrayBox>& MultiBlockLevelData::physIntegrationPointsGhostedAll()
{
  CH_assert(m_gotIntegralStuff);
  return m_physIntegrationPointsGhostedAll;
}


// ---------------------------------------------------------
const LevelData<FArrayBox>& MultiBlockLevelData::cellVolumesValid()
{
  CH_assert(m_gotIntegralStuff);
  return m_cellVolumesValid;
}


// ---------------------------------------------------------
const LevelData<FArrayBox>& MultiBlockLevelData::cellVolumesGhosted()
{
  CH_assert(m_gotIntegralStuff);
  CH_assert(m_averageGhost);
  return m_cellVolumesGhosted;
}


// ---------------------------------------------------------
void MultiBlockLevelData::exchangeCopyOnly(LevelData<FArrayBox>&   a_data,
                                           bool                    a_refine)
{
  exchangeCopyOnly(a_data, a_data.interval(), a_refine);
}


// ---------------------------------------------------------
void MultiBlockLevelData::exchangeCopyOnly(LevelData<FArrayBox>&   a_data,
                                           const Interval&         a_intvl,
                                           bool                    a_refine)
{
  CH_TIME("MultiBlockLevelData::exchangeCopyOnly");
  const MappedDomain& md = (a_refine) ?
    m_mappedDomainRefinedForQuadrature : m_mappedDomain;
  const DisjointBoxLayout& grids = (a_refine) ?
    m_gridsRefinedForQuadrature : m_grids;
  IntVect theseGhosts = (a_refine) ?
    (m_integralLength * m_ghostMore) : m_ghostMore;

  a_data.exchange(a_intvl);
  int ncomps = a_intvl.size();
  Interval intvlNew(0, ncomps-1);
  int nblocks = md.numBlocks();

  for (int iblock = 0; iblock < nblocks; iblock++)
    { // Copy valid a_data within iblock to ghost cells of other blocks.
      const MappedBlock& mblock = md.block(iblock);
      const Box& blockBox = mblock.box();
      int faceID = 0;
      for (SideIterator sit; sit.ok(); sit.next())
        {
          Side::LoHiSide side = sit();
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              const BlockBoundary& bb = mblock.boundary(faceID);
              if (bb.type() == BlockBoundary::CONFORMAL)
                { // This face of this block is a CONFORMAL boundary.
                  int ghostInDir = theseGhosts[idir];
                  if (ghostInDir > 0)
                    { // Set boundaryCells to the valid cells in this block
                      // that are within ghostInDir of this boundary.
                      Box boundaryCells =
                        adjCellBox(blockBox, idir, side, ghostInDir);
                      boundaryCells.shift(idir, -sign(side)*ghostInDir); // within

                      const Box& destBlockBox =
                        md.block(bb.neighbor()).box();

                      // Set validSubboxes to the nonempty intersections
                      // of all boxes with boundaryCells.
                      // We could instead use the SINGLE box boundaryCells,
                      // but then if the source block isn't covered by boxes,
                      // we'd be copying undefined data around.
                      Vector<Box> validSubboxes;
                      for (LayoutIterator lit = grids.layoutIterator(); lit.ok(); ++lit)
                        {
                          const Box& bx = grids[lit];
                          Box boundaryOverlap = bx & boundaryCells;
                          if (!boundaryOverlap.isEmpty())
                            {
                              validSubboxes.push_back(boundaryOverlap);
                            }
                        }
                      if (validSubboxes.size() > 0)
                        {
                          Vector<int> procs;
                          int eekflag = LoadBalance(procs, validSubboxes);
                          if (eekflag != 0)
                            {
                              cerr << "LoadBalance returned error code " << eekflag
                                   << endl;
                              return;
                            }
                          DisjointBoxLayout validLayout(validSubboxes, procs);
                          LevelData<FArrayBox> validData(validLayout, ncomps);
                          // Copy original data to validData, which is on
                          // valid cells next to this face of this block.
                          a_data.copyTo(a_intvl, validData, intvlNew);

                          // ghostLayout is validLayout remapped
                          // according to the neighboring block.
                          // It consists of ghost cells of m_grids.
                          DisjointBoxLayout ghostLayout;
                          ghostLayout.deepCopy(validLayout);
                          ghostLayout.convertOldToNew(bb.getPermutation(), bb.getSign(), bb.getTranslation());
                          ghostLayout.close();
                          // Now copy validData to ghostData according to
                          // the boundary transformation.
                          // Copying is local because ghostLayout was
                          // constructed from validLayout.
                          LevelData<FArrayBox> ghostData(ghostLayout, ncomps);
                          for (DataIterator dit = ghostData.dataIterator(); dit.ok(); ++dit)
                            {
                              const Box& validBox = validLayout[dit];
                              FArrayBox& ghostFab = ghostData[dit];
                              const FArrayBox& validFab = validData[dit];
                              for (BoxIterator bitValid(validBox); bitValid.ok(); ++bitValid)
                                {
                                  IntVect ivOld = bitValid();
                                  IntVect ivNew = bb.convertOldToNew(ivOld);
                                  for (int comp = 0; comp < ncomps; comp++)
                                    {
                                      ghostFab(ivNew, comp) =
                                        validFab(ivOld, comp);
                                    }
                                }
                            }
                          // Finally copy ghostData to a_data.
                          // Note that ghostData has no ghost cells.
                          // Ghost cells of a_data will be filled in
                          // with ghostData.

                          Box boundaryCellsConverted(bb.convertOldToNew(boundaryCells));

                          // I wonder if you could do this just by
                          // taking the correct ghost IntVect in Copier.
                          MultiBlockCopier mbCopier(ghostLayout,
                                                    grids,
                                                    theseGhosts,
                                                    boundaryCellsConverted, // source
                                                    destBlockBox);
                          ghostData.copyTo(intvlNew, a_data, a_intvl,
                                           mbCopier);
                          // ghostData.copyTo(intvlNew, a_data, a_intvl);
                        } // if at least one valid box on this boundary
                    } // if (ghostInDir > 0)
                } // if boundary is conformal
              faceID++;
            } // iterate over dimensions
        } // iterate over sides
    }
}


// ---------------------------------------------------------
void MultiBlockLevelData::exchangeAll(LevelData<FArrayBox>&   a_calcCenter,
                                      LevelData<FArrayBox>*   a_condPtr,
                                      LevelData<FArrayBox>*   a_calcAvgPtr)
{
  CH_TIME("MultiBlockLevelData::exchangeAll");
  //  if (m_useAverage)
  //    exchangeCopyOnly(*a_calcAvgPtr);
  //  else
  //    exchangeCopyOnly(a_calcCenter);
  exchangeCopyOnly(a_calcCenter);
  if (a_calcAvgPtr != NULL) exchangeCopyOnly(*a_calcAvgPtr);

  fillGhosts(a_calcCenter, a_condPtr, a_calcAvgPtr);
}


// ---------------------------------------------------------
void MultiBlockLevelData::fillGhosts(/// function values at cell centers
                                     LevelData<FArrayBox>&   a_calcCenter,
                                     /// condition numbers from least squares
                                     LevelData<FArrayBox>*   a_condPtr,
                                     /// interpolated averages on ghost cells (optional)
                                     LevelData<FArrayBox>*   a_calcAvgPtr)
{
  CH_TIME("MultiBlockLevelData::fillGhosts");
  for (DataIterator dit = a_calcCenter.dataIterator(); dit.ok(); ++dit)
    {
      /*
        We will fill in calcCenterFab, and optionally condFab and calcAvgFab.
      */
      FArrayBox& calcCenterFab = a_calcCenter[dit];
      FArrayBox& condFab = (a_condPtr == NULL) ?
        calcCenterFab : // in this case, not used
        a_condPtr->operator[](dit);
      if (a_condPtr != NULL)
        condFab.setVal(0.);
      FArrayBox& calcAvgFab = (a_calcAvgPtr == NULL) ?
        calcCenterFab : // in this case, not used
        a_calcAvgPtr->operator[](dit);

      const Box& bxBase = m_grids[dit];

      // Save calcCenterFab in calcValidCenterFab:
      // we'll use it on ghost cells (valid data from other boxes).
      // We must save it because we'll overwrite calcCenterFab on ghost cells.
      //      FArrayBox calcValidCenterFab(calcCenterFab.box(), 1);
      //      calcValidCenterFab.copy(calcCenterFab);
      FArrayBox& calcUseFab = (m_useAverage) ? calcAvgFab : calcCenterFab;
      // Save a copy of calcUseFab because its ghost cells will be overwritten.
      FArrayBox calcValidFab(calcUseFab.box(), 1);
      calcValidFab.copy(calcUseFab);

      const IVSFAB<int>& neighborhoodSize = *m_neighborhoodSize[dit];
      const IVSFAB<Real>& moments = *m_moments[dit];
      const IVSFAB<Real>& ghostMoments = *m_ghostMoments[dit];
      const IVSFAB<IntVect>& neighborhood = *m_neighborhood[dit];

      const IntVectSet& ivsGhostBetweenBlock = m_ghostBetweenBlockCells[dit];
      for (IVSIterator ivsit(ivsGhostBetweenBlock); ivsit.ok(); ++ivsit)
        {
          IntVect ghostCell = ivsit();

          // CH_assert(bxGhostedInternal.contains(ghostCell));
          /*
            We will fill in calcCenterFab(ghostCell, 0), and optionally
            condFab(ghostCell, 0) and calcAvgFab(ghostCell, 0).
          */

          // int npoints = neighborhoodVec.size();
          int npoints = neighborhoodSize(ghostCell, 0);
          Vector<Real> momentsVec(npoints*m_nvars);
          Vector<Real> ghostMomentsVec(m_nvars);
          // int npoints = momentsVec.size() / m_nvars;
          Vector<Real> calcValidVec(npoints);
          {
            CH_TIME("MultiBlockLevelData::fillGhosts_shuffle");
            for (int comp = 0; comp < npoints * m_nvars; comp++)
              {
                momentsVec[comp] = moments(ghostCell, comp);
              }
            if (m_averageGhost)
              for (int comp = 0; comp < m_nvars; comp++)
                {
                  ghostMomentsVec[comp] = ghostMoments(ghostCell, comp);
                }
            for (int inbr = 0; inbr < npoints; inbr++)
              {
                const IntVect& ivNbr = neighborhood(ghostCell, inbr);
                calcValidVec[inbr] = calcValidFab(ivNbr, 0);
              }
          }
          // Find the average over ghostCell if and only if
          // m_averageGhost is set AND this is NOT a valid cell.
          // (When m_averageGhost is set, we still set the centers
          // on valid cells that are adjacent to block domain boundaries.)
          bool whetherAvgGhost =
            ((m_averageGhost) && ! (bxBase.contains(ghostCell)));
          int flagAvgGhost = whetherAvgGhost ? 1 : 0;

          // FORT_LSGHOST will return these.
          Real calcCenter, cond, calcAvg;
          FORT_LSGHOST(CHF_REAL(calcCenter),
                       CHF_REAL(cond),
                       CHF_REAL(calcAvg),
                       CHF_CONST_VR(momentsVec),
                       CHF_CONST_VR(ghostMomentsVec),
                       CHF_CONST_VR(calcValidVec),
                       CHF_CONST_INT(flagAvgGhost),
                       CHF_CONST_INT(npoints),
                       CHF_CONST_INT(m_nvars),
                       CHF_CONST_INT(m_lwork),
                       CHF_R1D(m_workArray, m_lwork));
          calcCenterFab(ghostCell, 0) = calcCenter;
          if (a_condPtr != NULL) condFab(ghostCell, 0) = cond;
          if (a_calcAvgPtr != NULL) calcAvgFab(ghostCell, 0) = calcAvg;
        } // end iteration over IntVectSet ivsGhostBetweenBlock
    } // end loop over grids
}



// ---------------------------------------------------------
void MultiBlockLevelData::setMoments(Vector<Real>&      a_moments,
                                     Vector<Real>&      a_ghostMoments,
                                     Vector<IntVect>&   a_neighborhood,
                                     int                a_block,
                                     const IntVect&     a_ghostCell,
                                     const RealVect&    a_ghostCellCenter,
                                     const FArrayBox&   a_cellCentersValidFab,
                                     const FArrayBox&   a_cellVolumesValidFab,
                                     const FArrayBox*   a_cellVolumesGhostedFabPtr,
                                     const FArrayBox&   a_weightedJacobianValidFab,
                                     const FArrayBox*   a_weightedJacobianGhostedFabPtr,
                                     FArrayBox&         a_physIntegrationPointsFab)
{
  CH_TIME("MultiBlockLevelData::setMoments");
  const Box& cellCentersValidBox = a_cellCentersValidFab.box();
  const Box& blockDomain = m_mappedDomain.block(a_block).box();

  // block containing ghostCellCenter
  int iblockOther = m_mappedDomain.findBlock(a_ghostCellCenter);

  const MappedBlock& mbOther = m_mappedDomain.block(iblockOther);
  const Box& domainOther = mbOther.box();
  const BlockMap& bmOther = mbOther.map();
  // Cartesian coordinates, in iblockOther, of ghost cell center
  RealVect cartesianOther(a_ghostCellCenter);
  {
    CH_TIME("MultiBlockLevelData::setMoments_realToCartesian");
    bmOther.realToCartesian(cartesianOther);
  }

  // cell in iblockOther containing ghost cell center,
  // in iblockOther's index space
  IntVect indicesCenterOther;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      indicesCenterOther[idir] = int(floor(cartesianOther[idir]));
      // Added 5 May 2008:  indicesCenterOther should be a valid cell.
      if (indicesCenterOther[idir] < domainOther.smallEnd(idir))
        indicesCenterOther[idir] = domainOther.smallEnd(idir);
      else if (indicesCenterOther[idir] > domainOther.bigEnd(idir))
        indicesCenterOther[idir] = domainOther.bigEnd(idir);
    }

  // cell in iblockOther containing ghost cell center,
  // exchanged into a_block's index space
  IntVect indicesInOther =
    m_mappedDomain.convertBetweenBlocks(iblockOther,
                                        a_block,
                                        indicesCenterOther);

  /*
    Get indices of cells in the neighborhood of indicesCenterOther,
    in block iblockOther.
  */
  DenseIntVectSet ivsOtherNbrs;
  IntVect centerNbrs;
  m_mappedDomain.getNeighborhood(ivsOtherNbrs,
                                 centerNbrs,
                                 indicesCenterOther,
                                 iblockOther,
                                 m_radius);

  /*
    a_neighborhood is ivsOtherNbrs converted from iblockOther to a_block.
    Then we can use a_neighborhood as indices of cells in the current box.

    Note that the cells of a_neighborhood are not necessarily in
    a box of radius m_radius like ivsOtherNbrs.
    Reason:  converting to central block, some cells in a_neighborhood
    may be in block to left and some may be in block below.
   */
  {
    CH_TIME("MultiBlockLevelData::setMoments_neighborhood");
    a_neighborhood.clear();
    for (DenseIntVectSetIterator ivsit(ivsOtherNbrs); ivsit.ok(); ++ivsit)
      {
        IntVect ivOtherNbr = ivsit();
        IntVect ivNbr; // return this in convertGhostBlocks
        bool found =
          m_mappedDomain.convertGhostBlocks(ivNbr, a_block, // new
                                            ivOtherNbr, iblockOther); //old
        if (found)
          {
            CH_TIME("MultiBlockLevelData::setMoments_neighborhood_add");
            // a_neighborhood |= ivNbr;
            bool alreadyHave = false;
            for (int inbr = 0; inbr < a_neighborhood.size(); inbr++)
              {
                if (a_neighborhood[inbr] == ivNbr)
                  {
                    alreadyHave = true;
                    break;
                  }
              }
            if (!alreadyHave)
              a_neighborhood.push_back(ivNbr);
          }
      }
  }

  /*
    Get magnitudes of gradients of X, Y, Z on ghostCell,
    where X, Y, Z are the maps on the valid block.
  */
  RealVect physicalCellWidth;
  {
    CH_TIME("MultiBlockLevelData::setMoments_physicalCellWidth");
    for (int idir = 0; idir < SpaceDim; idir++)
      {
        Real gradMag2 = 0.;
        for (int idirMap = 0; idirMap < SpaceDim; idirMap++)
          {
            Real fac = 0.5;

            IntVect nextLo = indicesInOther - BASISV(idirMap);
            if (! (cellCentersValidBox.contains(nextLo)) ||
                forbiddenZone(nextLo, blockDomain) )
              {
                nextLo = indicesInOther;
                fac = 1.;
              }

            IntVect nextHi = indicesInOther + BASISV(idirMap);
            if (! (cellCentersValidBox.contains(nextHi)) ||
                forbiddenZone(nextHi, blockDomain) )
              {
                nextHi = indicesInOther;
                fac = 1.;
              }

            Real diff = fac *
              (a_cellCentersValidFab(nextHi, idir) -
               a_cellCentersValidFab(nextLo, idir));
            gradMag2 += (diff*diff);
          }
        physicalCellWidth[idir] = sqrt(gradMag2);
      }
  }

  // Fill m_moments, and optionally a_ghostMoments.
  // setMoments(physicalCellWidth, blockNbrs, indicesNbrs);
  // We should really look into saving the moments in an IVSFAB<Real>.
  // FArrayBox momentsFab(bxNbrs, m_nvars);
  int npoints = a_neighborhood.size();
  if ((m_verbose >= 2) && (npoints < m_neqns))
    {
      pout() << "For block " << a_block
             << " cell " << a_ghostCell
             << " use only " << npoints
             << " equations." << endl;
    }

  // In a_moments, each contiguous section of m_nvars components
  // corresponds to a point in the neighborhood.
  if (!m_useAverage || m_order == 2)
    {
      CH_TIME("MultiBlockLevelData::setMoments_order_2");

      // physCoordsDiff will be on the order of 1 or 2,
      // because we're dividing by physical cell width.
      Vector<Real> physCoordsDiff(npoints * SpaceDim);
      // a_calcValidCenterVec.resize(npoints);
      for (int inbr = 0; inbr < npoints; inbr++)
        {
          const IntVect& ivNbr = a_neighborhood[inbr];
          for (int idir = 0; idir < SpaceDim; idir++)
            {
              // displacement from center of ghostCell
              // to center of neighboring valid cell
              physCoordsDiff[SpaceDim*inbr + idir] =
                (a_cellCentersValidFab(ivNbr, idir) -
                 a_ghostCellCenter[idir]) /
                physicalCellWidth[idir];
            }
          // a_calcValidCenterVec[inbr] = a_calcValidCenterFab(ivNbr, 0);
        }

      /*
        To second order, moments of displacement are merely
        products of powers of components of displacement.
        Displacements are from a_ghostCellCenter to a_cellCentersValidFab,
        divided by physicalCellWidth.
      */

      a_moments.resize(npoints * m_nvars);
      //   m_moments = new Real[m_neqns * m_nvars];
      FORT_SETMOMENTS(CHF_VR(a_moments),
                      CHF_CONST_VR(physCoordsDiff),
                      CHF_CONST_INT(m_degree),
                      CHF_CONST_INT(m_nvars),
                      CHF_BOX(m_degreeBox));
    } // end if order 2
  else if (m_order == 4)
    { // order 4
      CH_TIME("MultiBlockLevelData::setMoments_order_4");
      a_moments.clear();
      Vector<Real> momentsNbr(m_nvars);
      // m_integralCellBox indexes integration points on the cell [0].
      // The first SpaceDim coordinates of m_physPowFab
      // will list their physical coordinates,
      // translated from a_ghostCell to [0].
      for (int inbr = 0; inbr < npoints; inbr++)
        {
          IntVect ivNbr = a_neighborhood[inbr];

          // Recall we have physical coordinates RealVect ghostCellCenter.
          IntVect shiftCell = m_integralLength * ivNbr;

          // Shift from m_integralCellBox + shiftCell
          // to m_integralCellBox.
          a_physIntegrationPointsFab.shift(-shiftCell);

          // powers of displacement from ghostCellCenter to
          // integration points on this ghost cell
          // m_physPowFab(m_integralCellBox, SpaceDim*(m_degree + 1));
          FORT_SETDISPLACEMENTPOWERS(CHF_FRA(m_physPowFab),
                                     CHF_CONST_FRA(a_physIntegrationPointsFab),
                                     CHF_BOX(m_integralCellBox),
                                     CHF_CONST_REALVECT(a_ghostCellCenter),
                                     CHF_CONST_REALVECT(physicalCellWidth),
                                     CHF_CONST_INT(m_degree));
          a_physIntegrationPointsFab.shift(+shiftCell);

          Box shiftedIntegralCellBox = m_integralCellBox + shiftCell;
          m_physPowFab.shift(+shiftCell);

          Real volNbr = a_cellVolumesValidFab(ivNbr, 0);
          FORT_SETVECTORAVGMOMENTS1PLUS(CHF_VR(momentsNbr),
                                        CHF_CONST_FRA(m_physPowFab),
                                        CHF_CONST_FRA1(a_weightedJacobianValidFab, 0),
                                        CHF_CONST_REAL(volNbr),
                                        CHF_BOX(shiftedIntegralCellBox),
                                        CHF_CONST_INT(m_degree),
                                        CHF_CONST_INT(m_nvars),
                                        CHF_BOX(m_degreeBox));
          a_moments.append(momentsNbr);
          // return indices of m_physPowFab to original
          m_physPowFab.shift(-shiftCell);
        } // end loop over neighborhood
    } // end if order 4

  // Recall that we have valid data in a_calcValidCenterFab.

  a_ghostMoments.resize(m_nvars);

  // whether to return average over ghost cell
  if (m_averageGhost)
    {
      CH_TIME("MultiBlockLevelData::setMoments_averageGhost");
      // m_integralCellBox indexes integration points on the cell [0].
      // The first SpaceDim coordinates of m_physPowFab
      // will list their physical coordinates,
      // translated from a_ghostCell to [0].
      // Recall we have physical coordinates RealVect ghostCellCenter.
      IntVect shiftCell = m_integralLength * a_ghostCell;
      m_physPowFab.copy(a_physIntegrationPointsFab,
                        m_integralCellBox + shiftCell, // source box
                        0, // first source component
                        m_integralCellBox, // dest box
                        0, // first dest component
                        SpaceDim); // number of components
      for (int idir = 0; idir < SpaceDim; idir++)
        {
          m_physPowFab.plus(-a_ghostCellCenter[idir], idir);
          m_physPowFab.divide(physicalCellWidth[idir], idir);
        }
      // powers of displacement from ghostCellCenter to
      // integration points on this ghost cell
      // m_physPowFab(m_integralCellBox, SpaceDim*(m_degree + 1));
      FORT_SETPOWERS2PLUS(CHF_FRA(m_physPowFab),
                          CHF_BOX(m_integralCellBox),
                          CHF_CONST_INT(m_degree));

      Box shiftedIntegralCellBox = m_integralCellBox + shiftCell;
      m_physPowFab.shift(+shiftCell);

      Real volCell = a_cellVolumesGhostedFabPtr->operator()(a_ghostCell, 0);
      FORT_SETVECTORAVGMOMENTS1PLUS(CHF_VR(a_ghostMoments),
                                    CHF_CONST_FRA(m_physPowFab),
                                    CHF_CONST_FRA1((*a_weightedJacobianGhostedFabPtr), 0),
                                    CHF_CONST_REAL(volCell),
                                    CHF_BOX(shiftedIntegralCellBox),
                                    CHF_CONST_INT(m_degree),
                                    CHF_CONST_INT(m_nvars),
                                    CHF_BOX(m_degreeBox));
      // return indices of m_physPowFab to original
      m_physPowFab.shift(-shiftCell);
    }
}


// ---------------------------------------------------------
bool MultiBlockLevelData::forbiddenZone(const IntVect&   a_iv,
                                        const Box&       a_bx)
{
  int numOut = 0;
  for (int idir = 0; idir < SpaceDim; idir++)
    {
      if (a_iv[idir] < a_bx.smallEnd(idir) ||
          a_iv[idir] > a_bx.bigEnd(idir))
        numOut++;
    }
  bool forbidden = (numOut > 1);
  return forbidden;
}

#include "NamespaceFooter.H"
