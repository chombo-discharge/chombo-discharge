/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBRedistribution.cpp
  @brief  Implementation of CD_EBRedistribution.cpp
  @author Robert Marskar
*/

// Chombo includes
#include <CH_Timer.H>
#include <EBCellFactory.H>
#include <NeighborIterator.H>
#include <EBFastFR.cpp>

// Our includes
#include <CD_EBRedistribution.H>
#include <CD_VofUtils.H>
#include <CD_BoxLoops.H>
#include <CD_EBAddOp.H>
#include <CD_NamespaceHeader.H>

EBRedistribution::EBRedistribution() noexcept
{
  CH_TIME("EBRedistribution::EBRedistribution(weak)");

  m_isDefined = false;
}

EBRedistribution::EBRedistribution(const EBLevelGrid& a_eblgCoar,
                                   const EBLevelGrid& a_eblgCoarsened,
                                   const EBLevelGrid& a_eblg,
                                   const EBLevelGrid& a_eblgRefined,
                                   const EBLevelGrid& a_eblgFine,
                                   const int          a_refToCoar,
                                   const int          a_refToFine,
                                   const bool         a_redistributeOutside) noexcept
{
  CH_TIME("EBRedistribution::EBRedistribution(full)");

  this->define(a_eblgCoar,
               a_eblgCoarsened,
               a_eblg,
               a_eblgRefined,
               a_eblgFine,
               a_refToCoar,
               a_refToFine,
               a_redistributeOutside);
}

EBRedistribution::~EBRedistribution() noexcept { CH_TIME("EBRedistribution::~EBRedistribution"); }

void
EBRedistribution::define(const EBLevelGrid& a_eblgCoar,
                         const EBLevelGrid& a_eblgCoarsened,
                         const EBLevelGrid& a_eblg,
                         const EBLevelGrid& a_eblgRefined,
                         const EBLevelGrid& a_eblgFine,
                         const int          a_refToCoar,
                         const int          a_refToFine,
                         const bool         a_redistributeOutside) noexcept
{
  CH_TIME("EBRedistribution::define");

  CH_assert(a_eblg.isDefined());

  m_redistRadius        = 1;
  m_refToCoar           = -1;
  m_refToFine           = -1;
  m_hasCoar             = false;
  m_hasFine             = false;
  m_redistributeOutside = false; //a_redistributeOutside;

  if (a_eblgCoar.isDefined()) {
    CH_assert(a_refToCoar >= 2);
    CH_assert(a_refToCoar % 2 == 0);
    CH_assert(a_eblgCoarsened.isDefined());

    m_eblgCoar      = a_eblgCoar;
    m_eblgCoarsened = a_eblgCoarsened;

    m_refToCoar = a_refToCoar;
    m_hasCoar   = true;
  }

  m_eblg = a_eblg;

  if (a_eblgFine.isDefined()) {
    CH_assert(a_refToFine >= 2);
    CH_assert(a_refToFine % 2 == 0);
    CH_assert(a_eblgFine.isDefined());

    m_eblgFine  = a_eblgFine;
    m_refToFine = a_refToFine;
    m_hasFine   = true;

    // Gods how I hate filling EBISLayouts. But in this case we need it because if the user has asked for
    // a larger refinement ratio than we have ghost cells, we need to explicity fetch the necessary geometric
    // data in a larger radius than what is available in the input arguments.
    if (m_refToFine > a_eblgRefined.getGhost()) {
      DisjointBoxLayout dblRefined;

      refine(dblRefined, m_eblg.getDBL(), m_refToFine);
      m_eblgRefined.define(dblRefined, m_eblgFine.getDomain(), m_refToFine, m_eblg.getEBIS());
    }
    else {
      CH_assert(a_eblgRefined.isDefined());

      m_eblgRefined = a_eblgRefined;
    }
  }

  this->defineStencils();
  this->defineBuffers();    


  m_isDefined = true;
}

void
EBRedistribution::defineStencils() noexcept
{
  CH_TIMERS("EBRedistribution::defineStencils");

  // TLDR: This is a bit involved since we need to define stencils on the valid cut-cells on this level. If there's an EBCF interface we
  //       need to know about it because we don't want to:
  //       1) Redistribute from this level into ghost cells on the other side of the coarse-fine interface. That mass should go on the coarse level.
  //       2) Redistribute from this level into regions covered by the finer grid. That mass should go on the fine level instead.

  const DisjointBoxLayout& dbl   = m_eblg.getDBL();
  const EBISLayout&        ebisl = m_eblg.getEBISL();

  const ProblemDomain& domain     = m_eblg.getDomain();
  const ProblemDomain  domainCoar = m_hasCoar ? m_eblgCoar.getDomain() : ProblemDomain();
  const ProblemDomain  domainFine = m_hasFine ? m_eblgFine.getDomain() : ProblemDomain();

  const Real dx     = 1.0;
  const Real dxCoar = dx * m_refToCoar;
  const Real dxFine = dx / m_refToFine;

  const Real vol     = std::pow(dx, SpaceDim);
  const Real volCoar = std::pow(dxCoar, SpaceDim);
  const Real volFine = std::pow(dxFine, SpaceDim);

  // These are maps of the valid cells on this level, and cells that lie on the interface. The interfaceCells data is used to figure out
  // which cells we will redistribute to when we redistribute from a cut-cell and across the coarse-fine interface into a coarse-grid cell. The
  // validCells is used for 1) restricting which cells we redistribute from, and 2) which fine-grid cells redistribute to.
  LevelData<BaseFab<bool>> interfaceCellsLD;
  LevelData<BaseFab<bool>> validCellsLD;

  this->defineValidCells(validCellsLD);
  this->defineInterfaceCells(interfaceCellsLD);

  m_vofit.define(dbl);
  m_redistStencilsCoar.define(dbl);
  m_redistStencilsLevel.define(dbl);
  m_redistStencilsFine.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box            box            = dbl[dit()];
    const EBISBox&       ebisBox        = ebisl[dit()];
    const EBGraph&       ebgraph        = ebisBox.getEBGraph();
    const IntVectSet     irregIVS       = ebisBox.getIrregIVS(box);
    const BaseFab<bool>& validCells     = validCellsLD[dit()];
    const BaseFab<bool>& interfaceCells = interfaceCellsLD[dit()];

    EBISBox ebisBoxCoar;
    EBISBox ebisBoxFine;

    if (m_hasCoar) {
      ebisBoxCoar = m_eblgCoarsened.getEBISL()[dit()];
    }
    if (m_hasFine) {
      ebisBoxFine = m_eblgRefined.getEBISL()[dit()];
    }

    IntVectSet redistCells;
    for (IVSIterator ivsIt(irregIVS); ivsIt.ok(); ++ivsIt) {
      const IntVect& iv = ivsIt();
      if (validCells(iv)) {
        redistCells |= iv;
      }
    }

    VoFIterator&           vofit         = m_vofit[dit()];
    BaseIVFAB<VoFStencil>& stencilsCoar  = m_redistStencilsCoar[dit()];
    BaseIVFAB<VoFStencil>& stencilsLevel = m_redistStencilsLevel[dit()];
    BaseIVFAB<VoFStencil>& stencilsFine  = m_redistStencilsFine[dit()];

    vofit.define(redistCells, ebgraph);
    stencilsCoar.define(redistCells, ebgraph, 1);
    stencilsLevel.define(redistCells, ebgraph, 1);
    stencilsFine.define(redistCells, ebgraph, 1);

    for (vofit.reset(); vofit.ok(); ++vofit) {
      const VolIndex& vof = vofit();

      // Get the VoFs in the neighborhood of this VoF. When we redistribute, also permit self-redistribution back into this cell.
      const bool             includeSelf = true;
      const Vector<VolIndex> neighborVoFs =
        VofUtils::getVofsInRadius(vof, ebisBox, m_redistRadius, VofUtils::Connectivity::SimplyConnected, includeSelf);

      Vector<VolIndex> fineVoFs;
      Vector<VolIndex> levelVoFs;
      Vector<VolIndex> coarVoFs;

      Real totalVolume = 0.0;

      // Go through the cells in the redistribution neighborhood of the current cut-cell. If the cell lies on
      // the coarse-side of the interface, we fetch the corresponding coarse cell. If the cell lie underneath
      // the finer level, we refine the cell.
      for (int i = 0; i < neighborVoFs.size(); i++) {
        const VolIndex& curVoF = neighborVoFs[i];
        const IntVect&  curIV  = curVoF.gridIndex();

        // If we only redistribute to cells inside the domain, skip this cell.
        if (!(domain.contains(curIV)) && !m_redistributeOutside) {
          continue;
        }

        if (!(validCells(curIV))) {
          CH_assert(m_hasFine);

          const Vector<VolIndex>& refinedVoFs = ebisl.refine(curVoF, m_refToFine, dit());

          fineVoFs.append(refinedVoFs);

          for (int i = 0; i < refinedVoFs.size(); i++) {
            totalVolume += volFine * ebisBoxFine.volFrac(refinedVoFs[i]);
          }

          continue;
        }

        if (interfaceCells(curIV)) {
          CH_assert(m_hasCoar);

          const VolIndex& coarsenedVoF = ebisl.coarsen(curVoF, m_refToCoar, dit());

          coarVoFs.push_back(coarsenedVoF);

          totalVolume += volCoar * ebisBoxCoar.volFrac(coarsenedVoF);

          continue;
        }

        levelVoFs.push_back(curVoF);

        totalVolume += vol * ebisBox.volFrac(curVoF);
      }

      VoFStencil& coarStencil  = stencilsCoar(vof, 0);
      VoFStencil& levelStencil = stencilsLevel(vof, 0);
      VoFStencil& fineStencil  = stencilsFine(vof, 0);

      const Real inverseVolume = 1.0 / totalVolume;

      for (int i = 0; i < coarVoFs.size(); i++) {
        const VolIndex& coarVoF = coarVoFs[i];

        if (domainCoar.contains(coarVoF.gridIndex())) {
          coarStencil.add(coarVoF, inverseVolume);
        }
      }

      for (int i = 0; i < levelVoFs.size(); i++) {
        const VolIndex& levelVoF = levelVoFs[i];

        if (domain.contains(levelVoF.gridIndex())) {
          levelStencil.add(levelVoF, inverseVolume);
        }
      }

      for (int i = 0; i < fineVoFs.size(); i++) {
        const VolIndex& fineVoF = fineVoFs[i];

        if (domainFine.contains(fineVoF.gridIndex())) {
          fineStencil.add(fineVoF, inverseVolume);
        }
      }
    }
  }
}

void
EBRedistribution::defineValidCells(LevelData<BaseFab<bool>>& a_validCells) const noexcept
{
  CH_TIME("EBRedistribution::defineValidCells");

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  a_validCells.define(dbl, 1, m_redistRadius * IntVect::Unit);
  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    a_validCells[dit()].setVal(true);
  }

  // If there's a finer level we need to figure out which cells on the fine level overlap with this level. Then we set those cells
  // to 'false'. If there's no finer level then all cells on this level are valid cells.
  if (m_hasFine) {
    DisjointBoxLayout dblCoFi;
    coarsen(dblCoFi, m_eblgFine.getDBL(), m_refToFine);

    // Create some data = 0 on the coarse grid and = 1 on the fine grid.
    LevelData<FArrayBox> data(dbl, 1, m_redistRadius * IntVect::Unit);
    LevelData<FArrayBox> dataCoFi(dblCoFi, 1, m_redistRadius * IntVect::Unit);

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      data[dit()].setVal(0.0);
    }
    for (DataIterator dit(dblCoFi); dit.ok(); ++dit) {
      dataCoFi[dit()].setVal(1.0);
    }

    // Need a new Copier here. 
    Copier copier(dblCoFi, dbl);
    dataCoFi.copyTo(Interval(0,0), data, Interval(0,0), copier);

    // Go through the coarse grid and set cells to false wherever we find data > 0.0
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box cellBox = dbl[dit()];

      BaseFab<bool>&   validCells = a_validCells[dit()];
      const FArrayBox& mask       = data[dit()];

      auto kernel = [&](const IntVect& iv) -> void {
        if (mask(iv) > 0.0) {
          validCells(iv) = false;
        }
      };

      BoxLoops::loop(mask.box(), kernel);
    }
  }
}

void
EBRedistribution::defineInterfaceCells(LevelData<BaseFab<bool>>& a_interfaceCells) const noexcept
{
  CH_TIME("EBRedistribution::defineInterfaceCells");

  const DisjointBoxLayout& dbl = m_eblg.getDBL();
  a_interfaceCells.define(dbl, 1, m_redistRadius * IntVect::Unit);

  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    a_interfaceCells[dit()].setVal(false);
  }

  // If there's a coarse grid then some of the cells we would redistribute to overlap with the valid region
  // on the coarse grid.
  if (m_hasCoar) {
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box grownBox = grow(dbl[dit()], m_redistRadius);

      DenseIntVectSet divs(grownBox, true);

      NeighborIterator nit(dbl);
      for (nit.begin(dit()); nit.ok(); ++nit) {
        divs -= dbl[nit()];
      }
      divs -= dbl[dit()];

      // Flag the cells on the other side of the interface.
      BaseFab<bool>& interfaceCells = a_interfaceCells[dit()];
      for (DenseIntVectSetIterator ivsIt(divs); ivsIt.ok(); ++ivsIt) {
        interfaceCells(ivsIt()) = true;
      }
    }
  }
}

void
EBRedistribution::defineBuffers() noexcept
{
  CH_TIME("EBRedistribution::defineBuffers");

  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const ProblemDomain&     domain = m_eblg.getDomain();

  // ghostDefine implies that we copy from valid+ghost to valid. Which is what we want since we
  // also redistribute over patch boundaries (duh). Note that the ghost cells we define here MUST
  // match the ghost cells in the buffers that are defined in the actual redistribution routines.

  if (m_hasCoar) {
    const DisjointBoxLayout& dblCoar      = m_eblgCoar.getDBL();
    const DisjointBoxLayout& dblCoarsened = m_eblgCoarsened.getDBL();
    const ProblemDomain&     domainCoar   = m_eblgCoar.getDomain();

    m_coarCopier.ghostDefine(dblCoarsened, dblCoar, domainCoar, m_redistRadius * IntVect::Unit);
  }

  m_levelCopier.ghostDefine(dbl, dbl, domain, m_redistRadius * IntVect::Unit);

  if (m_hasFine) {
    const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
    const DisjointBoxLayout& dblRefined = m_eblgRefined.getDBL();
    const ProblemDomain&     domainFine = m_eblgFine.getDomain();

    m_fineCopier.ghostDefine(dblRefined, dblFine, domainFine, m_refToFine * m_redistRadius * IntVect::Unit);
  }
}

void
EBRedistribution::redistributeAMR(LevelData<EBCellFAB>*             a_phiCoar,
                                  LevelData<EBCellFAB>*             a_phi,
                                  LevelData<EBCellFAB>*             a_phiFine,
                                  const LevelData<BaseIVFAB<Real>>& a_deltaM,
                                  const Real                        a_scaleCoar,
                                  const Real                        a_scale,
                                  const Real                        a_scaleFine,
                                  const Interval&                   a_variables) const noexcept
{
  CH_TIME("EBRedistribution::redistributeAMR");

  CH_assert(m_isDefined);
  CH_assert(a_phi != nullptr);

  if (m_hasCoar) {
    CH_assert(a_phiCoar != nullptr);

    this->redistributeCoar(*a_phiCoar, a_deltaM, a_scaleCoar, a_variables);
  }

  this->redistributeLevel(*a_phi, a_deltaM, a_scale, a_variables);

  if (m_hasFine) {
    CH_assert(a_phiFine != nullptr);

    this->redistributeFine(*a_phiFine, a_deltaM, a_scaleFine, a_variables);
  }
}

void
EBRedistribution::redistributeCoar(LevelData<EBCellFAB>&             a_phiCoar,
                                   const LevelData<BaseIVFAB<Real>>& a_deltaM,
                                   const Real&                       a_scaleCoar,
                                   const Interval&                   a_variables) const noexcept
{
  CH_TIMERS("EBRedistribution::redistributeCoar");
  CH_TIMER("EBRedistribution::redistributeCoar::irreg_cells", t1);

  CH_assert(m_isDefined);
  CH_assert(m_hasCoar);
  CH_assert(a_phiCoar.isDefined());
  CH_assert(a_deltaM.isDefined());
  CH_assert(a_deltaM.disjointBoxLayout() == m_eblg.getDBL());
  CH_assert(a_phiCoar.disjointBoxLayout() == m_eblgCoar.getDBL());
  CH_assert(a_phiCoar.nComp() > a_variables.end());
  CH_assert(a_deltaM.nComp() > a_variables.end());

  const DisjointBoxLayout& dblCoar    = m_eblgCoarsened.getDBL();
  const EBISLayout&        ebislCoar  = m_eblgCoarsened.getEBISL();
  const ProblemDomain&     domainCoar = m_eblgCoarsened.getDomain();

  LevelData<EBCellFAB> coarBuffer(dblCoar, 1, m_redistRadius * IntVect::Unit, EBCellFactory(ebislCoar));

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {
    for (DataIterator dit(dblCoar); dit.ok(); ++dit) {
      const BaseIVFAB<Real>&       deltaM   = a_deltaM[dit()];
      const BaseIVFAB<VoFStencil>& stencils = m_redistStencilsCoar[dit()];

      EBCellFAB& buffer = coarBuffer[dit()];
      buffer.setVal(0.0);

      // Apply stencil into buffer.
      CH_START(t1);
      VoFIterator& vofit = m_vofit[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit) {
        const VolIndex&   vof      = vofit();
        const VoFStencil& stencil  = stencils(vof, 0);
        const Real&       massDiff = deltaM(vof, ivar);

        for (int i = 0; i < stencil.size(); i++) {
          buffer(stencil.vof(i), 0) += a_scaleCoar * stencil.weight(i) * massDiff;
        }
      }
      CH_STOP(t1);
    }

    // Increment a_phi by the result.
    const Interval srcInterv = Interval(0, 0);
    const Interval dstInterv = Interval(ivar, ivar);

    coarBuffer.copyTo(srcInterv, a_phiCoar, dstInterv, m_coarCopier, EBAddOp());
  }
}
void
EBRedistribution::redistributeLevel(LevelData<EBCellFAB>&             a_phi,
                                    const LevelData<BaseIVFAB<Real>>& a_deltaM,
                                    const Real&                       a_scale,
                                    const Interval&                   a_variables) const noexcept
{
  CH_TIMERS("EBRedistribution::redistributeLevel");
  CH_TIMER("EBRedistribution::redistributeLevel::irreg_cells", t1);

  CH_assert(m_isDefined);
  CH_assert(a_phi.isDefined());
  CH_assert(a_deltaM.isDefined());
  CH_assert(a_deltaM.disjointBoxLayout() == m_eblg.getDBL());
  CH_assert(a_phi.disjointBoxLayout() == m_eblg.getDBL());
  CH_assert(a_phi.nComp() > a_variables.end());
  CH_assert(a_deltaM.nComp() > a_variables.end());

  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const EBISLayout&        ebisl  = m_eblg.getEBISL();
  const ProblemDomain&     domain = m_eblg.getDomain();

  LevelData<EBCellFAB> levelBuffer(dbl, 1, m_redistRadius * IntVect::Unit, EBCellFactory(ebisl));

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const BaseIVFAB<Real>&       deltaM   = a_deltaM[dit()];
      const BaseIVFAB<VoFStencil>& stencils = m_redistStencilsLevel[dit()];

      EBCellFAB& buffer = levelBuffer[dit()];
      buffer.setVal(0.0);

      // Apply stencil into buffer.
      CH_START(t1);
      VoFIterator& vofit = m_vofit[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit) {
        const VolIndex&   vof      = vofit();
        const VoFStencil& stencil  = stencils(vof, 0);
        const Real&       massDiff = deltaM(vof, ivar);

        for (int i = 0; i < stencil.size(); i++) {
          buffer(stencil.vof(i), 0) += a_scale * stencil.weight(i) * massDiff;
        }
      }
      CH_STOP(t1);
    }

    // Increment a_phi by the result.
    const Interval srcInterv = Interval(0, 0);
    const Interval dstInterv = Interval(ivar, ivar);

    levelBuffer.copyTo(srcInterv, a_phi, dstInterv, m_levelCopier, EBAddOp());
  }
}

void
EBRedistribution::redistributeFine(LevelData<EBCellFAB>&             a_phiFine,
                                   const LevelData<BaseIVFAB<Real>>& a_deltaM,
                                   const Real&                       a_scaleFine,
                                   const Interval&                   a_variables) const noexcept
{
  CH_TIMERS("EBRedistribution::redistributeFine");
  CH_TIMER("EBRedistribution::redistributeFine::irreg_cells", t1);

  CH_assert(m_isDefined);
  CH_assert(m_hasFine);
  CH_assert(a_phiFine.isDefined());
  CH_assert(a_deltaM.isDefined());
  CH_assert(a_deltaM.disjointBoxLayout() == m_eblg.getDBL());
  CH_assert(a_phiFine.disjointBoxLayout() == m_eblgFine.getDBL());
  CH_assert(a_phiFine.nComp() > a_variables.end());
  CH_assert(a_deltaM.nComp() > a_variables.end());

  const DisjointBoxLayout& dblFine   = m_eblgRefined.getDBL();
  const EBISLayout&        ebislFine = m_eblgRefined.getEBISL();

  LevelData<EBCellFAB> fineBuffer(dblFine, 1, m_refToFine * m_redistRadius * IntVect::Unit, EBCellFactory(ebislFine));

  for (int ivar = a_variables.begin(); ivar <= a_variables.end(); ivar++) {
    for (DataIterator dit(dblFine); dit.ok(); ++dit) {
      const BaseIVFAB<Real>&       deltaM   = a_deltaM[dit()];
      const BaseIVFAB<VoFStencil>& stencils = m_redistStencilsFine[dit()];

      EBCellFAB& buffer = fineBuffer[dit()];
      buffer.setVal(0.0);

      // Apply stencil into buffer.
      CH_START(t1);
      VoFIterator& vofit = m_vofit[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit) {
        const VolIndex&   vof      = vofit();
        const VoFStencil& stencil  = stencils(vof, 0);
        const Real&       massDiff = deltaM(vof, ivar);

        for (int i = 0; i < stencil.size(); i++) {
          buffer(stencil.vof(i), 0) += a_scaleFine * stencil.weight(i) * massDiff;
        }
      }
      CH_STOP(t1);
    }

    // Increment a_phi by the result.
    const Interval srcInterv = Interval(0, 0);
    const Interval dstInterv = Interval(ivar, ivar);

    fineBuffer.copyTo(srcInterv, a_phiFine, dstInterv, m_fineCopier, EBAddOp());
  }
}

#include <CD_NamespaceFooter.H>
