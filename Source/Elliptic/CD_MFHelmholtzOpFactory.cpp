/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_MFHelmholtzOpFactory.cpp
 * @brief  Implementation of CD_MFHelmholtzOpFactory.H
 * @author Robert Marskar
 */

// Std includes
#include <sstream>

// Chombo includes
#include <ParmParse.H>
#include <BaseIVFactory.H>
#include <CH_Timer.H>

// Our includes
#include <CD_MFHelmholtzOpFactory.H>
#include <CD_MultifluidAlias.H>
#include <CD_DataOps.H>
#include <CD_MFBaseIVFAB.H>
#include <CD_NamespaceHeader.H>

constexpr int MFHelmholtzOpFactory::m_comp;
constexpr int MFHelmholtzOpFactory::m_nComp;
constexpr int MFHelmholtzOpFactory::m_mainPhase;

MFHelmholtzOpFactory::MFHelmholtzOpFactory(const MFIS&                  a_mfis,
                                           const Location::Cell         a_dataLocation,
                                           const Real&                  a_alpha,
                                           const Real&                  a_beta,
                                           const RealVect&              a_probLo,
                                           const AmrLevelGrids&         a_amrLevelGrids,
                                           const AmrMask&               a_validCells,
                                           const AmrInterpolators&      a_amrInterpolators,
                                           const AmrFluxRegisters&      a_amrFluxRegisters,
                                           const AmrCoarseners&         a_amrCoarseners,
                                           const AmrRefRatios&          a_amrRefRatios,
                                           const AmrResolutions&        a_amrResolutions,
                                           const AmrCellData&           a_amrAcoef,
                                           const AmrFluxData&           a_amrBcoef,
                                           const AmrIrreData&           a_amrBcoefIrreg,
                                           const DomainBCFactory&       a_domainBcFactory,
                                           const EBBCFactory&           a_ebBcFactory,
                                           const JumpBCFactory&         a_jumpBcFactory,
                                           const IntVect&               a_ghostPhi,
                                           const IntVect&               a_ghostRhs,
                                           const Smoother&              a_smoother,
                                           const Real&                  a_relaxFactor,
                                           const int&                   a_chebyOrder,
                                           const Real&                  a_chebyEigRatio,
                                           const int&                   a_rasInnerSweeps,
                                           const int&                   a_jumpOrder,
                                           const int&                   a_jumpWeight,
                                           const int&                   a_preCondSmooth,
                                           const Vector<AmrLevelGrids>& a_mgLevelGrids,
                                           const Vector<AmrLevelGrids>& a_mgLevelGridsCoFi,
                                           const Vector<AmrCoarseners>& a_mgCoarseners,
                                           const AmrLevelGrids&         a_amrLevelGridsCoFi,
                                           const ProblemDomain&         a_bottomDomain,
                                           const bool                   a_refluxFree)
  : m_mfis(a_mfis),
    m_dataLocation(a_dataLocation),
    m_smoother(a_smoother),
    m_chebyOrder(a_chebyOrder),
    m_chebyEigRatio(a_chebyEigRatio),
    m_rasInnerSweeps(a_rasInnerSweeps),
    m_numPreCondSmooth(a_preCondSmooth),
    m_ghostPhi(a_ghostPhi),
    m_ghostRhs(a_ghostRhs),
    m_alpha(a_alpha),
    m_beta(a_beta),
    m_relaxFactor(a_relaxFactor),
    m_refluxFree(a_refluxFree),
    m_probLo(a_probLo),
    m_amrLevelGrids(a_amrLevelGrids),
    m_validCells(a_validCells),
    m_amrInterpolators(a_amrInterpolators),
    m_amrFluxRegisters(a_amrFluxRegisters),
    m_amrCoarseners(a_amrCoarseners),
    m_amrRefRatios(a_amrRefRatios),
    m_amrResolutions(a_amrResolutions),
    m_amrAcoef(a_amrAcoef),
    m_amrBcoef(a_amrBcoef),
    m_amrBcoefIrreg(a_amrBcoefIrreg),
    m_domainBcFactory(a_domainBcFactory),
    m_ebBcFactory(a_ebBcFactory),
    m_jumpBcFactory(a_jumpBcFactory),
    m_bottomDomain(a_bottomDomain),
    m_amrLevelGridsCoFi(a_amrLevelGridsCoFi),
    m_jumpOrder(a_jumpOrder),
    m_jumpWeight(a_jumpWeight)
{
  CH_TIME("MFHelmholtzOpFactory::MFHelmholtzOpFactory()");

  m_numAmrLevels = static_cast<int>(m_amrLevelGrids.size());

  // Asking multigrid to do the bottom solve at a refined AMR level is classified as bad input.
  if (ChomboDischarge::MFHelmholtzOpFactory::isFiner(m_bottomDomain, m_amrLevelGrids[0].getDomain())) {
    MayDay::Abort("MFHelmholtzOpFactory -- bottomsolver domain can't be larger than the base AMR domain!");
  }

  // Define the jump data and the multigrid levels.
  this->defineJump();
  this->defineMultigridLevels(a_mgLevelGrids, a_mgLevelGridsCoFi, a_mgCoarseners);
  this->setJump(0.0, 0.0);
}

MFHelmholtzOpFactory::~MFHelmholtzOpFactory()
{
  CH_TIME("MFHelmholtzOpFactory::~MFHelmholtzOpFactory()");
}

void
MFHelmholtzOpFactory::setJump(const EBAMRIVData& a_sigma, const Real& a_scale)
{
  CH_TIME("MFHelmholtzOpFactory::setJump(EBAMRIVData, Real)");

  // TLDR: This routine sets the jump coefficient on each level. Data on deeper levels
  //       are coarsened from finer ones.

  const Interval interv(m_comp, m_comp);

  for (int lvl = 0; lvl < m_numAmrLevels; lvl++) {
    const DisjointBoxLayout& dbl = a_sigma[lvl]->disjointBoxLayout();
    const DataIterator&      dit = dbl.dataIterator();

    const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box      box = dbl[din];
      const Interval interv(m_comp, m_comp);

      (*m_amrJump[lvl])[din].copy(box, interv, box, (*a_sigma[lvl])[din], interv);
    }

    LayoutData<VoFIterator> vofIter(dbl);
#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex&       din = dit[mybox];
      const BaseIVFAB<Real>& ivf = (*m_amrJump[lvl])[din];
      vofIter[din]               = VoFIterator(ivf.getIVS(), ivf.getEBGraph());
    }

    DataOps::scale(*m_amrJump[lvl], a_scale, vofIter);

    m_amrJump[lvl]->exchange();
  }

  // Average down on AMR levels.
  for (int lvl = m_numAmrLevels - 1; lvl > 0; lvl--) {
    const RefCountedPtr<EBCoarAve>& aveOp = m_amrCoarseners[lvl].getAveOp(m_mainPhase);

    const Average average = Average::Arithmetic;

    aveOp->averageData(*m_amrJump[lvl - 1], *m_amrJump[lvl], Interval(m_comp, m_comp), average);
  }
}

void
MFHelmholtzOpFactory::setJump(const Real& a_sigma, const Real& a_scale)
{
  CH_TIME("MFHelmholtzOpFactory::setJump(Real, Real)");

  // TLDR: This setse the jump coefficient on each level.
  DataOps::setValue(m_amrJump, a_sigma * a_scale);

  for (int i = 0; i < m_mgJump.size(); i++) {
    DataOps::setValue(m_mgJump[i], a_sigma * a_scale);
  }
}

const EBAMRIVData&
MFHelmholtzOpFactory::getSigma() const
{
  CH_TIME("MFHelmholtzOpFactory::getSigma()");

  return m_amrJump;
}

void
MFHelmholtzOpFactory::defineJump()
{
  CH_TIME("MFHelmholtzOpFactory::defineJump()");

  // clang-format off
  // TLDR: This defines m_amrJump on the first phase (gas phase). This is irregular data intended to be interfaced into the
  //       boundary condition class. When we match the BC we get the data from here (the gas phase). Note that we define
  //       m_amrJump on all irregular cells, but the operators will do matching on a subset of them.
  // clang-format on

  m_amrJump.resize(m_numAmrLevels);

  for (int lvl = 0; lvl < m_numAmrLevels; lvl++) {
    const DisjointBoxLayout& dbl   = m_amrLevelGrids[lvl].getGrids();
    const EBLevelGrid&       eblg  = m_amrLevelGrids[lvl].getEBLevelGrid(m_mainPhase);
    const EBISLayout&        ebisl = eblg.getEBISL();
    const DataIterator&      dit   = dbl.dataIterator();

    const int nbox = dit.size();

    LayoutData<IntVectSet> irregCells(dbl);

#pragma omp parallel for schedule(runtime)
    for (int mybox = 0; mybox < nbox; mybox++) {
      const DataIndex& din = dit[mybox];

      const Box      box     = dbl[din];
      const EBISBox& ebisbox = ebisl[din];

      irregCells[din] = ebisbox.getIrregIVS(box);
    }

    BaseIVFactory<Real> fact(ebisl, irregCells);
    m_amrJump[lvl] = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
      new LevelData<BaseIVFAB<Real>>(dbl, m_nComp, IntVect::Zero, fact));
  }
}

void
MFHelmholtzOpFactory::defineMultigridLevels(const Vector<AmrLevelGrids>& a_mgLevelGrids,
                                            const Vector<AmrLevelGrids>& a_mgLevelGridsCoFi,
                                            const Vector<AmrCoarseners>& a_mgCoarseners)
{
  CH_TIME("MFHelmholtzOpFactory::defineMultigridLevels(Vector<AmrLevelGrids>, Vector<AmrLevelGrids>, "
          "Vector<AmrCoarseners>)");

  // TLDR: This routine defines what is needed for making the multigrid levels. This includes the intermediate
  // levels (if you run with refinement factor 4) as well as the deeper multigrid levels that are coarsenings of the
  // base AMR level. Recall that in Chombo-speak a multigrid level is a grid level completely covered by another grid
  // level. The realm builds the grids; this factory takes the entries that are not coarser than its bottom domain and
  // defines the coefficients on them.

  CH_assert(static_cast<int>(a_mgLevelGrids.size()) == m_numAmrLevels);
  CH_assert(static_cast<int>(a_mgLevelGridsCoFi.size()) == m_numAmrLevels);
  CH_assert(static_cast<int>(a_mgCoarseners.size()) == m_numAmrLevels);

  m_mgLevelGrids.resize(m_numAmrLevels);
  m_mgAcoef.resize(m_numAmrLevels);
  m_mgBcoef.resize(m_numAmrLevels);
  m_mgBcoefIrreg.resize(m_numAmrLevels);
  m_mgJump.resize(m_numAmrLevels);
  m_mgCoarseners.resize(m_numAmrLevels);
  m_mgLevelGridsCoFi.resize(m_numAmrLevels);
  m_hasMgLevels.resize(m_numAmrLevels);

  for (int amrLevel = 0; amrLevel < m_numAmrLevels; amrLevel++) {
    const AmrLevelGrids& realmLevels = a_mgLevelGrids[amrLevel];

    CH_assert(realmLevels.size() >= 1);
    CH_assert(realmLevels[0].getDomain() == m_amrLevelGrids[amrLevel].getDomain());

    // Use the realm's levels down to, and including, the last one that is not coarser than the bottom domain.
    int numMgLevels = 1;

    while (numMgLevels < static_cast<int>(realmLevels.size()) &&
           !isCoarser(realmLevels[numMgLevels].getDomain(), m_bottomDomain)) {
      numMgLevels++;
    }

    m_hasMgLevels[amrLevel] = numMgLevels > 1;

    // The realm builds down to the coarsest bottom domain registered on it, so it stops before this factory's bottom
    // domain only if the grids could not be coarsened further, or if the bottom cell count was lowered after
    // registration.
    if (amrLevel == 0 && isCoarser(m_bottomDomain, realmLevels.back().getDomain())) {
      std::stringstream ss;

      ss << "MFHelmholtzOpFactory::defineMultigridLevels - the multigrid levels stop at "
         << realmLevels.back().getDomain() << ", which is finer than the bottom domain " << m_bottomDomain;

      MayDay::Warning(ss.str().c_str());
    }

    m_mgLevelGrids[amrLevel].resize(0);
    m_mgAcoef[amrLevel].resize(0);
    m_mgBcoef[amrLevel].resize(0);
    m_mgBcoefIrreg[amrLevel].resize(0);
    m_mgJump[amrLevel].resize(0);
    m_mgCoarseners[amrLevel].resize(0);
    m_mgLevelGridsCoFi[amrLevel].resize(0);

    if (m_hasMgLevels[amrLevel]) {
      m_mgLevelGrids[amrLevel].push_back(m_amrLevelGrids[amrLevel]);
      m_mgAcoef[amrLevel].push_back(m_amrAcoef[amrLevel]);
      m_mgBcoef[amrLevel].push_back(m_amrBcoef[amrLevel]);
      m_mgBcoefIrreg[amrLevel].push_back(m_amrBcoefIrreg[amrLevel]);
      m_mgJump[amrLevel].push_back(m_amrJump[amrLevel]);

      for (int img = 1; img < numMgLevels; img++) {
        const MFLevelGrid& mgMflgCoar = realmLevels[img];

        for (int i = 0; i < mgMflgCoar.numPhases(); i++) {
          CH_assert(mgMflgCoar.getEBLevelGrid(i).getGhost() >= m_ghostPhi.max());
        }

        const DisjointBoxLayout& dblCoar = mgMflgCoar.getGrids();

        // The Chombo multifluid things need Vector<EBISLayout> for the factories.
        Vector<int>        ebislComps;
        Vector<EBISLayout> ebislCoar;

        for (int i = 0; i < mgMflgCoar.numPhases(); i++) {
          ebislComps.push_back(m_nComp);
          ebislCoar.push_back(mgMflgCoar.getEBLevelGrid(i).getEBISL());
        }

        // Factories for making coarse stuff. Need one ghost cell because we interpolate b*grad(phi) to face
        // centroids.
        const int           nghost = 1;
        MFCellFactory       cellFact(ebislCoar, ebislComps);
        MFFluxFactory       fluxFact(ebislCoar, ebislComps);
        MFBaseIVFABFactory  ivFact(ebislCoar, ebislComps);
        BaseIVFactory<Real> irregFact(mgMflgCoar.getEBLevelGrid(m_mainPhase).getEBISL());

        // Multifluid a bit special -- number of components come in through the factory.
        const int dummy = 1;

        auto coarAcoef = RefCountedPtr<LevelData<MFCellFAB>>(
          new LevelData<MFCellFAB>(dblCoar, dummy, m_ghostPhi, cellFact));
        auto coarBcoef = RefCountedPtr<LevelData<MFFluxFAB>>(
          new LevelData<MFFluxFAB>(dblCoar, dummy, m_ghostPhi, fluxFact));
        auto coarBcoefIrreg = RefCountedPtr<LevelData<MFBaseIVFAB>>(
          new LevelData<MFBaseIVFAB>(dblCoar, dummy, m_ghostPhi, ivFact));
        auto coarJump = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
          new LevelData<BaseIVFAB<Real>>(dblCoar, m_nComp, nghost * IntVect::Zero, irregFact));

        DataOps::setValue(*coarJump, 0.0);

        m_mgLevelGrids[amrLevel].push_back(mgMflgCoar);
        m_mgAcoef[amrLevel].push_back(coarAcoef);
        m_mgBcoef[amrLevel].push_back(coarBcoef);
        m_mgBcoefIrreg[amrLevel].push_back(coarBcoefIrreg);
        m_mgJump[amrLevel].push_back(coarJump);
        m_mgCoarseners[amrLevel].push_back(a_mgCoarseners[amrLevel][img - 1]);
        m_mgLevelGridsCoFi[amrLevel].push_back(a_mgLevelGridsCoFi[amrLevel][img - 1]);
      }
    }
  }

  // Coarsen coefficients for the deeper multigrid levels.
  this->coarsenCoefficientsMG();
}

void
MFHelmholtzOpFactory::coarsenCoefficientsMG()
{
  CH_TIME("MFHelmholtzOpFactory::coarsenCoefficientsMG");

  for (int amrLevel = 0; amrLevel < m_numAmrLevels; amrLevel++) {

    if (m_hasMgLevels[amrLevel]) {
      // In these vectors, mgAco[0] is the AMR level, mgAco[1] is a refinement 2 coarsening of mgAco[0], mgAco[2] is the
      // coarsening of mgAco[1] and so on.
      const AmrLevelGrids mgGrids = m_mgLevelGrids[amrLevel];

      AmrCellData& mgAco      = m_mgAcoef[amrLevel];
      AmrFluxData& mgBco      = m_mgBcoef[amrLevel];
      AmrIrreData& mgBcoIrreg = m_mgBcoefIrreg[amrLevel];

      CH_assert(mgAco.size() == (size_t)mgGrids.size());
      CH_assert(mgBco.size() == (size_t)mgGrids.size());
      CH_assert(mgBcoIrreg.size() == (size_t)mgGrids.size());

      const Interval interv(m_comp, m_comp);
      const Average  average = Average::Arithmetic;

      for (int mgLevel = 0; mgLevel < mgGrids.size() - 1; mgLevel++) {
        const MFCoarAve& mfAveOp = m_mgCoarseners[amrLevel][mgLevel];

        // Average down on each phase.
        for (int i = 0; i < mgGrids[mgLevel].numPhases(); i++) {
          LevelData<EBCellFAB>       coarAco;
          LevelData<EBFluxFAB>       coarBco;
          LevelData<BaseIVFAB<Real>> coarBcoIrreg;

          LevelData<EBCellFAB>       fineAco;
          LevelData<EBFluxFAB>       fineBco;
          LevelData<BaseIVFAB<Real>> fineBcoIrreg;

          MultifluidAlias::aliasMF(coarAco, i, *mgAco[mgLevel + 1]);
          MultifluidAlias::aliasMF(coarBco, i, *mgBco[mgLevel + 1]);
          MultifluidAlias::aliasMF(coarBcoIrreg, i, *mgBcoIrreg[mgLevel + 1]);

          MultifluidAlias::aliasMF(fineAco, i, *mgAco[mgLevel]);
          MultifluidAlias::aliasMF(fineBco, i, *mgBco[mgLevel]);
          MultifluidAlias::aliasMF(fineBcoIrreg, i, *mgBcoIrreg[mgLevel]);

          const EBCoarAve& aveOp = *mfAveOp.getAveOp(i);

          aveOp.averageData(coarAco, fineAco, interv, average);
          aveOp.averageData(coarBco, fineBco, interv, average);
          aveOp.averageData(coarBcoIrreg, fineBcoIrreg, interv, average);

          coarAco.exchange();
          coarBco.exchange();
          coarBcoIrreg.exchange();
        }
      }
    }
  }
}

MFHelmholtzOp*
MFHelmholtzOpFactory::MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool /*a_homogeneousOnly*/)
{
  CH_TIME("EBHelmholtzOpFactory::MGnewOp(ProblemDomain, int, bool)");

  MFHelmholtzOp* mgOp = nullptr;

  // First, check if multigrid should be turned off completely.
  ParmParse pp;
  bool      turnOffMG = false;
  pp.query("turn_off_multigrid", turnOffMG);
  if (turnOffMG) {
    pout() << "turning off multigrid for MFHelmholtzOpFactory because 'turn_off_multigrid' = true" << endl;
    return mgOp;
  }

  const int amrLevel = this->findAmrLevel(
    a_fineDomain); // Run-time abort if a_fineDomain is not found in anhy amr level.

  const int mgRefRat = 2;

  // Things that are needed for defining the operator. This might seem weird but
  // the multigrid operators only do relaxation and there's no fine-to-coar stuff. So hasFine = hasCoar = false.
  // But we might need to restrict and interpolate another, even coarser multigrid level so in that case we need
  // to define that.
  MFLevelGrid mflg;
  MFLevelGrid mflgMgCoFi;
  MFLevelGrid mflgMgCoar;

  MFMultigridInterpolator interpolator; // Only if defined on an AMR level
  MFReflux                fluxReg;      // Only if defined on an AMR level
  MFCoarAve               coarsener;    // Only if defined on an AMR level

  bool hasMGObjects;

  RefCountedPtr<LevelData<MFCellFAB>>       Acoef;
  RefCountedPtr<LevelData<MFFluxFAB>>       Bcoef;
  RefCountedPtr<LevelData<MFBaseIVFAB>>     BcoefIrreg;
  RefCountedPtr<LevelData<BaseIVFAB<Real>>> jump;
  RefCountedPtr<LevelData<BaseFab<bool>>>   validCells(nullptr);

  bool foundMgLevel = false;

  if (a_depth == 0) { // Asking for the AMR level.
    mflg       = m_amrLevelGrids[amrLevel];
    Acoef      = m_amrAcoef[amrLevel];
    Bcoef      = m_amrBcoef[amrLevel];
    BcoefIrreg = m_amrBcoefIrreg[amrLevel];
    jump       = m_amrJump[amrLevel];
    validCells = m_validCells[amrLevel];

    interpolator = m_amrInterpolators[amrLevel];
    fluxReg      = m_amrFluxRegisters[amrLevel];
    coarsener    = m_amrCoarseners[amrLevel];

    hasMGObjects = m_hasMgLevels[amrLevel];

    if (hasMGObjects) {
      mflgMgCoFi = m_mgLevelGridsCoFi[amrLevel][0];
      mflgMgCoar = m_mgLevelGrids[amrLevel][1];
    }

    foundMgLevel = true;
  }
  else { // Asking for a coarsening. No interp or flux reg object here.
    // TLDR: Go through the coarsened levels for the specified amr level and see if we find a coarsening at the
    //       specified depth.
    const ProblemDomain coarDomain = coarsen(a_fineDomain, static_cast<int>(std::pow(mgRefRat, a_depth)));

    // These are the things that live below the AMR level corresponding to a_fineDomain.
    const AmrLevelGrids& mgLevelGrids = m_mgLevelGrids[amrLevel];
    const AmrCellData&   mgAcoef      = m_mgAcoef[amrLevel];
    const AmrFluxData&   mgBcoef      = m_mgBcoef[amrLevel];
    const AmrIrreData&   mgBcoefIrreg = m_mgBcoefIrreg[amrLevel];

    // See if we have a corresponding multigrid level.
    int mgLevel = 0;

    for (int img = 0; img < mgLevelGrids.size(); img++) {
      if (mgLevelGrids[img].getDomain() == coarDomain) {
        mgLevel      = img;
        foundMgLevel = true;

        break;
      }
    }

    // Found the multigrid level. We can define the operator.
    if (foundMgLevel) {
      Acoef      = mgAcoef[mgLevel];
      Bcoef      = mgBcoef[mgLevel];
      BcoefIrreg = mgBcoefIrreg[mgLevel];
      mflg       = mgLevelGrids[mgLevel];
      jump       = m_mgJump[amrLevel][mgLevel];

      hasMGObjects = (mgLevel < mgLevelGrids.size() - 1); // This just means that mgLevel was not the last entry in
                                                          // mgLevelGrids so there's even coarser stuff below.
      if (hasMGObjects) {
        mflgMgCoFi = m_mgLevelGridsCoFi[amrLevel][mgLevel];
        mflgMgCoar = mgLevelGrids[mgLevel + 1];
      }
    }
  }

  // Make the operator
  if (foundMgLevel) {
    const Real dx = m_amrResolutions[amrLevel] * std::pow(mgRefRat, a_depth); //

    constexpr int bogusRef = 2;

    mgOp = new MFHelmholtzOp(m_dataLocation,
                             MFLevelGrid(),
                             mflg,
                             MFLevelGrid(),
                             MFLevelGrid(),
                             mflgMgCoFi,
                             mflgMgCoar,
                             interpolator,
                             fluxReg,
                             coarsener,
                             validCells,
                             m_domainBcFactory,
                             m_ebBcFactory,
                             m_jumpBcFactory,
                             m_probLo,
                             dx,
                             bogusRef,
                             bogusRef,
                             false,
                             false,
                             hasMGObjects,
                             true,
                             m_alpha,
                             m_beta,
                             Acoef,
                             Bcoef,
                             BcoefIrreg,
                             m_ghostPhi,
                             m_ghostRhs,
                             m_jumpOrder,
                             m_jumpWeight,
                             m_numPreCondSmooth,
                             m_smoother,
                             m_relaxFactor,
                             m_chebyOrder,
                             m_chebyEigRatio,
                             m_rasInnerSweeps,
                             m_refluxFree);

    mgOp->setJump(jump);
  }

  return mgOp;
}

MFHelmholtzOp*
MFHelmholtzOpFactory::AMRnewOp(const ProblemDomain& a_domain)
{
  CH_TIME("MFHelmholtzOpFactory::AMRnewOp(ProblemDomain)");

  const int amrLevel = this->findAmrLevel(a_domain);

  const bool hasFine = amrLevel < m_numAmrLevels - 1;
  const bool hasCoar = amrLevel > 0;

  MFLevelGrid mflgFine;
  MFLevelGrid mflg;
  MFLevelGrid mflgCoFi;
  MFLevelGrid mflgCoar;
  MFLevelGrid mflgCoFiMG;
  MFLevelGrid mflgCoarMG;

  Real dx;

  mflg = m_amrLevelGrids[amrLevel];
  dx   = m_amrResolutions[amrLevel];

  int refToCoar = 1;
  int refToFine = 1;

  if (hasCoar) {
    mflgCoar  = m_amrLevelGrids[amrLevel - 1];
    refToCoar = m_amrRefRatios[amrLevel - 1];
  }

  if (hasFine) {
    refToFine = m_amrRefRatios[amrLevel];
    mflgFine  = m_amrLevelGrids[amrLevel + 1];
  }

  const bool hasMGObjects = m_hasMgLevels[amrLevel];
  if (hasMGObjects) {
    mflgCoFiMG = m_mgLevelGridsCoFi[amrLevel][0];
    mflgCoarMG = m_mgLevelGrids[amrLevel][1];
  }

  if (hasCoar) {
    mflgCoFi = m_amrLevelGridsCoFi[amrLevel - 1];
  }

  auto* op = new MFHelmholtzOp(m_dataLocation,
                               mflgFine,
                               mflg,
                               mflgCoFi,
                               mflgCoar,
                               mflgCoFiMG,
                               mflgCoarMG,
                               m_amrInterpolators[amrLevel],
                               m_amrFluxRegisters[amrLevel],
                               m_amrCoarseners[amrLevel],
                               m_validCells[amrLevel],
                               m_domainBcFactory,
                               m_ebBcFactory,
                               m_jumpBcFactory,
                               m_probLo,
                               dx,
                               refToFine,
                               refToCoar,
                               hasFine,
                               hasCoar,
                               hasMGObjects,
                               false, // Not an MG operator.
                               m_alpha,
                               m_beta,
                               m_amrAcoef[amrLevel],
                               m_amrBcoef[amrLevel],
                               m_amrBcoefIrreg[amrLevel],
                               m_ghostPhi,
                               m_ghostRhs,
                               m_jumpOrder,
                               m_jumpWeight,
                               m_numPreCondSmooth,
                               m_smoother,
                               m_relaxFactor,
                               m_chebyOrder,
                               m_chebyEigRatio,
                               m_rasInnerSweeps,
                               m_refluxFree);

  // Give the operator access by reference to the jump data.
  op->setJump(m_amrJump[amrLevel]);

  return op;
}

bool
MFHelmholtzOpFactory::isCoarser(const ProblemDomain& A, const ProblemDomain& B)
{
  CH_TIME("MFHelmholtzOpFactory::isCoarser(ProblemDomain, ProblemDomain)");

  return A.domainBox().numPts() < B.domainBox().numPts();
}

bool
MFHelmholtzOpFactory::isFiner(const ProblemDomain& A, const ProblemDomain& B)
{
  CH_TIME("MFHelmholtzOpFactory::isFiner(ProblemDomain, ProblemDomain)");

  return A.domainBox().numPts() > B.domainBox().numPts();
}

int
MFHelmholtzOpFactory::refToFiner(const ProblemDomain& a_domain) const
{
  CH_TIME("MFHelmholtzOpFactory::refToFiner(ProblemDomain)");

  int ref = -1;

  for (int ilev = 0; ilev < m_amrLevelGrids.size(); ilev++) {
    if (m_amrLevelGrids[ilev].getDomain() == a_domain) {
      ref = m_amrRefRatios[ilev];
    }
  }

  return ref;
}

int
MFHelmholtzOpFactory::findAmrLevel(const ProblemDomain& a_domain) const
{
  CH_TIME("MFHelmholtzOpFactory::findAmrLevel(ProblemDomain)");

  int amrLevel = -1;
  for (int lvl = 0; lvl < m_amrLevelGrids.size(); lvl++) {
    if (m_amrLevelGrids[lvl].getDomain() == a_domain) {
      amrLevel = lvl;
      break;
    }
  }

  if (amrLevel < 0) {
    MayDay::Error("MFHelmholtzOpFactory::findAmrLevel - no corresponding amr level found!");
  }

  return amrLevel;
}

#include <CD_NamespaceFooter.H>
