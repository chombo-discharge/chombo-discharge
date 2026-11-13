/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBHelmholtzOpFactory.cpp
 * @brief  Implementation of CD_EBHelmholtzOpFactory.H
 * @author Robert Marskar
 */

// Chombo includes
#include <ParmParse.H>
#include <BaseIVFactory.H>
#include <CH_Timer.H>
#include <EBCellFactory.H>
#include <EBFluxFactory.H>

// Our includes
#include <CD_EBHelmholtzOpFactory.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzOpFactory::EBHelmholtzOpFactory(const Location::Cell         a_dataLocation,
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
                                           const EBBCFactory&           a_ebbcFactory,
                                           const IntVect&               a_ghostPhi,
                                           const IntVect&               a_ghostRhs,
                                           const Smoother&              a_smoother,
                                           const Real&                  a_relaxFactor,
                                           const int&                   a_chebyOrder,
                                           const Real&                  a_chebyEigRatio,
                                           const int&                   a_rasInnerSweeps,
                                           const Vector<AmrLevelGrids>& a_mgLevelGrids,
                                           const Vector<AmrLevelGrids>& a_mgLevelGridsCoFi,
                                           const Vector<AmrCoarseners>& a_mgCoarseners,
                                           const AmrLevelGrids&         a_amrLevelGridsCoFi,
                                           const ProblemDomain&         a_bottomDomain,
                                           const bool                   a_refluxFree)
  : m_dataLocation(a_dataLocation),
    m_smoother(a_smoother),
    m_chebyOrder(a_chebyOrder),
    m_chebyEigRatio(a_chebyEigRatio),
    m_rasInnerSweeps(a_rasInnerSweeps),
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
    m_ebBcFactory(a_ebbcFactory),
    m_bottomDomain(a_bottomDomain),
    m_amrLevelGridsCoFi(a_amrLevelGridsCoFi)
{
  CH_TIME("EBHelmholtzOpFactory::EBHelmholtzOpFactory(...)");

  // Define constructor arguments.

  m_numAmrLevels = static_cast<int>(m_amrLevelGrids.size());

  // Asking multigrid to do the bottom solve at a refined AMR level is classified as bad input.
  if (ChomboDischarge::EBHelmholtzOpFactory::isFiner(m_bottomDomain, m_amrLevelGrids[0]->getDomain())) {
    MayDay::Error("EBHelmholtzOpFactory -- bottomsolver domain can't be larger than the base AMR domain!");
  }

  // Define the multigrid levels.
  this->defineMultigridLevels(a_mgLevelGrids, a_mgLevelGridsCoFi, a_mgCoarseners);
}

EBHelmholtzOpFactory::~EBHelmholtzOpFactory()
{
  CH_TIME("EBHelmholtzOpFactory::~EBHelmholtzOpFactory()");
}

void
EBHelmholtzOpFactory::defineMultigridLevels(const Vector<AmrLevelGrids>& a_mgLevelGrids,
                                            const Vector<AmrLevelGrids>& a_mgLevelGridsCoFi,
                                            const Vector<AmrCoarseners>& a_mgCoarseners)
{
  CH_TIME("EBHelmholtzOpFactory::defineMultigridLevels(Vector<AmrLevelGrids>, Vector<AmrLevelGrids>, "
          "Vector<AmrCoarseners>)");

  // TLDR: This routine defines what is needed for making the multigrid levels. This includes the intermediate
  // levels (if you run with refinement factor 4) as well as the deeper multigrid levels that are coarsenings of the
  // base AMR level. Recall that in Chombo-speak a multigrid level is a grid level completely covered by another grid
  // level. The realm builds the grids; this factory takes the entries that are not coarser than its bottom domain and
  // defines the coefficients on them.

  CH_assert(static_cast<int>(a_mgLevelGrids.size()) == m_numAmrLevels);
  CH_assert(static_cast<int>(a_mgLevelGridsCoFi.size()) == m_numAmrLevels);
  CH_assert(static_cast<int>(a_mgCoarseners.size()) == m_numAmrLevels);

  // Resize data holders.
  m_mgLevelGrids.resize(m_numAmrLevels);
  m_mgAcoef.resize(m_numAmrLevels);
  m_mgBcoef.resize(m_numAmrLevels);
  m_mgBcoefIrreg.resize(m_numAmrLevels);
  m_mgCoarseners.resize(m_numAmrLevels);
  m_mgLevelGridsCoFi.resize(m_numAmrLevels);
  m_hasMgLevels.resize(m_numAmrLevels);

  for (int amrLevel = 0; amrLevel < m_numAmrLevels; amrLevel++) {
    const AmrLevelGrids& realmLevels = a_mgLevelGrids[amrLevel];

    CH_assert(realmLevels.size() >= 1);
    CH_assert(realmLevels[0]->getDomain() == m_amrLevelGrids[amrLevel]->getDomain());

    // Use the realm's levels down to, and including, the last one that is not coarser than the bottom domain.
    int numMgLevels = 1;

    while (numMgLevels < static_cast<int>(realmLevels.size()) &&
           !(ChomboDischarge::EBHelmholtzOpFactory::isCoarser(realmLevels[numMgLevels]->getDomain(), m_bottomDomain))) {
      numMgLevels++;
    }

    m_hasMgLevels[amrLevel] = numMgLevels > 1;

    // The realm builds down to the coarsest bottom domain registered on it, so it stops before this factory's bottom
    // domain only if the grids could not be coarsened further, or if the bottom cell count was lowered after
    // registration.
    if (amrLevel == 0 &&
        ChomboDischarge::EBHelmholtzOpFactory::isCoarser(m_bottomDomain, realmLevels.back()->getDomain())) {
      pout() << "EBHelmholtzOpFactory::defineMultigridLevels - the multigrid levels stop at "
             << realmLevels.back()->getDomain() << ", which is finer than the bottom domain " << m_bottomDomain << endl;
    }

    m_mgLevelGrids[amrLevel].resize(0);
    m_mgAcoef[amrLevel].resize(0);
    m_mgBcoef[amrLevel].resize(0);
    m_mgBcoefIrreg[amrLevel].resize(0);
    m_mgCoarseners[amrLevel].resize(0);
    m_mgLevelGridsCoFi[amrLevel].resize(0);

    if (m_hasMgLevels[amrLevel]) {

      // Note that m_mgLevelGrids[amrLevel] corresponds to the multigrid levels below amrLevel. With this indexing the
      // starting index is the finest level, and elements later in the vector are coarser levels.
      m_mgLevelGrids[amrLevel].push_back(m_amrLevelGrids[amrLevel]);
      m_mgAcoef[amrLevel].push_back(m_amrAcoef[amrLevel]);
      m_mgBcoef[amrLevel].push_back(m_amrBcoef[amrLevel]);
      m_mgBcoefIrreg[amrLevel].push_back(m_amrBcoefIrreg[amrLevel]);

      for (int img = 1; img < numMgLevels; img++) {
        const RefCountedPtr<EBLevelGrid>& mgEblgCoar = realmLevels[img];

        const EBLevelGrid& eblgCoar = *mgEblgCoar;

        CH_assert(eblgCoar.getGhost() >= m_ghostPhi.max());

        const EBISLayout& ebislCoar = eblgCoar.getEBISL();

        const DisjointBoxLayout& gridsCoar = eblgCoar.getDBL();

        const ProblemDomain& domainCoar = eblgCoar.getDomain();

        // Make the irregular sets.
        constexpr int nghost = 1;

        LayoutData<IntVectSet> irregSets(gridsCoar);

        const DataIterator& dit = gridsCoar.dataIterator();

        const int nbox = dit.size();
#pragma omp parallel for schedule(runtime)
        for (int mybox = 0; mybox < nbox; mybox++) {
          const DataIndex& din = dit[mybox];

          Box bx = gridsCoar[din];
          bx.grow(nghost);
          bx &= domainCoar;

          const EBISBox& ebisbox = ebislCoar[din];

          irregSets[din] = ebisbox.getIrregIVS(bx);
        }

        // Factories for generating A and B-coefficients on the multigrid levels.
        EBCellFactory       cellFactory(ebislCoar);
        EBFluxFactory       fluxFactory(ebislCoar);
        BaseIVFactory<Real> irreFactory(ebislCoar, irregSets);

        // Define coarsened coefficients
        RefCountedPtr<LevelData<EBCellFAB>> coarAcoef(
          new LevelData<EBCellFAB>(gridsCoar, m_nComp, nghost * IntVect::Unit, cellFactory));
        RefCountedPtr<LevelData<EBFluxFAB>> coarBcoef(
          new LevelData<EBFluxFAB>(gridsCoar, m_nComp, nghost * IntVect::Unit, fluxFactory));
        RefCountedPtr<LevelData<BaseIVFAB<Real>>> coarBcoefIrreg(
          new LevelData<BaseIVFAB<Real>>(gridsCoar, m_nComp, nghost * IntVect::Unit, irreFactory));

        m_mgLevelGrids[amrLevel].push_back(mgEblgCoar);
        m_mgAcoef[amrLevel].push_back(coarAcoef);
        m_mgBcoef[amrLevel].push_back(coarBcoef);
        m_mgBcoefIrreg[amrLevel].push_back(coarBcoefIrreg);
        m_mgCoarseners[amrLevel].push_back(a_mgCoarseners[amrLevel][img - 1]);
        m_mgLevelGridsCoFi[amrLevel].push_back(a_mgLevelGridsCoFi[amrLevel][img - 1]);
      }
    }
  }

  this->coarsenCoefficientsMG();
}

void
EBHelmholtzOpFactory::coarsenCoefficientsMG()
{
  CH_TIME("EBHelmholtzOpFactory::coarsenCoefficientsMG");

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
        LevelData<EBCellFAB>&       coarAcoef      = *mgAco[mgLevel + 1];
        LevelData<EBFluxFAB>&       coarBcoef      = *mgBco[mgLevel + 1];
        LevelData<BaseIVFAB<Real>>& coarBcoefIrreg = *mgBcoIrreg[mgLevel + 1];

        const LevelData<EBCellFAB>&       fineAcoef      = *mgAco[mgLevel];
        const LevelData<EBFluxFAB>&       fineBcoef      = *mgBco[mgLevel];
        const LevelData<BaseIVFAB<Real>>& fineBcoefIrreg = *mgBcoIrreg[mgLevel];

        const EBCoarAve& averageOp = *m_mgCoarseners[amrLevel][mgLevel];

        averageOp.averageData(coarAcoef, fineAcoef, interv, average);
        averageOp.averageData(coarBcoef, fineBcoef, interv, average);
        averageOp.averageData(coarBcoefIrreg, fineBcoefIrreg, interv, average);

        coarAcoef.exchange();
        coarBcoef.exchange();
        coarBcoefIrreg.exchange();
      }
    }
  }
}

bool
EBHelmholtzOpFactory::isCoarser(const ProblemDomain& A, const ProblemDomain& B)
{
  CH_TIME("EBHelmholtzOpFactory::isCoarser(ProblemDomain, ProblemDomain)");

  return A.domainBox().numPts() < B.domainBox().numPts();
}

bool
EBHelmholtzOpFactory::isFiner(const ProblemDomain& A, const ProblemDomain& B)
{
  CH_TIME("EBHelmholtzOpFactory::isFiner(ProblemDomain, ProblemDomain)");

  return A.domainBox().numPts() > B.domainBox().numPts();
}

EBHelmholtzOp*
EBHelmholtzOpFactory::MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool /*a_homogeneousOnly*/)
{
  CH_TIME("EBHelmholtzOpFactory::MGnewOp(ProblemDomain, int, bool)");

  EBHelmholtzOp* mgOp = nullptr;

  // First, check if multigrid should be turned off completely.
  ParmParse pp;
  bool      turnOffMG = false;
  pp.query("turn_off_multigrid", turnOffMG);
  if (turnOffMG) {
    pout() << "turning off multigrid for EBHelmholtzOpFactory because 'turn_off_multigrid' = true" << endl;
    return mgOp;
  }

  // Find the AMR level corresponding to a coarsen of a_fineDomain by a factor 2^a_depth. This will issue a run-time
  // error if the domain is not found!
  const int amrLevel = this->findAmrLevel(a_fineDomain);

  constexpr int mgRefRat = 2;

  // Things that are needed for defining the operator. This might seem weird but
  // the multigrid operators only do relaxation and there's no fine-to-coar stuff. So hasFine = hasCoar = false.
  // But we might need to restrict and interpolate another, even coarser multigrid level so in that case we need
  // to define that.
  EBLevelGrid eblg;
  EBLevelGrid eblgMgCoFi;
  EBLevelGrid eblgMgCoar;

  bool hasMGObjects;

  RefCountedPtr<EBMultigridInterpolator> interpolator; // Only if defined on an AMR level
  RefCountedPtr<EBReflux>                fluxReg;      // Only if defined on an AMR level
  RefCountedPtr<EBCoarAve>               coarsener;    // Only if defined on an AMR level

  RefCountedPtr<LevelData<BaseFab<bool>>>   validCells; // Might be nullptr
  RefCountedPtr<LevelData<EBCellFAB>>       Acoef;      // Always defined.
  RefCountedPtr<LevelData<EBFluxFAB>>       Bcoef;      // Always defined.
  RefCountedPtr<LevelData<BaseIVFAB<Real>>> BcoefIrreg; // Always defined.

  bool foundMgLevel = false;

  if (a_depth == 0) { // Asking for the AMR level.
    eblg       = *m_amrLevelGrids[amrLevel];
    Acoef      = m_amrAcoef[amrLevel];
    Bcoef      = m_amrBcoef[amrLevel];
    BcoefIrreg = m_amrBcoefIrreg[amrLevel];

    validCells   = m_validCells[amrLevel];
    interpolator = m_amrInterpolators[amrLevel];
    fluxReg      = m_amrFluxRegisters[amrLevel];
    coarsener    = m_amrCoarseners[amrLevel];

    hasMGObjects = m_hasMgLevels[amrLevel];

    if (hasMGObjects) {
      eblgMgCoFi = *m_mgLevelGridsCoFi[amrLevel][0];
      eblgMgCoar = *m_mgLevelGrids[amrLevel][1];
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
      if (mgLevelGrids[img]->getDomain() == coarDomain) {
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
      eblg       = *mgLevelGrids[mgLevel];

      // This just means that mgLevel was not the last entry in mgLevelGrids so there's even coarser stuff below.
      hasMGObjects = (mgLevel < mgLevelGrids.size() - 1);
      if (hasMGObjects) {
        eblgMgCoFi = *m_mgLevelGridsCoFi[amrLevel][mgLevel];
        eblgMgCoar = *mgLevelGrids[mgLevel + 1];
      }
    }
  }

  // Make the operator
  if (foundMgLevel) {
    const Real dx = m_amrResolutions[amrLevel] * std::pow(mgRefRat, a_depth); //

    auto domBC = m_domainBcFactory->create();
    auto ebBC  = m_ebBcFactory->create();

    mgOp = new EBHelmholtzOp(m_dataLocation,
                             EBLevelGrid(), // Multigrid operator, so no fine.
                             eblg,
                             EBLevelGrid(), // Multigrid operator, so no cofi.
                             EBLevelGrid(), // Multigrid operator, so no coarse.
                             eblgMgCoFi,
                             eblgMgCoar,
                             validCells,   // Defined if an amr level
                             interpolator, // Defined if an amr level
                             fluxReg,      // Defined if an amr level
                             coarsener,    // Defined if an amr level
                             domBC,
                             ebBC,
                             m_probLo,
                             dx,    // Set from depth
                             1,     // Multigrid operator. Set to 1 in operator anyways.
                             1,     // Multigrid operator. Set to 1 in operator anyways.
                             false, // Multigrid operator, so false.
                             false, // Multigrid operator, so false.
                             hasMGObjects,
                             m_alpha,
                             m_beta,
                             Acoef,
                             Bcoef,
                             BcoefIrreg,
                             m_ghostPhi,
                             m_ghostRhs,
                             m_smoother,
                             m_relaxFactor,
                             m_chebyOrder,
                             m_chebyEigRatio,
                             m_rasInnerSweeps,
                             m_refluxFree);
  }

  return mgOp;
}

EBHelmholtzOp*
EBHelmholtzOpFactory::AMRnewOp(const ProblemDomain& a_domain)
{
  CH_TIME("EBHelmholtzOpFactory::AMRnewOp(ProblemDomain)");

  EBHelmholtzOp* op;

  const int amrLevel = this->findAmrLevel(a_domain);

  const bool hasFine = amrLevel < m_numAmrLevels - 1;
  const bool hasCoar = amrLevel > 0;

  EBLevelGrid eblgFine;
  EBLevelGrid eblg;
  EBLevelGrid eblgCoFi;
  EBLevelGrid eblgCoar;
  EBLevelGrid eblgCoFiMG;
  EBLevelGrid eblgCoarMG;

  Real dx;

  eblg = *m_amrLevelGrids[amrLevel];
  dx   = m_amrResolutions[amrLevel];

  int refToCoar = 1;
  int refToFine = 1;

  if (hasCoar) {
    eblgCoar  = *m_amrLevelGrids[amrLevel - 1];
    refToCoar = m_amrRefRatios[amrLevel - 1];
  }

  if (hasFine) {
    eblgFine  = *m_amrLevelGrids[amrLevel + 1];
    refToFine = m_amrRefRatios[amrLevel];
  }

  const bool hasMGObjects = m_hasMgLevels[amrLevel];
  if (hasMGObjects) {
    eblgCoFiMG = *m_mgLevelGridsCoFi[amrLevel][0];
    eblgCoarMG = *m_mgLevelGrids[amrLevel][1];
    CH_assert(eblgCoarMG.isDefined());
  }

  if (hasCoar) {
    eblgCoFi = *m_amrLevelGridsCoFi[amrLevel - 1];
    CH_assert(eblgCoFi.isDefined());
  }

  auto domainBC = m_domainBcFactory->create();
  auto ebBC     = m_ebBcFactory->create();

  op = new EBHelmholtzOp(m_dataLocation,
                         eblgFine,
                         eblg,
                         eblgCoFi,
                         eblgCoar,
                         eblgCoFiMG,
                         eblgCoarMG,
                         m_validCells[amrLevel],
                         m_amrInterpolators[amrLevel],
                         m_amrFluxRegisters[amrLevel],
                         m_amrCoarseners[amrLevel],
                         domainBC,
                         ebBC,
                         m_probLo,
                         dx,
                         refToFine,
                         refToCoar,
                         hasFine,
                         hasCoar,
                         hasMGObjects,
                         m_alpha,
                         m_beta,
                         m_amrAcoef[amrLevel],
                         m_amrBcoef[amrLevel],
                         m_amrBcoefIrreg[amrLevel],
                         m_ghostPhi,
                         m_ghostRhs,
                         m_smoother,
                         m_relaxFactor,
                         m_chebyOrder,
                         m_chebyEigRatio,
                         m_rasInnerSweeps,
                         m_refluxFree);

  return op;
}

int
EBHelmholtzOpFactory::refToFiner(const ProblemDomain& a_domain) const
{
  CH_TIME("EBHelmholtzOpFactory::refToFiner(ProblemDomain)");

  int  ref   = -1;
  bool found = false;

  for (int ilev = 0; ilev < m_amrLevelGrids.size(); ilev++) {
    if (m_amrLevelGrids[ilev]->getDomain() == a_domain) {
      ref   = m_amrRefRatios[ilev];
      found = true;
    }
  }

  // I will call this an error
  if (!found) {
    MayDay::Error("EBHelmholtzOpFactory::refToFiner - Domain not found in the AMR hierarchy");
  }

  return ref;
}

int
EBHelmholtzOpFactory::findAmrLevel(const ProblemDomain& a_domain) const
{
  CH_TIME("EBHelmholtzOpFactory::findAmrLevel(ProblemDomain)");

  int amrLevel = -1;
  for (int lvl = 0; lvl < m_amrLevelGrids.size(); lvl++) {
    if (m_amrLevelGrids[lvl]->getDomain() == a_domain) {
      amrLevel = lvl;
      break;
    }
  }

  if (amrLevel < 0) {
    MayDay::Error("EBHelmholtzOpFactory::findAmrLevel - no corresponding amr level found!");
  }

  return amrLevel;
}

#include <CD_NamespaceFooter.H>
