/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOpFactory.cpp
  @brief  Implementation of CD_MFHelmholtzOpFactory.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <BRMeshRefine.H>
#include <LoadBalance.H>
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

MFHelmholtzOpFactory::MFHelmholtzOpFactory(const MFIS&             a_mfis,
                                           const Location::Cell    a_dataLocation,
                                           const Real&             a_alpha,
                                           const Real&             a_beta,
                                           const RealVect&         a_probLo,
                                           const AmrLevelGrids&    a_amrLevelGrids,
                                           const AmrInterpolators& a_amrInterpolators,
                                           const AmrFluxRegisters& a_amrFluxRegisters,
                                           const AmrCoarseners&    a_amrCoarseners,
                                           const AmrRefRatios&     a_amrRefRatios,
                                           const AmrResolutions&   a_amrResolutions,
                                           const AmrCellData&      a_amrAcoef,
                                           const AmrFluxData&      a_amrBcoef,
                                           const AmrIrreData&      a_amrBcoefIrreg,
                                           const DomainBCFactory&  a_domainBcFactory,
                                           const EBBCFactory&      a_ebBcFactory,
                                           const JumpBCFactory&    a_jumpBcFactory,
                                           const IntVect&          a_ghostPhi,
                                           const IntVect&          a_ghostRhs,
                                           const Smoother&         a_smoother,
                                           const ProblemDomain&    a_bottomDomain,
                                           const int&              a_jumpOrder,
                                           const int&              a_jumpWeight,
                                           const int&              a_blockingFactor,
                                           const AmrLevelGrids&    a_deeperLevelGrids)
{
  CH_TIME("MFHelmholtzOpFactory::MFHelmholtzOpFactory()");

  m_mfis         = a_mfis;
  m_dataLocation = a_dataLocation;

  m_alpha  = a_alpha;
  m_beta   = a_beta;
  m_probLo = a_probLo;

  m_amrLevelGrids    = a_amrLevelGrids;
  m_amrInterpolators = a_amrInterpolators;
  m_amrFluxRegisters = a_amrFluxRegisters;
  m_amrCoarseners    = a_amrCoarseners;
  m_amrRefRatios     = a_amrRefRatios;
  m_amrResolutions   = a_amrResolutions;

  m_amrAcoef      = a_amrAcoef;
  m_amrBcoef      = a_amrBcoef;
  m_amrBcoefIrreg = a_amrBcoefIrreg;

  m_domainBcFactory = a_domainBcFactory;
  m_ebBcFactory     = a_ebBcFactory;
  m_jumpBcFactory   = a_jumpBcFactory;

  m_ghostPhi = a_ghostPhi;
  m_ghostRhs = a_ghostRhs;

  m_smoother     = a_smoother;
  m_bottomDomain = a_bottomDomain;

  m_jumpOrder  = a_jumpOrder;
  m_jumpWeight = a_jumpWeight;

  m_mgBlockingFactor = a_blockingFactor;

  m_deeperLevelGrids = a_deeperLevelGrids;

  m_numAmrLevels = m_amrLevelGrids.size();

  // Asking multigrid to do the bottom solve at a refined AMR level is classified as bad input.
  if (this->isFiner(m_bottomDomain, m_amrLevelGrids[0].getDomain())) {
    MayDay::Abort("MFHelmholtzOpFactory -- bottomsolver domain can't be larger than the base AMR domain!");
  }

  // Define the jump data and the multigrid levels.
  this->defineJump();
  this->defineMultigridLevels();
  this->setJump(0.0, 0.0);
}

MFHelmholtzOpFactory::~MFHelmholtzOpFactory() { CH_TIME("MFHelmholtzOpFactory::~MFHelmholtzOpFactory()"); }

void
MFHelmholtzOpFactory::setJump(const EBAMRIVData& a_sigma, const Real& a_scale)
{
  CH_TIME("MFHelmholtzOpFactory::setJump(EBAMRIVData, Real)");

  // TLDR: This routine sets the jump coefficient on each level. Data on deeper levels
  //       are coarsened from finer ones.

  const Interval interv(m_comp, m_comp);

  for (int lvl = 0; lvl < m_numAmrLevels; lvl++) {
    const DisjointBoxLayout& dbl = a_sigma[lvl]->disjointBoxLayout();

    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box      box = dbl[dit()];
      const Interval interv(m_comp, m_comp);
      (*m_amrJump[lvl])[dit()].copy(box, interv, box, (*a_sigma[lvl])[dit()], interv);
    }

    DataOps::scale(*m_amrJump[lvl], a_scale);

    m_amrJump[lvl]->exchange();
  }

  // Average down on AMR levels.
  for (int lvl = m_numAmrLevels - 1; lvl > 0; lvl--) {
    const RefCountedPtr<EBCoarAve>& aveOp = m_amrCoarseners[lvl].getAveOp(m_mainPhase);

    const Average average = Average::Arithmetic;

    aveOp->averageData(*m_amrJump[lvl - 1], *m_amrJump[lvl], Interval(m_comp, m_comp), average);
  }

  // Average down on MG levels. I'm not really sure if we need to do this. Also,
  // note the weird reversed order for the levels.
#if 0
  // R.M. Nov. 2021: This broke when I compiled with DEBUG=TRUE. Not sure what causes it but it shouldn't matter because
  // the data is not used anywhere anyways (I think).
  for(int amrLevel = 0; amrLevel < m_numAmrLevels; amrLevel++){

    if(m_hasMgLevels[amrLevel]){
      EBAMRIVData& jumpMG = m_mgJump[amrLevel];

      const int finestMGLevel   = 0;
      const int coarsestMGLevel = jumpMG.size() - 1;
      for (int img = finestMGLevel+1; img <= coarsestMGLevel; img++){
	m_mgAveOp[amrLevel][img]->average(*jumpMG[img], *jumpMG[img-1], interv);
      }
    }
  }
#endif
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

  // TLDR: This defines m_amrJump on the first phase (gas phase). This is irregular data intended to be interfaced into the
  //       boundary condition class. When we match the BC we get the data from here (the gas phase). Note that we define
  //       m_amrJump on all irregular cells, but the operators will do matching on a subset of them.

  m_amrJump.resize(m_numAmrLevels);

  for (int lvl = 0; lvl < m_numAmrLevels; lvl++) {
    const DisjointBoxLayout& dbl   = m_amrLevelGrids[lvl].getGrids();
    const EBLevelGrid&       eblg  = m_amrLevelGrids[lvl].getEBLevelGrid(m_mainPhase);
    const EBISLayout&        ebisl = eblg.getEBISL();

    LayoutData<IntVectSet> irregCells(dbl);
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box      box     = dbl[dit()];
      const EBISBox& ebisbox = ebisl[dit()];

      irregCells[dit()] = ebisbox.getIrregIVS(box);
    }

    BaseIVFactory<Real> fact(ebisl, irregCells);
    m_amrJump[lvl] =
      RefCountedPtr<LevelData<BaseIVFAB<Real>>>(new LevelData<BaseIVFAB<Real>>(dbl, m_nComp, IntVect::Zero, fact));
  }
}

void
MFHelmholtzOpFactory::defineMultigridLevels()
{
  CH_TIME("MFHelmholtzOpFactory::defineMultigridLevels()");

  // TLDR: This routine defines what is needed for making the multigrid levels. This includes the intermediate
  // levels (if you run with refinement factor 4) as well as the deeper multigrid levels that are coarsenings of the base AMR level. Recall that
  // in Chombo-speak a multigrid level is a grid level completely covered by another grid level.

  m_mgLevelGrids.resize(m_numAmrLevels);
  m_mgAcoef.resize(m_numAmrLevels);
  m_mgBcoef.resize(m_numAmrLevels);
  m_mgBcoefIrreg.resize(m_numAmrLevels);
  m_mgJump.resize(m_numAmrLevels);
  m_hasMgLevels.resize(m_numAmrLevels);
  m_mgAveOp.resize(m_numAmrLevels);

  for (int amrLevel = 0; amrLevel < m_numAmrLevels; amrLevel++) {
    m_hasMgLevels[amrLevel] = false;

    // We can have a multigrid level either if the refinement factor to the coarse level is larger than two, or we are at the bottom
    // of the AMR hierarchy.
    if (amrLevel == 0 && this->isCoarser(m_bottomDomain, m_amrLevelGrids[amrLevel].getDomain())) {
      m_hasMgLevels[amrLevel] = true;
    }

    if (amrLevel > 0) {
      if (m_amrRefRatios[amrLevel - 1] > 2) {
        m_hasMgLevels[amrLevel] = true;
      }
    }

    // Create MG levels.
    if (m_hasMgLevels[amrLevel]) {

      constexpr int mgRefRatio = 2;

      m_mgLevelGrids[amrLevel].resize(0);
      m_mgAcoef[amrLevel].resize(0);
      m_mgBcoef[amrLevel].resize(0);
      m_mgBcoefIrreg[amrLevel].resize(0);
      m_mgJump[amrLevel].resize(0);
      m_mgAveOp[amrLevel].resize(0);

      m_mgLevelGrids[amrLevel].push_back(m_amrLevelGrids[amrLevel]);
      m_mgAcoef[amrLevel].push_back(m_amrAcoef[amrLevel]);
      m_mgBcoef[amrLevel].push_back(m_amrBcoef[amrLevel]);
      m_mgBcoefIrreg[amrLevel].push_back(m_amrBcoefIrreg[amrLevel]);
      m_mgJump[amrLevel].push_back(m_amrJump[amrLevel]);
      m_mgAveOp[amrLevel].push_back(RefCountedPtr<EBCoarAve>(nullptr));

      bool hasCoarser = true;

      while (hasCoarser) {
        const int          curMgLevels = m_mgLevelGrids[amrLevel].size();
        const MFLevelGrid& mgMflgFine  = m_mgLevelGrids[amrLevel].back();

        // This is the one we will define
        MFLevelGrid mgMflgCoar;

        // This is an overriding option where we use the pre-defined coarsenings in m_deeperMultigridLevels. This is only valid for coarsenings of
        // the base AMR level, hence amrLevel == 0. Once those levels are exhausted we begin with direct coarsening.
        if (amrLevel == 0 && curMgLevels < m_deeperLevelGrids.size()) {
          hasCoarser = true; // Note that m_deeperLevelGrids[0] should be a factor 2 coarsening of the
          mgMflgCoar = m_deeperLevelGrids[curMgLevels - 1]; // coarsest AMR level. So curMgLevels-1 is correct.
        }
        else {
          // Let the operator factory do the coarsening this time.
          hasCoarser = this->getCoarserLayout(mgMflgCoar, mgMflgFine, mgRefRatio, m_mgBlockingFactor);
        }

        // Do not coarsen further if we end up with a domain smaller than m_bottomDomain. In this case
        // we will terminate the coarsening and let AMRMultiGrid do the bottom solve.
        if (hasCoarser) {
          if (this->isCoarser(mgMflgCoar.getDomain(), m_bottomDomain)) {
            hasCoarser = false;
          }
          else {
            // Not so sure about this one, will we ever be asked to make an coarsened MG level which is also an AMR level? If not, this code
            // will reduce the coarsening efforts.
            for (int iamr = 0; iamr < m_numAmrLevels; iamr++) {
              if (mgMflgCoar.getDomain() == m_amrLevelGrids[iamr].getDomain()) {
                hasCoarser = false;
              }
            }
          }
        }

        // Ok, we can coarsen the domain define by mgMflgCoar. Set up that domain as a multigrid level.
        if (hasCoarser) {

          const DisjointBoxLayout& dblCoar = mgMflgCoar.getGrids();

          // The Chombo multifluid things need Vector<EBISLayout> for the factories.
          Vector<int>        ebislComps;
          Vector<EBISLayout> ebislCoar;
          for (int i = 0; i < mgMflgCoar.numPhases(); i++) {
            ebislComps.push_back(m_nComp);
            ebislCoar.push_back(mgMflgCoar.getEBLevelGrid(i).getEBISL());
          }

          // Factories for making coarse stuff. Need one ghost cell because we interpolate b*grad(phi) to face centroids.
          const int           nghost = 1;
          MFCellFactory       cellFact(ebislCoar, ebislComps);
          MFFluxFactory       fluxFact(ebislCoar, ebislComps);
          MFBaseIVFABFactory  ivFact(ebislCoar, ebislComps);
          BaseIVFactory<Real> irregFact(mgMflgCoar.getEBLevelGrid(m_mainPhase).getEBISL());

          // Multifluid a bit special -- number of components come in through the factory.
          const int dummy = 1;

          auto coarAcoef = RefCountedPtr<LevelData<MFCellFAB>>(
            new LevelData<MFCellFAB>(dblCoar, dummy, nghost * IntVect::Unit, cellFact));
          auto coarBcoef = RefCountedPtr<LevelData<MFFluxFAB>>(
            new LevelData<MFFluxFAB>(dblCoar, dummy, nghost * IntVect::Unit, fluxFact));
          auto coarBcoefIrreg = RefCountedPtr<LevelData<MFBaseIVFAB>>(
            new LevelData<MFBaseIVFAB>(dblCoar, dummy, nghost * IntVect::Unit, ivFact));
          auto coarJump = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
            new LevelData<BaseIVFAB<Real>>(dblCoar, m_nComp, nghost * IntVect::Zero, irregFact));

          const LevelData<MFCellFAB>&   fineAcoef      = *m_mgAcoef[amrLevel].back();
          const LevelData<MFFluxFAB>&   fineBcoef      = *m_mgBcoef[amrLevel].back();
          const LevelData<MFBaseIVFAB>& fineBcoefIrreg = *m_mgBcoefIrreg[amrLevel].back();

          // Coarsening coefficients.
          this->coarsenCoefficients(*coarAcoef,
                                    *coarBcoef,
                                    *coarBcoefIrreg,
                                    fineAcoef,
                                    fineBcoef,
                                    fineBcoefIrreg,
                                    mgMflgCoar,
                                    mgMflgFine,
                                    mgRefRatio);
          DataOps::setValue(*coarJump, 0.0);

          // This is a special object for coarsening jump data between MG levels.
          const EBLevelGrid&       eblgFine = mgMflgFine.getEBLevelGrid(m_mainPhase);
          const EBLevelGrid&       eblgCoar = mgMflgCoar.getEBLevelGrid(m_mainPhase);
          RefCountedPtr<EBCoarAve> aveOp(new EBCoarAve(eblgFine.getDBL(),
                                                       eblgCoar.getDBL(),
                                                       eblgFine.getEBISL(),
                                                       eblgCoar.getEBISL(),
                                                       eblgCoar.getDomain(),
                                                       mgRefRatio,
                                                       eblgCoar.getEBIS()));

          // Append. Phew.
          m_mgLevelGrids[amrLevel].push_back(mgMflgCoar);
          m_mgAcoef[amrLevel].push_back(coarAcoef);
          m_mgBcoef[amrLevel].push_back(coarBcoef);
          m_mgBcoefIrreg[amrLevel].push_back(coarBcoefIrreg);
          m_mgJump[amrLevel].push_back(coarJump);
          m_mgAveOp[amrLevel].push_back(aveOp);
        }
      }
    }

    if (m_mgLevelGrids[amrLevel].size() <= 1) {
      m_hasMgLevels[amrLevel] = false;
    }
  }
}

bool
MFHelmholtzOpFactory::getCoarserLayout(MFLevelGrid&       a_coarMflg,
                                       const MFLevelGrid& a_fineMflg,
                                       const int          a_refRat,
                                       const int          a_blockingFactor) const
{
  CH_TIME("MFHelmholtzOpFactory::getCoarserLayout(MFLevelGrid, MFLevelGrid, int, int)");

  // TLDR: This creates a coarsening of a_fineGrid with refinement factor 2. The strategy is to first coarsen directly, leaving the box-to-rank array intact
  //       but where the boxes are coarsened by a factor of two. If that does not work we will try domain decomposition with the blocking factor.
  //
  //       This returns true if the fine grid fully covers the domain. The nature of this makes it
  //       always true for the "deeper" multigridlevels,  but not so for the intermediate levels.

  bool hasCoarser = false;

  // This returns true if the fine grid fully covers the domain. The nature of this makes it
  // always true for the "deeper" multigridlevels,  but not so for the intermediate levels.
  auto isFullyCovered = [&](const MFLevelGrid& a_mflg) -> bool {
    unsigned long long       numPtsLeft = a_mflg.getDomain().domainBox().numPts(); // Number of grid points in
    const DisjointBoxLayout& dbl        = a_mflg.getGrids();

    for (LayoutIterator lit = dbl.layoutIterator(); lit.ok(); ++lit) {
      numPtsLeft -= dbl[lit()].numPts();
    }

    return (numPtsLeft == 0);
  };

  const ProblemDomain      fineDomain = a_fineMflg.getDomain();
  const ProblemDomain      coarDomain = coarsen(fineDomain, a_refRat);
  const DisjointBoxLayout& fineDbl    = a_fineMflg.getGrids();
  DisjointBoxLayout        coarDbl;

  // Check if we can get a coarsenable domain. Don't want to coarsen to 1x1 so hence the factor of 2 in the test here.
  ProblemDomain test = fineDomain;
  if (refine(coarsen(test, 2 * a_refRat), 2 * a_refRat) == fineDomain) {

    // Use coarsening if we can
    if (fineDbl.coarsenable(2 * a_refRat)) {
      coarsen(coarDbl, fineDbl, a_refRat);

      hasCoarser = true;
    }
    else { // if(false){// { // Check if we can use box aggregation
      if (isFullyCovered(a_fineMflg)) {
        Vector<Box> boxes;
        Vector<int> procs;

        // We could have use our load balancing here, but I don't see why.
        domainSplit(coarDomain, boxes, a_blockingFactor);
        mortonOrdering(boxes);
        LoadBalance(procs, boxes);

        coarDbl.define(boxes, procs, coarDomain);

        hasCoarser = true;
      }
    }

    // Ok, found a coarsened layout.
    if (hasCoarser) {
      a_coarMflg = MFLevelGrid(coarDbl, coarDomain, m_ghostPhi.max(), m_mfis);
    }
  }
  else {
    hasCoarser = false;
  }

  return hasCoarser;
}

MFHelmholtzOp*
MFHelmholtzOpFactory::MGnewOp(const ProblemDomain& a_fineDomain, int a_depth, bool a_homogeneousOnly)
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

  const int amrLevel =
    this->findAmrLevel(a_fineDomain); // Run-time abort if a_fineDomain is not found in anhy amr level.

  const int mgRefRat = 2;

  // Things that are needed for defining the operator. This might seem weird but
  // the multigrid operators only do relaxation and there's no fine-to-coar stuff. So hasFine = hasCoar = false.
  // But we might need to restrict and interpolate another, even coarser multigrid level so in that case we need
  // to define that.
  MFLevelGrid mflg;
  MFLevelGrid mflgMgCoar;

  MFMultigridInterpolator interpolator; // Only if defined on an AMR level
  MFFluxReg               fluxReg;      // Only if defined on an AMR level
  MFCoarAve               coarsener;    // Only if defined on an AMR level

  bool hasMGObjects;

  RefCountedPtr<LevelData<MFCellFAB>>       Acoef;
  RefCountedPtr<LevelData<MFFluxFAB>>       Bcoef;
  RefCountedPtr<LevelData<MFBaseIVFAB>>     BcoefIrreg;
  RefCountedPtr<LevelData<BaseIVFAB<Real>>> jump;

  bool foundMgLevel = false;

  if (a_depth == 0) { // Asking for the AMR level.
    mflg       = m_amrLevelGrids[amrLevel];
    Acoef      = m_amrAcoef[amrLevel];
    Bcoef      = m_amrBcoef[amrLevel];
    BcoefIrreg = m_amrBcoefIrreg[amrLevel];
    jump       = m_amrJump[amrLevel];

    interpolator = m_amrInterpolators[amrLevel];
    fluxReg      = m_amrFluxRegisters[amrLevel];
    coarsener    = m_amrCoarseners[amrLevel];

    hasMGObjects = m_hasMgLevels[amrLevel];

    if (hasMGObjects) {
      mflgMgCoar = m_mgLevelGrids[amrLevel][1];
    }

    foundMgLevel = true;
  }
  else { // Asking for a coarsening. No interp or flux reg object here.
    // TLDR: Go through the coarsened levels for the specified amr level and see if we find a coarsening at the
    //       specified depth.
    const ProblemDomain coarDomain = coarsen(a_fineDomain, std::pow(mgRefRat, a_depth));

    // These are the things that live below the AMR level corresponding to a_fineDomain.
    const AmrLevelGrids& mgLevelGrids = m_mgLevelGrids[amrLevel];
    const AmrCellData&   mgAcoef      = m_mgAcoef[amrLevel];
    const AmrFluxData&   mgBcoef      = m_mgBcoef[amrLevel];
    const AmrIrreData&   mgBcoefIrreg = m_mgBcoefIrreg[amrLevel];

    // See if we have a corresponding multigrid level.
    int mgLevel;
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

      hasMGObjects =
        (mgLevel <
         mgLevelGrids.size() -
           1); // This just means that mgLevel was not the last entry in mgLevelGrids so there's even coarser stuff below.
      if (hasMGObjects) {
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
                             mflgMgCoar,
                             interpolator,
                             fluxReg,
                             coarsener,
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
                             m_smoother);

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
    mflgCoarMG = m_mgLevelGrids[amrLevel][1];
  }

  if (hasCoar) { // Make a coarser layout
    const bool gotCoarse = this->getCoarserLayout(mflgCoFi, mflg, refToCoar, m_mgBlockingFactor);
    mflgCoFi.setMaxRefinementRatio(refToCoar);
    CH_assert(gotCoarse);
  }

  MFHelmholtzOp* op = new MFHelmholtzOp(m_dataLocation,
                                        mflgFine,
                                        mflg,
                                        mflgCoFi,
                                        mflgCoar,
                                        mflgCoarMG,
                                        m_amrInterpolators[amrLevel],
                                        m_amrFluxRegisters[amrLevel],
                                        m_amrCoarseners[amrLevel],
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
                                        m_smoother);

  // Give the operator access by reference to the jump data.
  op->setJump(m_amrJump[amrLevel]);

  return op;
}

bool
MFHelmholtzOpFactory::isCoarser(const ProblemDomain& A, const ProblemDomain& B) const
{
  CH_TIME("MFHelmholtzOpFactory::isCoarser(ProblemDomain, ProblemDomain)");

  return A.domainBox().numPts() < B.domainBox().numPts();
}

bool
MFHelmholtzOpFactory::isFiner(const ProblemDomain& A, const ProblemDomain& B) const
{
  CH_TIME("MFHelmholtzOpFactory::isFiner(ProblemDomain, ProblemDomain)");

  return A.domainBox().numPts() > B.domainBox().numPts();
}

int
MFHelmholtzOpFactory::refToFiner(const ProblemDomain& a_domain) const
{
  CH_TIME("MFHelmholtzOpFactory::refToFiner(ProblemDomain)");

  int  ref   = -1;
  bool found = false;

  for (int ilev = 0; ilev < m_amrLevelGrids.size(); ilev++) {
    if (m_amrLevelGrids[ilev].getDomain() == a_domain) {
      found = true;
      ref   = m_amrRefRatios[ilev];
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

void
MFHelmholtzOpFactory::coarsenCoefficients(LevelData<MFCellFAB>&         a_coarAcoef,
                                          LevelData<MFFluxFAB>&         a_coarBcoef,
                                          LevelData<MFBaseIVFAB>&       a_coarBcoefIrreg,
                                          const LevelData<MFCellFAB>&   a_fineAcoef,
                                          const LevelData<MFFluxFAB>&   a_fineBcoef,
                                          const LevelData<MFBaseIVFAB>& a_fineBcoefIrreg,
                                          const MFLevelGrid&            a_mflgCoar,
                                          const MFLevelGrid&            a_mflgFine,
                                          const int                     a_refRat)
{
  CH_TIME("MFHelmholtzOpFactory::coarsenCoefficients(...)");

  const Interval interv(m_comp, m_comp);

  if (a_refRat == 1) {
    a_fineAcoef.copyTo(a_coarAcoef);
    a_fineBcoef.copyTo(a_coarBcoef);
    a_fineBcoefIrreg.copyTo(a_coarBcoefIrreg);
  }
  else {

    // Average down on each phase.
    for (int i = 0; i < a_mflgCoar.numPhases(); i++) {
      const EBLevelGrid& eblgCoar = a_mflgCoar.getEBLevelGrid(i);
      const EBLevelGrid& eblgFine = a_mflgFine.getEBLevelGrid(i);

      EBCoarAve aveOp(eblgFine.getDBL(),
                      eblgCoar.getDBL(),
                      eblgFine.getEBISL(),
                      eblgCoar.getEBISL(),
                      eblgCoar.getDomain(),
                      a_refRat,
                      eblgCoar.getEBIS());

      LevelData<EBCellFAB>       coarAco;
      LevelData<EBFluxFAB>       coarBco;
      LevelData<BaseIVFAB<Real>> coarBcoIrreg;

      LevelData<EBCellFAB>       fineAco;
      LevelData<EBFluxFAB>       fineBco;
      LevelData<BaseIVFAB<Real>> fineBcoIrreg;

      MultifluidAlias::aliasMF(coarAco, i, a_coarAcoef);
      MultifluidAlias::aliasMF(coarBco, i, a_coarBcoef);
      MultifluidAlias::aliasMF(coarBcoIrreg, i, a_coarBcoefIrreg);

      MultifluidAlias::aliasMF(fineAco, i, a_fineAcoef);
      MultifluidAlias::aliasMF(fineBco, i, a_fineBcoef);
      MultifluidAlias::aliasMF(fineBcoIrreg, i, a_fineBcoefIrreg);

      const Average average = Average::Arithmetic;

      aveOp.averageData(coarAco, fineAco, interv, average);
      aveOp.averageData(coarBco, fineBco, interv, average);
      aveOp.averageData(coarBcoIrreg, fineBcoIrreg, interv, average);

      coarAco.exchange();
      coarBco.exchange();
      coarBcoIrreg.exchange();
    }
  }
}

#include <CD_NamespaceFooter.H>
