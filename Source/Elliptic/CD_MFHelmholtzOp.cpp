/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzOp.cpp
  @brief  Implementation of CD_MFHelmholtzOp.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <CH_Timer.H>
#include <EBLevelDataOps.H>

// Our includes
#include <CD_Timer.H>
#include <CD_MFHelmholtzOp.H>
#include <CD_MultifluidAlias.H>
#include <CD_NamespaceHeader.H>

constexpr int MFHelmholtzOp::m_comp;
constexpr int MFHelmholtzOp::m_nComp;

MFHelmholtzOp::MFHelmholtzOp(const Location::Cell                             a_dataLocation,
                             const MFLevelGrid&                               a_mflgFine,
                             const MFLevelGrid&                               a_mflg,
                             const MFLevelGrid&                               a_mflgCoFi,
                             const MFLevelGrid&                               a_mflgCoar,
                             const MFLevelGrid&                               a_mflgCoarMG,
                             const MFMultigridInterpolator&                   a_interpolator,
                             const MFReflux&                                  a_fluxReg,
                             const MFCoarAve&                                 a_coarAve,
                             const RefCountedPtr<MFHelmholtzDomainBCFactory>& a_domainBcFactory,
                             const RefCountedPtr<MFHelmholtzEBBCFactory>&     a_ebBcFactory,
                             const RefCountedPtr<MFHelmholtzJumpBCFactory>&   a_jumpBcFactory,
                             const RealVect&                                  a_probLo,
                             const Real&                                      a_dx,
                             const int&                                       a_refToFine,
                             const int&                                       a_refToCoar,
                             const bool&                                      a_hasFine,
                             const bool&                                      a_hasCoar,
                             const bool&                                      a_hasMGObjects,
                             const bool&                                      a_isMGOperator,
                             const Real&                                      a_alpha,
                             const Real&                                      a_beta,
                             const RefCountedPtr<LevelData<MFCellFAB>>&       a_Acoef,
                             const RefCountedPtr<LevelData<MFFluxFAB>>&       a_Bcoef,
                             const RefCountedPtr<LevelData<MFBaseIVFAB>>&     a_BcoefIrreg,
                             const IntVect&                                   a_ghostPhi,
                             const IntVect&                                   a_ghostRhs,
                             const int&                                       a_jumpOrder,
                             const int&                                       a_jumpWeight,
                             const Smoother&                                  a_relaxType)
{
  CH_TIME("MFHelmholtzOp::MFHelmholtzOp(...)");

  CH_assert(!a_Acoef.isNull());
  CH_assert(!a_Bcoef.isNull());
  CH_assert(!a_BcoefIrreg.isNull());

  m_dataLocation = a_dataLocation;
  m_mflg         = a_mflg;
  m_numPhases    = m_mflg.numPhases();
  m_multifluid   = m_numPhases > 1;
  m_hasMGObjects = a_hasMGObjects;
  m_refToCoar    = a_refToCoar;
  m_smoother     = a_relaxType;
  m_hasCoar      = a_hasCoar;
  m_hasFine      = a_hasFine;
  m_Acoef        = a_Acoef;
  m_Bcoef        = a_Bcoef;
  m_BcoefIrreg   = a_BcoefIrreg;
  m_ghostPhi     = a_ghostPhi;
  m_ghostRhs     = a_ghostRhs;

  if (a_hasCoar) {
    m_mflgCoFi = a_mflgCoFi;
    m_mflgCoar = a_mflgCoar;
  }

  if (m_hasFine) {
    m_coarAve = a_coarAve;
  }

  m_interpolator = a_interpolator;

  if (m_hasMGObjects) {
    m_mflgCoarMG = a_mflgCoarMG;
  }

  EBArith::getMultiColors(m_colors);

  // Instantiate jump bc object.
  const int ghostCF = a_hasCoar ? a_interpolator.getGhostCF() : 1;
  m_jumpBC =
    a_jumpBcFactory
      ->create(m_dataLocation, m_mflg, a_BcoefIrreg, a_dx, a_jumpOrder, a_jumpWeight, a_jumpOrder, ghostCF, a_ghostPhi);

  // Make the operators on eachphase.
  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    EBLevelGrid eblg = a_mflg.getEBLevelGrid(iphase);
    EBLevelGrid dummy;
    EBLevelGrid eblgFine   = a_hasFine ? a_mflgFine.getEBLevelGrid(iphase) : dummy;
    EBLevelGrid eblgCoFi   = a_hasCoar ? a_mflgCoFi.getEBLevelGrid(iphase) : dummy;
    EBLevelGrid eblgCoar   = a_hasCoar ? a_mflgCoar.getEBLevelGrid(iphase) : dummy;
    EBLevelGrid eblgCoarMG = a_hasMGObjects ? a_mflgCoarMG.getEBLevelGrid(iphase) : dummy;

    RefCountedPtr<EBMultigridInterpolator> interpolator = RefCountedPtr<EBMultigridInterpolator>(nullptr);
    RefCountedPtr<EBReflux>                fluxRegister = RefCountedPtr<EBReflux>(nullptr);
    RefCountedPtr<EBCoarAve>               coarsener    = RefCountedPtr<EBCoarAve>(nullptr);

    if (!a_isMGOperator) {
      if (a_hasFine) {
        fluxRegister = a_fluxReg.getFluxRegPointer(iphase);
      }

      coarsener = a_coarAve.getAveOp(iphase);

      if (a_hasCoar) {
        interpolator = a_interpolator.getInterpolator(iphase);
      }
    }

    auto domainBC = a_domainBcFactory->create(iphase);
    auto ebBC     = a_ebBcFactory->create(iphase, m_jumpBC);

    // Alias the multifluid-coefficients onto a single phase.
    RefCountedPtr<LevelData<EBCellFAB>>       Acoef = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());
    RefCountedPtr<LevelData<EBFluxFAB>>       Bcoef = RefCountedPtr<LevelData<EBFluxFAB>>(new LevelData<EBFluxFAB>());
    RefCountedPtr<LevelData<BaseIVFAB<Real>>> BcoefIrreg = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
      new LevelData<BaseIVFAB<Real>>());

    MultifluidAlias::aliasMF(*Acoef, iphase, *a_Acoef);
    MultifluidAlias::aliasMF(*Bcoef, iphase, *a_Bcoef);
    MultifluidAlias::aliasMF(*BcoefIrreg, iphase, *a_BcoefIrreg);

    EBHelmholtzOp::Smoother ebHelmRelax;

    switch (a_relaxType) {
    case MFHelmholtzOp::Smoother::PointJacobi: {
      ebHelmRelax = EBHelmholtzOp::Smoother::PointJacobi;

      break;
    }
    case MFHelmholtzOp::Smoother::GauSaiRedBlack: {
      ebHelmRelax = EBHelmholtzOp::Smoother::GauSaiRedBlack;

      break;
    }
    case MFHelmholtzOp::Smoother::GauSaiMultiColor: {
      ebHelmRelax = EBHelmholtzOp::Smoother::GauSaiMultiColor;

      break;
    }
    default: {
      MayDay::Error("MFHelmholtzOp::MFHelmholtzOp - unsupported relaxation method requested");

      break;
    }
    }

    RefCountedPtr<EBHelmholtzOp> oper = RefCountedPtr<EBHelmholtzOp>(new EBHelmholtzOp(m_dataLocation,
                                                                                       eblgFine,
                                                                                       eblg,
                                                                                       eblgCoFi,
                                                                                       eblgCoar,
                                                                                       eblgCoarMG,
                                                                                       interpolator,
                                                                                       fluxRegister,
                                                                                       coarsener,
                                                                                       domainBC,
                                                                                       ebBC,
                                                                                       a_probLo,
                                                                                       a_dx,
                                                                                       a_refToFine,
                                                                                       a_refToCoar,
                                                                                       a_hasFine,
                                                                                       a_hasCoar,
                                                                                       a_hasMGObjects,
                                                                                       a_alpha,
                                                                                       a_beta,
                                                                                       Acoef,
                                                                                       Bcoef,
                                                                                       BcoefIrreg,
                                                                                       a_ghostPhi,
                                                                                       a_ghostRhs,
                                                                                       ebHelmRelax));

    m_helmOps.insert({iphase, oper});
  }
}

MFHelmholtzOp::~MFHelmholtzOp()
{
  CH_TIME("MFHelmholtzOp()::~MFHelmholtzOp()");

  m_helmOps.clear();
}

void
MFHelmholtzOp::setAcoAndBco(const RefCountedPtr<LevelData<MFCellFAB>>&   a_Acoef,
                            const RefCountedPtr<LevelData<MFFluxFAB>>&   a_Bcoef,
                            const RefCountedPtr<LevelData<MFBaseIVFAB>>& a_BcoefIrreg)
{
  CH_TIME("MFHelmholtzOp::setAcoAndBco");

  // Make the operators on eachphase.
  for (int iphase = 0; iphase < m_numPhases; iphase++) {

    // Alias the multifluid-coefficients onto a single phase.
    RefCountedPtr<LevelData<EBCellFAB>>       Acoef = RefCountedPtr<LevelData<EBCellFAB>>(new LevelData<EBCellFAB>());
    RefCountedPtr<LevelData<EBFluxFAB>>       Bcoef = RefCountedPtr<LevelData<EBFluxFAB>>(new LevelData<EBFluxFAB>());
    RefCountedPtr<LevelData<BaseIVFAB<Real>>> BcoefIrreg = RefCountedPtr<LevelData<BaseIVFAB<Real>>>(
      new LevelData<BaseIVFAB<Real>>());

    MultifluidAlias::aliasMF(*Acoef, iphase, *a_Acoef);
    MultifluidAlias::aliasMF(*Bcoef, iphase, *a_Bcoef);
    MultifluidAlias::aliasMF(*BcoefIrreg, iphase, *a_BcoefIrreg);

    m_helmOps.at(iphase)->setAcoAndBco(Acoef, Bcoef, BcoefIrreg);
  }

  // Jump BC object also needs to update coefficients.
  m_jumpBC->setBco(a_BcoefIrreg);
}

const RefCountedPtr<LevelData<MFCellFAB>>&
MFHelmholtzOp::getAcoef()
{
  return m_Acoef;
}

const RefCountedPtr<LevelData<MFFluxFAB>>&
MFHelmholtzOp::getBcoef()
{
  return m_Bcoef;
}

const RefCountedPtr<LevelData<MFBaseIVFAB>>&
MFHelmholtzOp::getBcoefIrreg()
{
  return m_BcoefIrreg;
}

void
MFHelmholtzOp::setJump(RefCountedPtr<LevelData<BaseIVFAB<Real>>>& a_jump)
{
  CH_TIME("MFHelmholtzOp::setJump");

  m_jump = a_jump;
}

int
MFHelmholtzOp::refToCoarser()
{
  CH_TIME("MFHelmholtzOp::refToCoarser");

  return m_refToCoar;
}

void
MFHelmholtzOp::setAlphaAndBeta(const Real& a_alpha, const Real& a_beta)
{
  CH_TIME("MFHelmholtzOp::setAlphaAndBeta");

  for (auto& op : m_helmOps) {
    op.second->setAlphaAndBeta(a_alpha, a_beta);
  }
}

void
MFHelmholtzOp::divideByIdentityCoef(LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFHelmholtzOp::divideByIdentityCoef");

  LevelData<EBCellFAB> rhs;

  for (auto& op : m_helmOps) {
    MultifluidAlias::aliasMF(rhs, op.first, a_rhs);

    op.second->divideByIdentityCoef(rhs);
  }
}

void
MFHelmholtzOp::applyOpNoBoundary(LevelData<MFCellFAB>& a_ans, const LevelData<MFCellFAB>& a_phi)
{
  CH_TIME("MFHelmholtzOp::applyOpNoBoundary");

  LevelData<EBCellFAB> ans;
  LevelData<EBCellFAB> phi;

  for (auto& op : m_helmOps) {
    MultifluidAlias::aliasMF(ans, op.first, a_ans);
    MultifluidAlias::aliasMF(phi, op.first, a_phi);

    op.second->applyOpNoBoundary(ans, phi);
  }
}

void
MFHelmholtzOp::fillGrad(const LevelData<MFCellFAB>& a_phi)
{
  CH_TIME("MFHelmholtzOp::fillGrad(LD<MFCellFAB>)");

  LevelData<EBCellFAB> phi;

  for (auto& op : m_helmOps) {
    MultifluidAlias::aliasMF(phi, op.first, a_phi);

    op.second->fillGrad(phi);
  }
}

void
MFHelmholtzOp::getFlux(MFFluxFAB&                  a_flux,
                       const LevelData<MFCellFAB>& a_data,
                       const Box&                  a_grid,
                       const DataIndex&            a_dit,
                       Real                        a_scale)
{
  MayDay::Warning("MFHelmholtzOp::getFlux - not implemented (yet)");
}

void
MFHelmholtzOp::incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, Real a_scale)
{
  CH_TIME("MFHelmholtzOp::incr");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()].plus(a_rhs[dit()], a_scale);
  }
}

void
MFHelmholtzOp::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale)
{
  CH_TIME("MFHelmholtzOp::scale");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()] *= a_scale;
  }
}

void
MFHelmholtzOp::setToZero(LevelData<MFCellFAB>& a_lhs)
{
  CH_TIME("MFHelmholtzOp::setToZero)");

  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    a_lhs[dit()].setVal(0.0);
  }
}

void
MFHelmholtzOp::assign(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFHelmholtzOp::assign");

  if (!(m_copier.isDefined())) {
    m_copier.define(a_rhs.disjointBoxLayout(), a_lhs.disjointBoxLayout(), m_ghostPhi);
  }
  a_rhs.copyTo(a_lhs, m_copier);
}

void
MFHelmholtzOp::assignLocal(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFHelmholtzOp::assignLocal");

  a_rhs.localCopyTo(a_lhs);
}

Real
MFHelmholtzOp::norm(const LevelData<MFCellFAB>& a_lhs, int a_order)
{
  CH_TIME("MFHelmholtzOp::norm");

  Real norm = 0.0;
  for (auto& op : m_helmOps) {
    LevelData<EBCellFAB> lhs;

    MultifluidAlias::aliasMF(lhs, op.first, a_lhs);

    const Real curNorm = op.second->norm(lhs, a_order);

    norm = std::max(norm, curNorm);
  }

  return norm;
}

Real
MFHelmholtzOp::dotProduct(const LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFHelmholtzOp::dotProduct)");

  Real ret = 0.0;

  Real accum = 0.0;
  Real volum = 0.0;
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit) {
    const MFCellFAB& lhs = a_lhs[dit()];

    Real phaseVolume;
    for (int i = 0; i < lhs.numPhases(); i++) {
      const EBCellFAB& data1   = a_lhs[dit()].getPhase(i);
      const EBCellFAB& data2   = a_rhs[dit()].getPhase(i);
      const EBISBox&   ebisbox = data1.getEBISBox();

      if (!ebisbox.isAllCovered()) {
        const Box box = a_lhs.disjointBoxLayout()[dit()];

        accum += EBLevelDataOps::sumKappaDotProduct(phaseVolume,
                                                    data1,
                                                    data2,
                                                    box,
                                                    EBLEVELDATAOPS_ALLVOFS,
                                                    m_mflg.getDomain());
        volum += phaseVolume;
      }
    }
  }

#ifdef CH_MPI
  Real recv;
  MPI_Allreduce(&accum, &recv, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  accum = recv;
  MPI_Allreduce(&volum, &recv, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm);
  volum = recv;
#endif

  if (volum > 0.0) {
    ret = accum / volum;
  }

  return ret;
}

void
MFHelmholtzOp::create(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFHelmholtzOp::create");

  Vector<EBISLayout> layouts;
  Vector<int>        comps;
  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    layouts.push_back(m_mflg.getEBLevelGrid(iphase).getEBISL());
    comps.push_back(m_nComp);
  }

  MFCellFactory cellFact(layouts, comps);
  a_lhs.define(m_mflg.getGrids(), m_nComp, a_rhs.ghostVect(), cellFact);
}

void
MFHelmholtzOp::createCoarser(LevelData<MFCellFAB>& a_coarse, const LevelData<MFCellFAB>& a_fine, bool a_ghosted)
{
  CH_TIME("MFHelmholtzOp::createCoarser");

  Vector<EBISLayout> layouts;
  Vector<int>        comps;
  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    layouts.push_back(m_mflgCoarMG.getEBLevelGrid(iphase).getEBISL());
    comps.push_back(m_nComp);
  }

  MFCellFactory cellFact(layouts, comps);
  a_coarse.define(m_mflgCoarMG.getGrids(), m_nComp, a_fine.ghostVect(), cellFact);
}

void
MFHelmholtzOp::createCoarsened(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const int& a_refRat)
{
  CH_TIME("MFHelmholtzOp::createCoarsened");

  Vector<EBISLayout> layouts;
  Vector<int>        comps;
  for (int iphase = 0; iphase < m_numPhases; iphase++) {
    layouts.push_back(m_mflgCoFi.getEBLevelGrid(iphase).getEBISL());
    comps.push_back(m_nComp);
  }

  MFCellFactory cellFact(layouts, comps);
  a_lhs.define(m_mflgCoFi.getGrids(), m_nComp, a_rhs.ghostVect(), cellFact);
}

void
MFHelmholtzOp::preCond(LevelData<MFCellFAB>& a_corr, const LevelData<MFCellFAB>& a_residual)
{
  CH_TIME("MFHelmholtzOp::preCond");
#if 1
  this->relax(a_corr, a_residual, 40);
#else
  m_jumpBC->resetBC();

  for (auto& op : m_helmOps) {
    LevelData<EBCellFAB> corr;
    LevelData<EBCellFAB> resi;

    MultifluidAlias::aliasMF(corr, op.first, a_corr);
    MultifluidAlias::aliasMF(resi, op.first, a_residual);

    op.second->preCond(corr, resi);
  }
#endif
}

void
MFHelmholtzOp::applyOp(LevelData<MFCellFAB>& a_Lphi, const LevelData<MFCellFAB>& a_phi, bool a_homogeneousPhysBC)
{
  CH_TIME("MFHelmholtzOp::applyOp");

  constexpr bool homogeneousCFBC = true;

  this->applyOp(a_Lphi, a_phi, nullptr, a_homogeneousPhysBC, homogeneousCFBC);
}

void
MFHelmholtzOp::computeOperatorLoads(LevelData<MFCellFAB>& a_phi, TimedDataIterator& a_timeDit)
{
  CH_TIME("MFHelmholtzOp::computeOperatorLoads");

  // TLDR: This routine estimates the time spent in each grid patch for a typical relaxation step. This includes
  //       coarse-fine interpolation, BC matching, and applying the operator.

  LevelData<MFCellFAB> Lphi;
  this->create(Lphi, a_phi);

  for (a_timeDit.reset(); a_timeDit.ok(); ++a_timeDit) {

    // Interpolation time with coarser
    if (m_hasCoar) {
      for (auto& op : m_helmOps) {
        const int iphase = op.first;

        RefCountedPtr<EBMultigridInterpolator>& phaseInterpolator = m_interpolator.getInterpolator(iphase);

        EBCellFAB& phi = (EBCellFAB&)a_phi[a_timeDit()].getPhase(iphase);

        phaseInterpolator->coarseFineInterpH(phi, Interval(m_comp, m_comp), a_timeDit());
      }
    }

    // Matching time
    if (m_multifluid) {
      m_jumpBC->matchBC((*m_jump)[a_timeDit()], a_phi[a_timeDit()], true, a_timeDit());
    }

    // Apply operator application
    for (auto& op : m_helmOps) {
      const Box cellBox = Lphi.disjointBoxLayout()[a_timeDit()];

      const int iphase = op.first;

      EBCellFAB& Lph = Lphi[a_timeDit()].getPhase(iphase);
      EBCellFAB& phi = a_phi[a_timeDit()].getPhase(iphase);

      op.second->applyOp(Lph, phi, cellBox, a_timeDit(), true);
    }
  }
}

void
MFHelmholtzOp::applyOp(LevelData<MFCellFAB>&             a_Lphi,
                       const LevelData<MFCellFAB>&       a_phi,
                       const LevelData<MFCellFAB>* const a_phiCoar,
                       const bool                        a_homogeneousPhysBC,
                       const bool                        a_homogeneousCFBC)
{
  CH_TIME("MFHelmholtzOp::applyOp");

  // We need updated ghost cells since both the operator stencil and the "jump" stencil
  // reach into ghost regions.
  this->exchangeGhost(a_phi);
  this->interpolateCF(a_phi, a_phiCoar, a_homogeneousCFBC);
  this->updateJumpBC(a_phi, a_homogeneousPhysBC);

  // Now apply the operator on each patch.
  const DisjointBoxLayout& dbl = m_mflg.getGrids();
  for (DataIterator dit(dbl); dit.ok(); ++dit) {
    const Box cellBox = dbl[dit()];

    for (auto& op : m_helmOps) {
      const int iphase = op.first;

      // Doing the nasty, but applyOp will only monkey with ghost cells in a_phi.
      EBCellFAB& Lph = (EBCellFAB&)a_Lphi[dit()].getPhase(iphase);
      EBCellFAB& phi = (EBCellFAB&)a_phi[dit()].getPhase(iphase);

      op.second->applyOp(Lph, phi, cellBox, dit(), a_homogeneousPhysBC);
    }
  }
}

void
MFHelmholtzOp::interpolateCF(const LevelData<MFCellFAB>& a_phi,
                             const LevelData<MFCellFAB>* a_phiCoar,
                             const bool                  a_homogeneousCF)
{
  CH_TIME("MFHelmholtzOp::interpolateCF");

  // TLDR: This is a wrapper for interpolating ghost cells on each phase. The user can put a_homogeneousCF = false if he wants inhomogeneous interpolation. This routine
  //       was written so that we avoid calling Multifluid::aliasMF, since that tends to be expensive to call during every smoothing step.

  if (m_hasCoar) {
    if (a_homogeneousCF) {
      // The homogeneous version will be called on every relaxation so we use a format which avoid having to alias data (which can be expensive).
      for (DataIterator dit = a_phi.dataIterator(); dit.ok(); ++dit) {
        for (auto& op : m_helmOps) {
          const int iphase = op.first;

          RefCountedPtr<EBMultigridInterpolator>& phaseInterpolator = m_interpolator.getInterpolator(iphase);

          EBCellFAB& phi = (EBCellFAB&)a_phi[dit()].getPhase(iphase);

          phaseInterpolator->coarseFineInterpH(phi, Interval(m_comp, m_comp), dit());
        }
      }
    }
    else { // Inhomogeneous coarse-fine is not called as often, so performance is not critical here.
      for (auto& op : m_helmOps) {
        LevelData<EBCellFAB> phi;
        LevelData<EBCellFAB> phiCoar;

        // I will call this an error because there's nothing to interpolate from.
        if (a_phiCoar == nullptr) {
          MayDay::Error(
            "MFHelmholtzOp::interpolateCF -- calling inhomogeneousCFInterp with nullptr coarse is an error.");
        }

        MultifluidAlias::aliasMF(phi, op.first, (LevelData<MFCellFAB>&)a_phi);
        MultifluidAlias::aliasMF(phiCoar, op.first, *a_phiCoar);

        op.second->inhomogeneousCFInterp(phi, phiCoar);
      }
    }
  }
}

void
MFHelmholtzOp::residual(LevelData<MFCellFAB>&       a_residual,
                        const LevelData<MFCellFAB>& a_phi,
                        const LevelData<MFCellFAB>& a_rhs,
                        const bool                  a_homogeneousPhysBC)
{
  CH_TIME("MFHelmholtzOp::residual(LD<MFCellFAB>, LD<MFCellFAB>, LD<MFCellFAB>, bool)");

  // Compute a_residual = rhs - L(phi)
  this->applyOp(a_residual, a_phi, a_homogeneousPhysBC);
  this->scale(a_residual, -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void
MFHelmholtzOp::axby(LevelData<MFCellFAB>&       a_lhs,
                    const LevelData<MFCellFAB>& a_x,
                    const LevelData<MFCellFAB>& a_y,
                    const Real                  a,
                    const Real                  b)
{
  CH_TIME("MFHelmholtzOp::axby");

  for (auto& op : m_helmOps) {
    LevelData<EBCellFAB> lhs;
    LevelData<EBCellFAB> x;
    LevelData<EBCellFAB> y;

    MultifluidAlias::aliasMF(lhs, op.first, a_lhs);
    MultifluidAlias::aliasMF(x, op.first, a_x);
    MultifluidAlias::aliasMF(y, op.first, a_y);

    op.second->axby(lhs, x, y, a, b);
  }
}

void
MFHelmholtzOp::updateJumpBC(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneousPhysBC)
{
  CH_TIME("MFHelmholtzOp::updateJumpBC");

  m_jumpBC->matchBC(*m_jump, a_phi, a_homogeneousPhysBC);
}

void
MFHelmholtzOp::exchangeGhost(const LevelData<MFCellFAB>& a_phi) const
{
  CH_TIME("MFHelmholtzOp::exchangeGhost");

  LevelData<MFCellFAB>& phi = (LevelData<MFCellFAB>&)a_phi;

  phi.exchange();
}

void
MFHelmholtzOp::relax(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, int a_iterations)
{
  CH_TIME("MFHelmholtzOp::relax");

  // This function performs relaxation. The user can switch between various kernels.

  switch (m_smoother) {
  case Smoother::PointJacobi: {
    this->relaxPointJacobi(a_correction, a_residual, a_iterations);

    break;
  }
  case Smoother::GauSaiRedBlack: {
    this->relaxGSRedBlack(a_correction, a_residual, a_iterations);

    break;
  }
  case Smoother::GauSaiMultiColor: {
    this->relaxGSMultiColor(a_correction, a_residual, a_iterations);

    break;
  }
  default: {
    MayDay::Error("MFHelmholtzOp::relax - bogus relaxation method requested");

    break;
  }
  };
}

void
MFHelmholtzOp::relaxPointJacobi(LevelData<MFCellFAB>&       a_correction,
                                const LevelData<MFCellFAB>& a_residual,
                                const int                   a_iterations)
{
  CH_TIME("MFHelmholtzOp::relaxPointJacobi");

  // TLDR: This function performs point Jacobi relaxation in the form phi^(k+1) = phi^k - (res - L(phi))/|diag(L)|. Here, diag(L) is captured
  //       in m_relCoef. For performance integration, EBHelmholtzOp has a public function for the kernel.

  LevelData<MFCellFAB> Lcorr;
  this->create(Lcorr, a_correction);

  const DisjointBoxLayout& dbl = m_mflg.getGrids();

  constexpr bool homogeneousCFBC   = true;
  constexpr bool homogeneousPhysBC = true;

  for (int i = 0; i < a_iterations; i++) {

    // Fill/interpolate ghost cells and match the BC.
    this->exchangeGhost(a_correction);
    this->interpolateCF(a_correction, nullptr, homogeneousCFBC);
    this->updateJumpBC(a_correction, homogeneousPhysBC);

    // Do relaxation on each patch.
    for (DataIterator dit(dbl); dit.ok(); ++dit) {
      const Box cellBox = dbl[dit()];

      for (auto& op : m_helmOps) {
        const int iphase = op.first;

        EBCellFAB&       Lph = Lcorr[dit()].getPhase(iphase);
        EBCellFAB&       phi = a_correction[dit()].getPhase(iphase);
        const EBCellFAB& res = a_residual[dit()].getPhase(iphase);

        op.second->pointJacobiKernel(Lph, phi, res, cellBox, dit());
      }
    }
  }
}

void
MFHelmholtzOp::relaxGSRedBlack(LevelData<MFCellFAB>&       a_correction,
                               const LevelData<MFCellFAB>& a_residual,
                               const int                   a_iterations)
{
  CH_TIME("MFHelmholtzOp::relaxGSRedBlack");

  // TLDR: This function performs red-black Gauss-Seidel relaxation. As always, this occurs in the form phi^(k+1) = phi^k - (res - L(phi))/|diag(L)| but
  //       for a red-black update pattern:
  //
  //       For performance integration, this calls the EBHelmholtzOp red-black kernel directly.

  LevelData<MFCellFAB> Lcorr;
  this->create(Lcorr, a_correction);

  const DisjointBoxLayout& dbl = m_mflg.getGrids();

  constexpr bool homogeneousCFBC   = true;
  constexpr bool homogeneousPhysBC = true;

  for (int i = 0; i < a_iterations; i++) {
    for (int redBlack = 0; redBlack <= 1; redBlack++) {

      // Fill/interpolate ghost cells and match the BC.
      this->exchangeGhost(a_correction);
      this->interpolateCF(a_correction, nullptr, homogeneousCFBC);
      this->updateJumpBC(a_correction, homogeneousPhysBC);

      // Do relaxation on each patch.
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const Box cellBox = dbl[dit()];

        for (auto& op : m_helmOps) {
          const int iphase = op.first;

          EBCellFAB&       Lph = Lcorr[dit()].getPhase(iphase);
          EBCellFAB&       phi = a_correction[dit()].getPhase(iphase);
          const EBCellFAB& res = a_residual[dit()].getPhase(iphase);

          op.second->gauSaiRedBlackKernel(Lph, phi, res, cellBox, dit(), redBlack);
        }
      }
    }
  }
}

void
MFHelmholtzOp::relaxGSMultiColor(LevelData<MFCellFAB>&       a_correction,
                                 const LevelData<MFCellFAB>& a_residual,
                                 const int                   a_iterations)
{
  CH_TIME("MFHelmholtzOp::relaxGSMultiColor");

  // TLDR: This function performs multi-colored Gauss-Seidel relaxation. As always, this occurs in the form phi^(k+1) = phi^k - (res - L(phi))/|diag(L)| but
  //       using more colors than just red-black. The update pattern here cycles through quadrants/octants in 2D/3D. This is just like red-black except that
  //       we have four/eight colors in 2D/3D.
  //
  //       For performance integration, this calls the EBHelmholtzOp multi-color kernel directly.

  LevelData<MFCellFAB> Lcorr;
  this->create(Lcorr, a_correction);

  const DisjointBoxLayout& dbl = m_mflg.getGrids();

  constexpr bool homogeneousCFBC   = true;
  constexpr bool homogeneousPhysBC = true;

  for (int i = 0; i < a_iterations; i++) {

    for (int icolor = 0; icolor < m_colors.size(); icolor++) {

      // Fill/interpolate ghost cells and match the BC.
      this->exchangeGhost(a_correction);
      this->interpolateCF(a_correction, nullptr, homogeneousCFBC);
      this->updateJumpBC(a_correction, homogeneousPhysBC);

      // Do relaxation on each patch
      for (DataIterator dit(dbl); dit.ok(); ++dit) {
        const Box cellBox = dbl[dit()];

        for (auto& op : m_helmOps) {
          const int iphase = op.first;

          EBCellFAB&       Lph = Lcorr[dit()].getPhase(iphase);
          EBCellFAB&       phi = a_correction[dit()].getPhase(iphase);
          const EBCellFAB& res = a_residual[dit()].getPhase(iphase);

          op.second->gauSaiMultiColorKernel(Lph, phi, res, cellBox, dit(), m_colors[icolor]);
        }
      }
    }
  }
}

void
MFHelmholtzOp::restrictResidual(LevelData<MFCellFAB>&       a_resCoar,
                                LevelData<MFCellFAB>&       a_phi,
                                const LevelData<MFCellFAB>& a_rhs)
{
  CH_TIME("MFHelmholtzOp::restrictResidual");

  constexpr bool homogeneousPhysBC = true;

  // EBHelmholtzOp::restrictResidual will call applyOp so we need to update the boundary condition
  // on multiphase cells first.
  this->exchangeGhost(a_phi);
  this->updateJumpBC(a_phi, homogeneousPhysBC);

  for (auto& op : m_helmOps) {
    LevelData<EBCellFAB> resCoar;
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> rhs;

    MultifluidAlias::aliasMF(resCoar, op.first, a_resCoar);
    MultifluidAlias::aliasMF(phi, op.first, a_phi);
    MultifluidAlias::aliasMF(rhs, op.first, a_rhs);

    op.second->restrictResidual(resCoar, phi, rhs);
  }
}

void
MFHelmholtzOp::prolongIncrement(LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_correctCoarse)
{
  CH_TIME("MFHelmholtzOp::prolongIncrement");

  for (auto& op : m_helmOps) {
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> correctCoarse;

    MultifluidAlias::aliasMF(phi, op.first, a_phi);
    MultifluidAlias::aliasMF(correctCoarse, op.first, a_correctCoarse);

    op.second->prolongIncrement(phi, correctCoarse);
  }
}

void
MFHelmholtzOp::AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
                                 const LevelData<MFCellFAB>& a_correction,
                                 const LevelData<MFCellFAB>& a_coarseCorrection)
{
  CH_TIME("MFHelmholtzOp::AMRUpdateResidual");

  constexpr bool homogeneousCFBC   = false;
  constexpr bool homogeneousPhysBC = true;

  // Need to update BC first!
  this->exchangeGhost(a_correction);
  this->interpolateCF((LevelData<MFCellFAB>&)a_correction, &a_coarseCorrection, homogeneousCFBC);
  this->updateJumpBC(a_correction, homogeneousPhysBC);

  for (auto& op : m_helmOps) {
    LevelData<EBCellFAB> residual;
    LevelData<EBCellFAB> correction;
    LevelData<EBCellFAB> coarseCorrection;

    MultifluidAlias::aliasMF(residual, op.first, a_residual);
    MultifluidAlias::aliasMF(correction, op.first, a_correction);
    MultifluidAlias::aliasMF(coarseCorrection, op.first, a_coarseCorrection);

    // Don't need to update ghost cells or exchange again.
    op.second->turnOffCFInterp();
    op.second->turnOffCoarsening();
    op.second->turnOffExchange();

    op.second->AMRUpdateResidual(residual, correction, coarseCorrection);

    op.second->turnOnCFInterp();
    op.second->turnOnCoarsening();
    op.second->turnOnExchange();
  }
}

void
MFHelmholtzOp::AMRRestrict(LevelData<MFCellFAB>&       a_residualCoarse,
                           const LevelData<MFCellFAB>& a_residual,
                           const LevelData<MFCellFAB>& a_correction,
                           const LevelData<MFCellFAB>& a_coarseCorrection,
                           bool                        a_skip_res)
{
  CH_TIME("MFHelmholtzOp::AMRRestrict");

  constexpr bool homogeneousCFBC   = false;
  constexpr bool homogeneousPhysBC = true;

  this->exchangeGhost(a_correction);
  this->interpolateCF((LevelData<MFCellFAB>&)a_correction, &a_coarseCorrection, homogeneousCFBC);
  this->updateJumpBC(a_correction, homogeneousPhysBC);

  for (auto& op : m_helmOps) {
    LevelData<EBCellFAB> residualCoarse;
    LevelData<EBCellFAB> residual;
    LevelData<EBCellFAB> correction;
    LevelData<EBCellFAB> coarseCorrection;

    MultifluidAlias::aliasMF(residualCoarse, op.first, a_residualCoarse);
    MultifluidAlias::aliasMF(residual, op.first, a_residual);
    MultifluidAlias::aliasMF(correction, op.first, a_correction);
    MultifluidAlias::aliasMF(coarseCorrection, op.first, a_coarseCorrection);

    // Don't need to update ghost cells or exchange again.
    op.second->turnOffCFInterp();
    op.second->turnOffCoarsening();
    op.second->turnOffExchange();

    op.second->AMRRestrict(residualCoarse, residual, correction, coarseCorrection, a_skip_res);

    op.second->turnOnCFInterp();
    op.second->turnOnCoarsening();
    op.second->turnOnExchange();
  }
}

void
MFHelmholtzOp::AMRProlong(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_coarseCorrection)
{
  CH_TIME("MFHelmholtzOp::AMRProlong");

  for (auto& op : m_helmOps) {
    LevelData<EBCellFAB> correction;
    LevelData<EBCellFAB> coarseCorrection;

    MultifluidAlias::aliasMF(correction, op.first, a_correction);
    MultifluidAlias::aliasMF(coarseCorrection, op.first, a_coarseCorrection);

    op.second->AMRProlong(correction, coarseCorrection);
  }
}

void
MFHelmholtzOp::AMRResidual(LevelData<MFCellFAB>&             a_residual,
                           const LevelData<MFCellFAB>&       a_phiFine,
                           const LevelData<MFCellFAB>&       a_phi,
                           const LevelData<MFCellFAB>&       a_phiCoar,
                           const LevelData<MFCellFAB>&       a_rhs,
                           bool                              a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<MFCellFAB>>* a_finerOp)
{
  CH_TIME("MFHelmholtzOp::AMRResidual");

  // Make residual = a_rhs - L(phi)
  this->AMROperator(a_residual, a_phiFine, a_phi, a_phiCoar, a_homogeneousPhysBC, a_finerOp);

  this->scale(a_residual, -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void
MFHelmholtzOp::AMRResidualNF(LevelData<MFCellFAB>&       a_residual,
                             const LevelData<MFCellFAB>& a_phi,
                             const LevelData<MFCellFAB>& a_phiCoar,
                             const LevelData<MFCellFAB>& a_rhs,
                             bool                        a_homogeneousPhysBC)
{
  CH_TIME("MFHelmholtzOp::AMRResidualNF");

  // Make residual = a_rhs - L(phi)
  this->AMROperatorNF(a_residual, a_phi, a_phiCoar, a_homogeneousPhysBC);

  this->scale(a_residual, -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void
MFHelmholtzOp::AMRResidualNC(LevelData<MFCellFAB>&             a_residual,
                             const LevelData<MFCellFAB>&       a_phiFine,
                             const LevelData<MFCellFAB>&       a_phi,
                             const LevelData<MFCellFAB>&       a_rhs,
                             bool                              a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<MFCellFAB>>* a_finerOp)
{
  CH_TIME("MFHelmholtzOp::AMRResidualNC");

  this->AMROperatorNC(a_residual, a_phiFine, a_phi, a_homogeneousPhysBC, a_finerOp);

  this->scale(a_residual, -1.0);
  this->incr(a_residual, a_rhs, 1.0);
}

void
MFHelmholtzOp::AMROperatorNF(LevelData<MFCellFAB>&       a_Lphi,
                             const LevelData<MFCellFAB>& a_phi,
                             const LevelData<MFCellFAB>& a_phiCoar,
                             bool                        a_homogeneousPhysBC)
{
  CH_TIME("MFHelmholtzOp::AMROperatorNF");

  constexpr bool homogeneousCFBC = false;

  // Note; There is no coarse-fine interpolation here because that will have been
  // done by AMROperator (which is called before this routine).
  this->exchangeGhost(a_phi);
  this->updateJumpBC(a_phi, a_homogeneousPhysBC);

  for (auto& op : m_helmOps) {
    LevelData<EBCellFAB> Lphi;
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> phiCoar;

    MultifluidAlias::aliasMF(Lphi, op.first, a_Lphi);
    MultifluidAlias::aliasMF(phi, op.first, a_phi);
    MultifluidAlias::aliasMF(phiCoar, op.first, a_phiCoar);

    op.second->turnOffCFInterp();
    op.second->turnOffCoarsening();
    op.second->turnOffExchange();

    op.second->AMROperatorNF(Lphi, phi, phiCoar, a_homogeneousPhysBC);

    op.second->turnOnCFInterp();
    op.second->turnOnCoarsening();
    op.second->turnOnExchange();
  }
}

void
MFHelmholtzOp::AMROperatorNC(LevelData<MFCellFAB>&             a_Lphi,
                             const LevelData<MFCellFAB>&       a_phiFine,
                             const LevelData<MFCellFAB>&       a_phi,
                             bool                              a_homogeneousPhysBC,
                             AMRLevelOp<LevelData<MFCellFAB>>* a_finerOp)
{
  CH_TIME("MFHelmholtzOp::AMROperatorNC");

  if (m_hasFine) {
    MFHelmholtzOp* finerOp = (MFHelmholtzOp*)a_finerOp;

    for (auto& op : finerOp->m_helmOps) {
      LevelData<EBCellFAB> phi;
      LevelData<EBCellFAB> phiFine;

      MultifluidAlias::aliasMF(phi, op.first, a_phi);
      MultifluidAlias::aliasMF(phiFine, op.first, a_phiFine);

      op.second->coarsenCell(phi, phiFine);
    }
  }

  this->exchangeGhost(a_phi);
  this->updateJumpBC(a_phi, a_homogeneousPhysBC);

  for (auto& op : m_helmOps) {
    LevelData<EBCellFAB> Lphi;
    LevelData<EBCellFAB> phiFine;
    LevelData<EBCellFAB> phi;

    MultifluidAlias::aliasMF(Lphi, op.first, a_Lphi);
    MultifluidAlias::aliasMF(phiFine, op.first, a_phiFine);
    MultifluidAlias::aliasMF(phi, op.first, a_phi);

    MFHelmholtzOp* finerOp = (MFHelmholtzOp*)(a_finerOp);

    // Don't need to update ghost cells again, coarsen, or exchange data.
    op.second->turnOffCFInterp();
    op.second->turnOffCoarsening();
    op.second->turnOffExchange();

    op.second->AMROperatorNC(Lphi, phiFine, phi, a_homogeneousPhysBC, (finerOp->m_helmOps).at(op.first));

    op.second->turnOnCFInterp();
    op.second->turnOnCoarsening();
    op.second->turnOnExchange();
  }
}

void
MFHelmholtzOp::AMROperator(LevelData<MFCellFAB>&             a_Lphi,
                           const LevelData<MFCellFAB>&       a_phiFine,
                           const LevelData<MFCellFAB>&       a_phi,
                           const LevelData<MFCellFAB>&       a_phiCoar,
                           const bool                        a_homogeneousPhysBC,
                           AMRLevelOp<LevelData<MFCellFAB>>* a_finerOp)
{
  CH_TIME("MFHelmholtzOp::AMROperator");

  constexpr bool homogeneousCFBC = false;

  if (m_hasFine) {
    MFHelmholtzOp* finerOp = (MFHelmholtzOp*)a_finerOp;

    for (auto& op : finerOp->m_helmOps) {
      LevelData<EBCellFAB> phi;
      LevelData<EBCellFAB> phiFine;

      MultifluidAlias::aliasMF(phi, op.first, a_phi);
      MultifluidAlias::aliasMF(phiFine, op.first, a_phiFine);

      op.second->coarsenCell(phi, phiFine);
    }
  }

  this->exchangeGhost(a_phi);
  this->updateJumpBC(a_phi, a_homogeneousPhysBC);

  for (auto& op : m_helmOps) {
    LevelData<EBCellFAB> Lphi;
    LevelData<EBCellFAB> phiFine;
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> phiCoar;

    MultifluidAlias::aliasMF(Lphi, op.first, a_Lphi);
    MultifluidAlias::aliasMF(phiFine, op.first, a_phiFine);
    MultifluidAlias::aliasMF(phi, op.first, a_phi);
    MultifluidAlias::aliasMF(phiCoar, op.first, a_phiCoar);

    MFHelmholtzOp* finerOp = (MFHelmholtzOp*)(a_finerOp);

    // Don't need to update ghost cells again, coarsen, or exchange data. Our ability to turn off coarse-fine interpolation comes from
    // the fact that EBHelmholtzOp will interpolate the ghost cells on the finer level during the reflux stage. When we enter this routine
    // we already have updated our ghost cells!
    op.second->turnOffCFInterp();
    op.second->turnOffCoarsening();
    op.second->turnOffExchange();

    op.second->AMROperator(Lphi, phiFine, phi, phiCoar, a_homogeneousPhysBC, (finerOp->m_helmOps).at(op.first));

    op.second->turnOnCFInterp();
    op.second->turnOnCoarsening();
    op.second->turnOnExchange();
  }
}

#include <CD_NamespaceFooter.H>
