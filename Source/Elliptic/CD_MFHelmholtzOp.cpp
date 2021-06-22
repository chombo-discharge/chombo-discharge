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

// Dev includes
#include <CD_EBHelmholtzDirichletEBBC.H>
#include <CD_EBHelmholtzDirichletDomainBC.H>

// Our includes
#include <CD_MFHelmholtzOp.H>
#include <CD_MultifluidAlias.H>
#include <CD_NamespaceHeader.H>

constexpr int MFHelmholtzOp::m_comp;
constexpr int MFHelmholtzOp::m_nComp;

MFHelmholtzOp::MFHelmholtzOp(){
  MayDay::Abort("MFHelmholtzOp - weak construction is not allowed");
}

MFHelmholtzOp::MFHelmholtzOp(const MFLevelGrid&                               a_mflgFine,
			     const MFLevelGrid&                               a_mflg,
			     const MFLevelGrid&                               a_mflgCoFi,
			     const MFLevelGrid&                               a_mflgCoar,
			     const MFLevelGrid&                               a_mflgCoarMG,
			     const MFMultigridInterpolator&                   a_interpolator,
			     const MFFluxReg&                                 a_fluxReg,
			     const MFCoarAve&                                 a_coarAve,
			     const RefCountedPtr<EBHelmholtzDomainBCFactory>& a_domainBcFactory,
			     const RealVect&                                  a_probLo,
			     const Real&                                      a_dx,
			     const int&                                       a_refToCoar,
			     const bool&                                      a_hasFine,
			     const bool&                                      a_hasCoar,			     
			     const bool&                                      a_hasMGObjects,
			     const Real&                                      a_alpha,
			     const Real&                                      a_beta,
			     const RefCountedPtr<LevelData<MFCellFAB> >&      a_Acoef,
			     const RefCountedPtr<LevelData<MFFluxFAB> >&      a_Bcoef,
			     const RefCountedPtr<LevelData<MFBaseIVFAB> >&    a_BcoefIrreg,
			     const IntVect&                                   a_ghostPhi,
			     const IntVect&                                   a_ghostRhs,
			     const int&                                       a_ebbcOrder,
			     const int&                                       a_jumpOrder,
			     const RelaxType&                                 a_relaxType){
  CH_TIME("MFHelmholtzOp::MFHelmholtzOp");
  m_debug = false;
  ParmParse pp ("MFHelmholtzOp");
  pp.query("debug", m_debug);
  //if(m_debug) pout() << "MFHelmholtzOp::MFHelmholtzOp -- begin" << endl;

  m_mflg         = a_mflg;
  m_numPhases    = m_mflg.numPhases();
  m_multifluid   = m_numPhases > 1;
  m_hasMGObjects = a_hasMGObjects;
  m_refToCoar    = a_refToCoar;

  if(a_hasCoar){
    m_refToCoar = a_refToCoar;
    m_mflgCoFi = a_mflgCoFi;
    m_mflgCoar = a_mflgCoar;
  }

  if(m_hasMGObjects){
    m_mflgCoarMG   = a_mflgCoarMG;
  }

  // Instantiate jump bc object.
  m_jumpBC = RefCountedPtr<JumpBC> (new JumpBC(m_mflg, *a_BcoefIrreg, a_dx, a_jumpOrder, a_jumpOrder, a_jumpOrder));

  // Make the operators on eachphase.
  for (int iphase = 0; iphase < m_numPhases; iphase++){

    EBLevelGrid dummy;
    EBLevelGrid eblgFine    = a_hasFine      ? a_mflgFine.getEBLevelGrid(iphase)   : dummy;
    EBLevelGrid eblg        = a_mflg.getEBLevelGrid(iphase);
    EBLevelGrid eblgCoFi    = a_hasCoar      ? a_mflgCoFi.getEBLevelGrid(iphase)   : dummy;
    EBLevelGrid eblgCoar    = a_hasCoar      ? a_mflgCoar.getEBLevelGrid(iphase)   : dummy;
    EBLevelGrid eblgCoarMG  = a_hasMGObjects ? a_mflgCoarMG.getEBLevelGrid(iphase) : dummy;


    // Development code. Give the operators Dirichlet EBBCs. This will be replaced by other stuff later. 
    auto domainBC = RefCountedPtr<EBHelmholtzDirichletDomainBC> (new EBHelmholtzDirichletDomainBC());
    auto ebBC     = RefCountedPtr<EBHelmholtzDirichletEBBC>     (new EBHelmholtzDirichletEBBC());

    domainBC->setValue(-1.0);
    ebBC->setValue(1.0);
    ebBC->setOrder(1);
    ebBC->setWeight(1);

    // Alias the multifluid-coefficients onto a single phase. 
    RefCountedPtr<LevelData<EBCellFAB> >        Acoef       = RefCountedPtr<LevelData<EBCellFAB> >        (new LevelData<EBCellFAB>());
    RefCountedPtr<LevelData<EBFluxFAB> >        Bcoef       = RefCountedPtr<LevelData<EBFluxFAB> >        (new LevelData<EBFluxFAB>());
    RefCountedPtr<LevelData<BaseIVFAB<Real> > > BcoefIrreg  = RefCountedPtr<LevelData<BaseIVFAB<Real> > > (new LevelData<BaseIVFAB<Real> >());

    MultifluidAlias::aliasMF(*Acoef,      iphase, *a_Acoef);
    MultifluidAlias::aliasMF(*Bcoef,      iphase, *a_Bcoef);
    MultifluidAlias::aliasMF(*BcoefIrreg, iphase, *a_BcoefIrreg);

    RefCountedPtr<EBHelmholtzOp> oper = RefCountedPtr<EBHelmholtzOp> (new EBHelmholtzOp(eblgFine,
											eblg,
											eblgCoFi,
											eblgCoar,
											eblgCoarMG,
											a_interpolator.getInterpolator(iphase),
											a_fluxReg.getFluxRegPointer(iphase),
											a_coarAve.getAveOp(iphase),
											domainBC,
											ebBC,
											a_probLo,
											a_dx,
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
											a_relaxType));


    m_helmOps.emplace(iphase, oper);
  }

  //if(m_debug) pout() << "MFHelmholtzOp::MFHelmholtzOp -- end" << endl;
}

MFHelmholtzOp::~MFHelmholtzOp(){
  CH_TIME("MFHelmholtzOp::~MFHelmholtzOp");
  //if(m_debug) pout() << "MFHelmholtzOp::~MFHelmholtzOp(begin)" << endl;
  //if(m_debug) pout() << "MFHelmholtzOp::~MFHelmholtzOp(end)" << endl;
}

void MFHelmholtzOp::setJump(const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_jump){
  CH_TIME("MFHelmholtzOp::setJump");
  //if(m_debug) pout() << "MFHelmholtzOp::setJump(begin)" << endl;
  m_jump = a_jump;
  //if(m_debug) pout() << "MFHelmholtzOp::setJump(end)" << endl;
}

int MFHelmholtzOp::refToCoarser() {
  CH_TIME("MFHelmholtzOp::refToCoarser");
  //if(m_debug) pout() << "MFHelmholtzOp::refToCoarser(begin)" << endl;
  //if(m_debug) pout() << "MFHelmholtzOp::refToCoarser(end)" << endl;
  return m_refToCoar;
}

unsigned int MFHelmholtzOp::orderOfAccuracy(void) const {
  CH_TIME("MFHelmholtzOp::orderOfAccuracy");
  //if(m_debug) pout() << "MFHelmholtzOp::orderOfAccuracy(begin)" << endl;
  //if(m_debug) pout() << "MFHelmholtzOp::orderOfAccuracy(end)" << endl;
  return 99;
}

void MFHelmholtzOp::enforceCFConsistency(LevelData<MFCellFAB>& a_coarCorr, const LevelData<MFCellFAB>& a_fineCorr) {
  CH_TIME("MFHelmholtzOp::enforceCFConsistency");
  //if(m_debug) pout() << "MFHelmholtzOp::enforceCFConsistency(begin)" << endl;
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> coarCorr;
    LevelData<EBCellFAB> fineCorr;
    
    MultifluidAlias::aliasMF(coarCorr, op.first, a_coarCorr);
    MultifluidAlias::aliasMF(fineCorr, op.first, a_fineCorr);

    op.second->enforceCFConsistency(coarCorr, fineCorr);
  }
  //if(m_debug) pout() << "MFHelmholtzOp::enforceCFConsistency(end)" << endl;
}

void MFHelmholtzOp::incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, Real a_scale){
  CH_TIME("MFHelmholtzOp::incr");
  //if(m_debug) pout() << "MFHelmholtzOp::incr(begin)" << endl;
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()].plus(a_rhs[dit()], a_scale);
  }
  //if(m_debug) pout() << "MFHelmholtzOp::incr(end)" << endl;
}

void MFHelmholtzOp::scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale) {
  CH_TIME("MFHelmholtzOp::scale");
  //if(m_debug) pout() << "MFHelmholtzOp::scale(begin)" << endl;
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()] *= a_scale;
  }
  //if(m_debug) pout() << "MFHelmholtzOp::scale(end)" << endl;
}

void MFHelmholtzOp::setToZero(LevelData<MFCellFAB>& a_lhs) {
  CH_TIME("MFHelmholtzOp::setToZero");
  //if(m_debug) pout() << "MFHelmholtzOp::setToZero(begin)" << endl;
#if 0
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    a_lhs[dit()].setVal(0.0);
  }
#else
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> lhs;

    MultifluidAlias::aliasMF(lhs, op.first, a_lhs);
    op.second->setToZero(lhs);
  }
#endif
  //if(m_debug) pout() << "MFHelmholtzOp::setToZero(end)" << endl;
}

void MFHelmholtzOp::assign(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs) {
  CH_TIME("MFHelmholtzOp::assign");
  //if(m_debug) pout() << "MFHelmholtzOp::setToZero(assign)" << endl;
  //  a_rhs.copyTo(a_lhs);
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> lhs;
    LevelData<EBCellFAB> rhs;

    MultifluidAlias::aliasMF(lhs, op.first, a_lhs);
    MultifluidAlias::aliasMF(rhs, op.first, a_rhs);

    op.second->assign(lhs, rhs);
  }

  // for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
  //   a_lhs[dit()].setVal(0.0);
  //   a_lhs[dit()] += a_rhs[dit()];
  // }
  //if(m_debug) pout() << "MFHelmholtzOp::setToZero(end)" << endl;
}

Real MFHelmholtzOp::norm(const LevelData<MFCellFAB>& a_lhs, int a_order){
  CH_TIME("MFHelmholtzOp::norm");
  //if(m_debug) pout() << "MFHelmholtzOp::norm(begin)" << endl;
  Real norm = 0.0;
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> lhs;

    MultifluidAlias::aliasMF(lhs, op.first, a_lhs);

    const Real curNorm = op.second->norm(lhs, a_order); 
    norm = std::max(norm, curNorm);
  }
  //if(m_debug) pout() << "MFHelmholtzOp::norm(end)" << endl;
  return norm;
}

Real MFHelmholtzOp::dotProduct(const LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs) {
  CH_TIME("MFHelmholtzOp::dotProduct");
  //if(m_debug) pout() << "MFHelmholtzOp::dotProduct(begin)" << endl;
  Real ret = 0.0;

  Real accum = 0.0;
  Real volum = 0.0;
  for (DataIterator dit = a_lhs.dataIterator(); dit.ok(); ++dit){
    const MFCellFAB& lhs = a_lhs[dit()];

    Real phaseVolume;
    for (int i = 0; i < lhs.numPhases(); i++){
      const EBCellFAB& data1 = a_lhs[dit()].getPhase(i);
      const EBCellFAB& data2 = a_rhs[dit()].getPhase(i);
      const Box box = a_lhs.disjointBoxLayout()[dit()];

      accum += EBLevelDataOps::sumKappaDotProduct(phaseVolume, data1, data2, box, EBLEVELDATAOPS_ALLVOFS, m_mflg.getDomain());
      volum += phaseVolume;
    }
  }

#ifdef CH_MPI
  Real recv;
  MPI_Allreduce(&accum, &recv, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm); accum = recv;
  MPI_Allreduce(&volum, &recv, 1, MPI_CH_REAL, MPI_SUM, Chombo_MPI::comm); volum = recv;
#endif


  if(volum > 0.0){
    ret = accum/volum;
  }

  //if(m_debug) pout() << "MFHelmholtzOp::dotProduct(end)" << endl;

  return ret;
}

void MFHelmholtzOp::create(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs){
  CH_TIME("MFHelmholtzOp::create");
  //if(m_debug) pout() << "MFHelmholtzOp::create(begin) on level = " << m_mflg.getDomain() << endl;

  Vector<EBISLayout> layouts;
  Vector<int> comps;
  for (int iphase = 0; iphase < m_numPhases; iphase++){
    layouts.push_back(m_mflg.getEBLevelGrid(iphase).getEBISL());
    comps.push_back(m_nComp);
  }

  MFCellFactory cellFact(layouts, comps);
  a_lhs.define(m_mflg.getGrids(), m_nComp, a_rhs.ghostVect(), cellFact);

  //if(m_debug) pout() << "MFHelmholtzOp::create(end)" << endl;
}

void MFHelmholtzOp::createCoarser(LevelData<MFCellFAB>& a_coarse, const LevelData<MFCellFAB>& a_fine, bool a_ghosted) {
  CH_TIME("MFHelmholtzOp::createCoarser");
  //if(m_debug) pout() << "MFHelmholtzOp::createCoarser(begin) on level = " << m_mflg.getDomain() << endl;
  
  Vector<EBISLayout> layouts;
  Vector<int> comps;
  for (int iphase = 0; iphase < m_numPhases; iphase++){
    layouts.push_back(m_mflgCoarMG.getEBLevelGrid(iphase).getEBISL());
    comps.push_back(m_nComp);
  }

  MFCellFactory cellFact(layouts, comps);
  a_coarse.define(m_mflgCoarMG.getGrids(), m_nComp, a_fine.ghostVect(), cellFact);
  //if(m_debug) pout() << "MFHelmholtzOp::createCoarser(end)" << endl;
}

void MFHelmholtzOp::createCoarsened(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const int& a_refRat) {
  CH_TIME("MFHelmholtzOp::createCoarsened");
  //if(m_debug) pout() << "MFHelmholtzOp::createCoarsened(begin)" << endl;

  Vector<EBISLayout> layouts;
  Vector<int> comps;
  for (int iphase = 0; iphase < m_numPhases; iphase++){
    layouts.push_back(m_mflgCoFi.getEBLevelGrid(iphase).getEBISL());
    comps.push_back(m_nComp);
  }

  MFCellFactory cellFact(layouts, comps);
  a_lhs.define(m_mflgCoFi.getGrids(), m_nComp, a_rhs.ghostVect(), cellFact);
  
  //if(m_debug) pout() << "MFHelmholtzOp::createCoarsened(end)" << endl;
}

void MFHelmholtzOp::preCond(LevelData<MFCellFAB>& a_corr, const LevelData<MFCellFAB>& a_residual) {
  CH_TIME("MFHelmholtzOp::preCond");
  //if(m_debug) pout() << "MFHelmholtzOp::preCond(begin)" << endl;
  //  this->relax(a_corr, a_residual, 40);
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> corr;
    LevelData<EBCellFAB> resi;

    MultifluidAlias::aliasMF(corr, op.first, a_corr);
    MultifluidAlias::aliasMF(resi, op.first, a_residual);

    op.second->preCond(corr, resi);
  }
  //if(m_debug) pout() << "MFHelmholtzOp::preCond(end)" << endl;
}

void MFHelmholtzOp::applyOp(LevelData<MFCellFAB>& a_Lphi, const LevelData<MFCellFAB>& a_phi, bool a_homogeneousPhysBC) {
  CH_TIME("MFHelmholtzOp::applyOp(small)");
  //if(m_debug) pout() << "MFHelmholtzOp::applyOp(begin)" << endl;
  this->applyOp(a_Lphi, a_phi, nullptr, a_homogeneousPhysBC, true);
  //if(m_debug) pout() << "MFHelmholtzOp::applyOp(end)" << endl;
}

void MFHelmholtzOp::applyOp(LevelData<MFCellFAB>&             a_Lphi,
			    const LevelData<MFCellFAB>&       a_phi,
			    const LevelData<MFCellFAB>* const a_phiCoar,
			    const bool                        a_homogeneousPhysBC,
			    const bool                        a_homogeneousCFBC){
  CH_TIME("MFHelmholtzOp::applyOp(big)");
  //if(m_debug) pout() << "MFHelmholtzOp::applyOp(big, begin)" << endl;
  
  // We MUST have updated ghost cells before applying the matching conditions. 
   this->interpolateCF(a_phi, a_phiCoar, a_homogeneousCFBC);
   //this->updateJumpBC(a_phi, a_homogeneousPhysBC);

  // Apply the operator on each level. 
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> Lphi;
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> phiCoar;

    MultifluidAlias::aliasMF(Lphi, op.first, a_Lphi);
    MultifluidAlias::aliasMF(phi,  op.first, a_phi);

    // applyOp usually interpolates ghost cells, but we did that above in interpolateCF(...). 
    op.second->turnOffBCs();
    if(a_phiCoar != nullptr){
      MultifluidAlias::aliasMF(phiCoar, op.first, *a_phiCoar);
      op.second->applyOp(Lphi, phi, &phiCoar, a_homogeneousPhysBC, a_homogeneousCFBC);
    }
    else{
      op.second->applyOp(Lphi, phi, nullptr, a_homogeneousPhysBC, a_homogeneousCFBC);
    }
    op.second->turnOnBCs();
  }
  //if(m_debug) pout() << "MFHelmholtzOp::applyOp(big, end)" << endl;  
}

void MFHelmholtzOp::interpolateCF(const LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>* a_phiCoar, const bool a_homogeneousCF){
  CH_TIME("MFHelmholtzOp::interpolateCF");
  //if(m_debug) pout() << "MFHelmholtzOp::interpolateCF(begin)" << endl;
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> phiCoar;
      
    if(a_homogeneousCF){
      MultifluidAlias::aliasMF(phi, op.first, (LevelData<MFCellFAB>&) a_phi);

      op.second->homogeneousCFInterp(phi);
    }
    else{
      if(a_phiCoar == nullptr) MayDay::Error("MFHelmholtzOp::interpolateCF -- calling inhomogeneousCFInterp with nullptr coarse is an error.");
      
      MultifluidAlias::aliasMF(phi,     op.first,  (LevelData<MFCellFAB>&) a_phi);
      MultifluidAlias::aliasMF(phiCoar, op.first, *a_phiCoar);

      op.second->inhomogeneousCFInterp(phi, phiCoar);
    }
  }

  //if(m_debug) pout() << "MFHelmholtzOp::interpolateCF(end)" << endl;
}

void MFHelmholtzOp::residual(LevelData<MFCellFAB>& a_residual, const LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_rhs, const bool a_homogeneousPhysBC){
  CH_TIME("MFHelmholtzOp::residual");
  //if(m_debug) pout() << "MFHelmholtzOp::residual(begin)" << endl;
  // Compute a_residual = rhs - L(phi)

  this->applyOp(a_residual, a_phi, a_homogeneousPhysBC);

  this->scale(a_residual, -1.0);
  this->incr(a_residual, a_rhs, 1.0);


  //if(m_debug) pout() << "MFHelmholtzOp::residual(end)" << endl;
}

void MFHelmholtzOp::axby(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_x, const LevelData<MFCellFAB>& a_y, const Real a, const Real b) {
  CH_TIME("MFHelmholtzOp::axby");
  //if(m_debug) pout() << "MFHelmholtzOp::axby(begin)" << endl;

  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> lhs;
    LevelData<EBCellFAB> x;
    LevelData<EBCellFAB> y;

    MultifluidAlias::aliasMF(lhs, op.first, a_lhs);
    MultifluidAlias::aliasMF(x,   op.first, a_x);
    MultifluidAlias::aliasMF(y,   op.first, a_y);

    op.second->axby(lhs, x, y, a, b);
  }
  //if(m_debug) pout() << "MFHelmholtzOp::axby(end)" << endl;
}

void MFHelmholtzOp::updateJumpBC(const LevelData<MFCellFAB>& a_phi, const bool a_homogeneousPhysBC){
  CH_TIME("MFHelmholtzOp::updateJumpBC");
  //if(m_debug) pout() << "MFHelmholtzOp::updateJumpBC - not implemented" << endl;
}

void MFHelmholtzOp::relax(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_residual, int a_iterations) {
  CH_TIME("MFHelmholtzOp::relax");
  //if(m_debug) pout() << "MFHelmholtzOp::relax(begin)" << endl;
  
  // Operators are uncoupled for now.
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> correction;
    LevelData<EBCellFAB> residual;

    MultifluidAlias::aliasMF(correction, op.first, a_correction);
    MultifluidAlias::aliasMF(residual,   op.first, a_residual);

    op.second->relax(correction, residual, a_iterations);
  }
  //if(m_debug) pout() << "MFHelmholtzOp::relax(end)" << endl;
}

void MFHelmholtzOp::restrictResidual(LevelData<MFCellFAB>& a_resCoar, LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_rhs) {
  CH_TIME("MFHelmholtzOp::restrictResidual");
  //if(m_debug) pout() << "MFHelmholtzOp::restrictResidual(begin)" << endl;
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> resCoar;
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> rhs;

    MultifluidAlias::aliasMF(resCoar, op.first, a_resCoar);
    MultifluidAlias::aliasMF(phi,     op.first, a_phi);
    MultifluidAlias::aliasMF(rhs,     op.first, a_rhs);

    op.second->restrictResidual(resCoar, phi, rhs);
  }
  //if(m_debug) pout() << "MFHelmholtzOp::restrictResidual(end)" << endl;
}

void MFHelmholtzOp::prolongIncrement(LevelData<MFCellFAB>& a_phi, const LevelData<MFCellFAB>& a_correctCoarse) {
  CH_TIME("MFHelmholtzOp::prolongIncrement");
  //if(m_debug) pout() << "MFHelmholtzOp::prolongIncrement(begin)" << endl;
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> correctCoarse;

    MultifluidAlias::aliasMF(phi,           op.first, a_phi);
    MultifluidAlias::aliasMF(correctCoarse, op.first, a_correctCoarse);

    op.second->prolongIncrement(phi, correctCoarse);
  }
  //if(m_debug) pout() << "MFHelmholtzOp::prolongIncrement(end)" << endl;
}

void MFHelmholtzOp::AMRUpdateResidual(LevelData<MFCellFAB>&       a_residual,
				      const LevelData<MFCellFAB>& a_correction,
				      const LevelData<MFCellFAB>& a_coarseCorrection){
  CH_TIME("MFHelmholtzOp::AMRUpdateResidual");
  //if(m_debug) pout() << "MFHelmholtzOp::AMRUpdateResidual(begin)" << endl;

  // Need to update BC first!
  //this->updateJumpBC(a_correction, true);

  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> residual;
    LevelData<EBCellFAB> correction;
    LevelData<EBCellFAB> coarseCorrection;

    MultifluidAlias::aliasMF(residual,         op.first, a_residual);
    MultifluidAlias::aliasMF(correction,       op.first, a_correction);
    MultifluidAlias::aliasMF(coarseCorrection, op.first, a_coarseCorrection);

    op.second->AMRUpdateResidual(residual, correction, coarseCorrection);
  }
  //if(m_debug) pout() << "MFHelmholtzOp::AMRUpdateResidual(end)" << endl;
}

void MFHelmholtzOp::AMRRestrict(LevelData<MFCellFAB>&       a_residualCoarse,
				const LevelData<MFCellFAB>& a_residual,
				const LevelData<MFCellFAB>& a_correction,
				const LevelData<MFCellFAB>& a_coarseCorrection,
				bool                        a_skip_res) {
  CH_TIME("MFHelmholtzOp::AMRRestrict");
  //if(m_debug) pout() << "MFHelmholtzOp::AMRRestrict(begin)" << endl;
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> residualCoarse;
    LevelData<EBCellFAB> residual;
    LevelData<EBCellFAB> correction;
    LevelData<EBCellFAB> coarseCorrection;
    
    MultifluidAlias::aliasMF(residualCoarse,   op.first, a_residualCoarse);
    MultifluidAlias::aliasMF(residual,         op.first, a_residual);
    MultifluidAlias::aliasMF(correction,       op.first, a_correction);
    MultifluidAlias::aliasMF(coarseCorrection, op.first, a_coarseCorrection);

    op.second->AMRRestrict(residualCoarse, residual, correction, coarseCorrection, a_skip_res);
  }
  //if(m_debug) pout() << "MFHelmholtzOp::AMRRestrict(end)" << endl;
}

void MFHelmholtzOp::AMRProlong(LevelData<MFCellFAB>& a_correction, const LevelData<MFCellFAB>& a_coarseCorrection) {
  CH_TIME("MFHelmholtzOp::AMRProlong");
  //if(m_debug) pout() << "MFHelmholtzOp::AMRProlong(begin)" << endl;
  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> correction;
    LevelData<EBCellFAB> coarseCorrection;
    
    MultifluidAlias::aliasMF(correction,       op.first, a_correction);
    MultifluidAlias::aliasMF(coarseCorrection, op.first, a_coarseCorrection);

    op.second->AMRProlong(correction, coarseCorrection);
  }
  //if(m_debug) pout() << "MFHelmholtzOp::AMRProlong(end)" << endl;
}

void MFHelmholtzOp::AMRResidual(LevelData<MFCellFAB>&              a_residual,
				const LevelData<MFCellFAB>&        a_phiFine,
				const LevelData<MFCellFAB>&        a_phi,
				const LevelData<MFCellFAB>&        a_phiCoar,
				const LevelData<MFCellFAB>&        a_rhs,
				bool                               a_homogeneousPhysBC,
				AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp) {
  CH_TIME("MFHelmholtzOp::AMRResidual");
  //if(m_debug) pout() << "MFHelmholtzOp::AMRResidual(begin)" << endl;
  // Make residual = a_rhs - L(phi)
  this->AMROperator(a_residual, a_phiFine, a_phi, a_phiCoar, a_homogeneousPhysBC, a_finerOp);

  this->scale(a_residual, -1.0);
  this->incr(a_residual, a_rhs, 1.0);

  //if(m_debug) pout() << "MFHelmholtzOp::AMRResidual(end)" << endl;
}

void MFHelmholtzOp::AMRResidualNF(LevelData<MFCellFAB>&       a_residual,
				  const LevelData<MFCellFAB>& a_phi,
				  const LevelData<MFCellFAB>& a_phiCoar,
				  const LevelData<MFCellFAB>& a_rhs,
				  bool                        a_homogeneousPhysBC) {
  CH_TIME("MFHelmholtzOp::AMRResidualNF");
  //if(m_debug) pout() << "MFHelmholtzOp::AMRResidualNF(begin)" << endl;
  // Make residual = a_rhs - L(phi)  
  this->AMROperatorNF(a_residual, a_phi, a_phiCoar, a_homogeneousPhysBC);

  this->scale(a_residual, -1.0);
  this->incr(a_residual, a_rhs, 1.0);

  //if(m_debug) pout() << "MFHelmholtzOp::AMRResidualNF(end)" << endl;
}

void MFHelmholtzOp::AMRResidualNC(LevelData<MFCellFAB>&              a_residual,
				  const LevelData<MFCellFAB>&        a_phiFine,
				  const LevelData<MFCellFAB>&        a_phi,
				  const LevelData<MFCellFAB>&        a_rhs,
				  bool                               a_homogeneousPhysBC,
				  AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp) {
  CH_TIME("MFHelmholtzOp::AMRResidualNC");
  //if(m_debug) pout() << "MFHelmholtzOp::AMResidualNC(begin)" << endl;
  this->AMROperatorNC(a_residual, a_phiFine, a_phi, a_homogeneousPhysBC, a_finerOp);

  this->scale(a_residual, -1.0);
  this->incr(a_residual, a_rhs, 1.0);

  //if(m_debug) pout() << "MFHelmholtzOp::AMRResidualNC(end)" << endl;
}

void MFHelmholtzOp::AMROperatorNF(LevelData<MFCellFAB>&       a_Lphi,
				  const LevelData<MFCellFAB>& a_phi,
				  const LevelData<MFCellFAB>& a_phiCoar,
				  bool                        a_homogeneousPhysBC) {
  CH_TIME("MFHelmholtzOp::AMROperatorNF");
  //if(m_debug) pout() << "MFHelmholtzOp::AMROperatorNF(begin)" << endl;
  this->interpolateCF(a_phi, &a_phiCoar, false);
  //this->updateJumpBC(a_phi, a_homogeneousPhysBC);

  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> Lphi;
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> phiCoar;

    MultifluidAlias::aliasMF(Lphi,    op.first, a_Lphi);
    MultifluidAlias::aliasMF(phi,     op.first, a_phi);
    MultifluidAlias::aliasMF(phiCoar, op.first, a_phiCoar);

    op.second->turnOffBCs(); // Don't need to interpolate ghost cells again. 
    op.second->AMROperatorNF(Lphi, phi, phiCoar, a_homogeneousPhysBC);
    op.second->turnOnBCs(); 
  }
  //if(m_debug) pout() << "MFHelmholtzOp::AMROperatorNF(end)" << endl;
}

void MFHelmholtzOp::AMROperatorNC(LevelData<MFCellFAB>&              a_Lphi,
				  const LevelData<MFCellFAB>&        a_phiFine,
				  const LevelData<MFCellFAB>&        a_phi,
				  bool                               a_homogeneousPhysBC,
				  AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp) {
  CH_TIME("MFHelmholtzOp::AMROperatorNC");
  //if(m_debug) pout() << "MFHelmholtzOp::AMROperatorNC(begin)" << endl;

  // Must update the jump BC first. 
  //this->updateJumpBC(a_phi, a_homogeneousPhysBC);

  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> Lphi;
    LevelData<EBCellFAB> phiFine;
    LevelData<EBCellFAB> phi;

    MultifluidAlias::aliasMF(Lphi,    op.first, a_Lphi   );
    MultifluidAlias::aliasMF(phiFine, op.first, a_phiFine);
    MultifluidAlias::aliasMF(phi,     op.first, a_phi    );

    MFHelmholtzOp* finerOp = (MFHelmholtzOp*) (a_finerOp);

    op.second->AMROperatorNC(Lphi, phiFine, phi, a_homogeneousPhysBC, (finerOp->m_helmOps).at(op.first));
  }
  
  //if(m_debug) pout() << "MFHelmholtzOp::AMROperatorNC(end)" << endl;
}

void MFHelmholtzOp::AMROperator(LevelData<MFCellFAB>&              a_Lphi,
				const LevelData<MFCellFAB>&        a_phiFine,
				const LevelData<MFCellFAB>&        a_phi,
				const LevelData<MFCellFAB>&        a_phiCoar,
				const bool                         a_homogeneousPhysBC,
				AMRLevelOp<LevelData<MFCellFAB> >* a_finerOp) {
  CH_TIME("MFHelmholtzOp::AMROperator");
  //if(m_debug) pout() << "MFHelmholtzOp::AMROperator(begin)" << endl;
  // Update ghost cells and jump conditions first.
  this->interpolateCF(a_phi, &a_phiCoar, false);
  //this->updateJumpBC(a_phi, a_homogeneousPhysBC);

  for (auto& op : m_helmOps){
    LevelData<EBCellFAB> Lphi;
    LevelData<EBCellFAB> phiFine;
    LevelData<EBCellFAB> phi;
    LevelData<EBCellFAB> phiCoar;

    MultifluidAlias::aliasMF(Lphi,    op.first, a_Lphi);
    MultifluidAlias::aliasMF(phiFine, op.first, a_phiFine);
    MultifluidAlias::aliasMF(phi,     op.first, a_phi);
    MultifluidAlias::aliasMF(phiCoar, op.first, a_phiCoar);

    MFHelmholtzOp* finerOp = (MFHelmholtzOp*) (a_finerOp);

    op.second->turnOffBCs(); // Don't need to interpolate ghost cells again.
    op.second->AMROperator(Lphi, phiFine, phi, phiCoar, a_homogeneousPhysBC, (finerOp->m_helmOps).at(op.first));
    op.second->turnOnBCs();
  }
  //if(m_debug) pout() << "MFHelmholtzOp::AMROperator(end)" << endl;
}

#include <CD_NamespaceFooter.H>
