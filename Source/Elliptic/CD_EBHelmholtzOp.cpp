/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzOp.cpp
  @brief  Implementation of CD_EBHelmholtzOp.H
  @author Robert Marskar
  @todo   Replace EBAMRPoissonOp::staticMaxNorm and don't use EBAMRPoissonOp dependencies
  @todo   When we redefine the EBBC, the getFluxStencil should return an empty stencil rather than a null pointer. Fix this in defineStencils() together with the other crap. 
*/

// Chombo includes
#include <EBLevelDataOps.H>
#include <EBCellFactory.H>
#include <EBAMRPoissonOp.H>
#include <EBConductivityOpF_F.H>

// Our includes
#include <CD_EBHelmholtzOp.H>
#include <CD_EBHelmholtzOpF_F.H>
#include <CD_NamespaceHeader.H>

EBHelmholtzOp::EBHelmholtzOp(const EBLevelGrid&                                 a_eblgFine,
			     const EBLevelGrid&                                 a_eblg,
			     const EBLevelGrid&                                 a_eblgCoFi,
			     const EBLevelGrid&                                 a_eblgCoar,
			     const EBLevelGrid&                                 a_eblgCoarMG,
			     const RefCountedPtr<EBMultigridInterpolator>&      a_interpolator,
			     const RefCountedPtr<EBFluxRegister>&               a_fluxReg,
			     const RefCountedPtr<EbCoarAve>&                    a_coarAve,			       
			     const RefCountedPtr<EBHelmholtzDomainBc>&          a_domainBc,
			     const RefCountedPtr<EBHelmholtzEbBc>&              a_ebBc,
			     const RealVect&                                    a_probLo,
			     const Real&                                        a_dx,
			     const int&                                         a_refToFine,
			     const int&                                         a_refToCoar,
			     const bool&                                        a_hasFine,
			     const bool&                                        a_hasCoar,
			     const bool&                                        a_hasMGObjects,
			     const Real&                                        a_alpha,
			     const Real&                                        a_beta,
			     const RefCountedPtr<LevelData<EBCellFAB> >&        a_Acoef,
			     const RefCountedPtr<LevelData<EBFluxFAB> >&        a_Bcoef,
			     const RefCountedPtr<LevelData<BaseIVFAB<Real> > >& a_BcoefIrreg,
			     const IntVect&                                     a_ghostPhi,
			     const IntVect&                                     a_ghostRhs,
			     const RelaxationMethod&                            a_relaxationMethod) :
  LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB>(false), // Time-independent
  m_eblgFine(),
  m_eblg(a_eblg),
  m_eblgCoFi(),
  m_eblgCoar(),
  m_eblgCoarMG(),
  m_interpolator(a_interpolator),
  m_fluxReg(a_fluxReg),
  m_coarAve(a_coarAve),
  m_domainBc(a_domainBc),
  m_ebBc(a_ebBc),
  m_probLo(a_probLo),
  m_dx(a_dx),
  m_refToFine(a_hasFine ? a_refToFine : 1),
  m_refToCoar(a_hasCoar ? a_refToCoar : 1),
  m_hasFine(a_hasFine),
  m_hasCoar(a_hasCoar),
  m_hasMGObjects(a_hasMGObjects),
  m_alpha(a_alpha),
  m_beta(a_beta),
  m_Acoef(a_Acoef),
  m_Bcoef(a_Bcoef),
  m_BcoefIrreg(a_BcoefIrreg),
  m_ghostPhi(a_ghostPhi),
  m_ghostRhs(a_ghostRhs),
  m_relaxationMethod(a_relaxationMethod) {

  // Default settings. Always solve for comp = 0. If you want something different, copy your
  // input two different data holders before you use AMRMultiGrid. 
  m_nComp      = 1;
  m_comp       = 0;
  m_turnOffBCs = false;
  m_vecDx      = m_dx*RealVect::Unit;
  m_interval   = Interval(m_comp, m_comp);

  if(m_hasFine){
    m_eblgFine = a_eblgFine;
  }

  if(m_hasCoar){
    m_eblgCoFi = a_eblgCoFi;
    m_eblgCoar = a_eblgCoar;

    m_ebInterp.define(m_eblg.getDBL(),   m_eblgCoar.getDBL(), m_eblg.getEBISL(), m_eblgCoar.getEBISL(), m_eblgCoar.getDomain(),
    		      m_refToCoar, m_nComp, m_eblg.getEBIS(), m_ghostPhi);

    m_ebAverage.define(m_eblg.getDBL(), m_eblgCoFi.getDBL(), m_eblg.getEBISL(), m_eblgCoFi.getEBISL(), m_eblgCoFi.getDomain(),
    		       m_refToCoar, m_nComp, m_eblg.getEBIS(), m_ghostRhs);
  }

  if(m_hasMGObjects){
    constexpr int mgRef = 2;
    
    m_eblgCoarMG = a_eblgCoarMG;

    m_ebInterpMG.define(m_eblg.getDBL(), m_eblgCoarMG.getDBL(), m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(), m_eblgCoarMG.getDomain(),
			mgRef, m_nComp, m_eblg.getEBIS(), m_ghostPhi);

    m_ebAverageMG.define(m_eblg.getDBL(), m_eblgCoarMG.getDBL(), m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(), m_eblgCoarMG.getDomain(),
			 mgRef, m_nComp, m_eblg.getEBIS(), m_ghostRhs);
  }

  // Define data holders and stencils
  this->defineStencils();
}

EBHelmholtzOp::~EBHelmholtzOp(){

}

LevelData<EBFluxFAB>& EBHelmholtzOp::getFlux(){
  return m_flux;
}

void EBHelmholtzOp::defineStencils(){
  // Basic defines. 
  EBCellFactory cellFact(m_eblg.getEBISL());
  EBFluxFactory fluxFact(m_eblg.getEBISL());
  
  m_relCoef.define(m_eblg.getDBL(), m_nComp, IntVect::Zero, cellFact); // Don't need ghost cells for this one. 
  m_flux.define(   m_eblg.getDBL(), m_nComp, IntVect::Unit, fluxFact); // One ghost cell needed for face center to face centroid interpolation. 

  m_opEBStencil.define(    m_eblg.getDBL());
  m_vofIterIrreg.define(   m_eblg.getDBL());
  m_vofIterMulti.define(   m_eblg.getDBL());
  m_alphaDiagWeight.define(m_eblg.getDBL());
  m_betaDiagWeight.define( m_eblg.getDBL());

  for (int dir = 0; dir < SpaceDim; dir++){
    m_vofIterDomLo [dir].define(m_eblg.getDBL());
    m_vofIterDomHi [dir].define(m_eblg.getDBL());
    m_interpolant  [dir].define(m_eblg.getDBL());
    m_interpStencil[dir].define(m_eblg.getDBL());
  }

  EBArith::getMultiColors(m_colors);

  // First strip of cells on the inside of the computational domain. I.e. the "domain wall" cells where we need BCs. 
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      Box domainBox = m_eblg.getDomain().domainBox();
      Box sidebox   = adjCellBox(domainBox, dir, sit(), 1);
      sidebox.shift(dir, -sign(sit()));
      m_sideBox.emplace(std::make_pair(dir, sit()), sidebox);
    }
  }

  // Fuck I hate the BC classes in Chombo.
  // TODO: Remember to replace this stuff...
  Real fakeBeta = 1.0;
  m_domainBc->setCoef(m_eblg, fakeBeta, m_Bcoef);
  m_ebBc->setCoef(    m_eblg, fakeBeta, m_BcoefIrreg);
  m_ebBc->define((*m_eblg.getCFIVS()), 1./m_dx);
  const LayoutData<BaseIVFAB<VoFStencil> >* const ebFluxStencil = m_ebBc->getFluxStencil(m_comp);

  // Define everything
  for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit){
    const Box cellBox      = m_eblg.getDBL()[dit()];
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();

    const IntVectSet irregIVS = ebisbox.getIrregIVS(cellBox);
    const IntVectSet multiIVS = ebisbox.getMultiCells(cellBox);

    m_vofIterIrreg[dit()].define(irregIVS, ebgraph);
    m_vofIterMulti[dit()].define(multiIVS, ebgraph);

    m_alphaDiagWeight[dit()].define(irregIVS, ebgraph, m_nComp);
    m_betaDiagWeight[dit()].define( irregIVS, ebgraph, m_nComp);

    for (int dir = 0; dir < SpaceDim; dir++){
      const IntVectSet loIrreg = irregIVS & m_sideBox.at(std::make_pair(dir, Side::Lo));
      const IntVectSet hiIrreg = irregIVS & m_sideBox.at(std::make_pair(dir, Side::Hi));

      m_vofIterDomLo[dir][dit()].define(loIrreg, ebgraph);
      m_vofIterDomHi[dir][dit()].define(hiIrreg, ebgraph);

      // Define the face interpolation stencils
      m_interpStencil[dir][dit()].define(irregIVS, ebgraph, dir, m_nComp);
      m_interpolant  [dir][dit()].define(irregIVS, ebgraph, dir, m_nComp);
      for (FaceIterator faceIt(irregIVS, ebgraph, dir, FaceStop::SurroundingNoBoundary); faceIt.ok(); ++faceIt){
	m_interpStencil[dir][dit()](faceIt(), m_comp) = EBArith::getInterpStencil(faceIt(), IntVectSet(), ebisbox, m_eblg.getDomain());
      }
    }


    // Build stencils and weights for irregular cells.
    VoFIterator& vofit = m_vofIterIrreg[dit()];
    BaseIVFAB<VoFStencil> opStencil(irregIVS, ebgraph, m_nComp);
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const IntVect& iv   = vof.gridIndex();

      VoFStencil& curStencil = opStencil(vof, m_comp);

      // Get stencil for this cell. 
      curStencil = this->getDivFStencil(vof, dit());
      if(ebFluxStencil != nullptr) curStencil += (*ebFluxStencil)[dit()](vof, m_comp);

      // Adjust the weight with domain boundary faces. 
      Real betaWeight = EBArith::getDiagWeight(curStencil, vof);
      for (int dir = 0; dir < SpaceDim; dir++){
	for (SideIterator sit; sit.ok(); ++sit){
	  const Box sidebox = m_sideBox.at(std::make_pair(dir, sit()));

	  if(sidebox.contains(iv)){
	    Real weightedAreaFrac = 0.0;
	    Vector<FaceIndex> faces = ebisbox.getFaces(vof, dir, sit());
	    for (auto& f : faces.stdVector()){
	      weightedAreaFrac += ebisbox.areaFrac(f)*(*m_Bcoef)[dit()][dir](f, m_comp)/(m_dx*m_dx);
	    }
	    betaWeight += -weightedAreaFrac;
	  }
	}
      }
      
      m_betaDiagWeight[dit()](vof, m_comp) = betaWeight;
    }

    m_opEBStencil[dit()] = RefCountedPtr<EBStencil>(new EBStencil(m_vofIterIrreg[dit()].getVector(), opStencil, cellBox, ebisbox, m_ghostPhi, m_ghostRhs, m_comp, true));
  }
    
  // Compute the alpha-weight and relaxation coefficient. 
  this->computeAlphaWeight();
  this->computeRelaxationCoefficient();
}

unsigned int EBHelmholtzOp::orderOfAccuracy(void) const {
  return 99;
}

void EBHelmholtzOp::enforceCFConsistency(LevelData<EBCellFAB>& a_coarCorr, const LevelData<EBCellFAB>& a_fineCorr){
  m_coarAve->average(a_coarCorr, a_fineCorr, a_coarCorr.interval());
}

void EBHelmholtzOp::setAlphaAndBeta(const Real& a_alpha, const Real& a_beta) {
  m_alpha = a_alpha;
  m_beta  = a_beta;

  // When we change alpha and beta we need to recompute relaxation coefficients...
  this->computeAlphaWeight(); 
  this->computeRelaxationCoefficient();
}

void EBHelmholtzOp::residual(LevelData<EBCellFAB>& a_residual, const LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_rhs, bool a_homogeneousPhysBC) {
  this->applyOp(a_residual, a_phi, nullptr, a_homogeneousPhysBC, true); // 
  this->axby(a_residual, a_residual, a_rhs, -1.0, 1.0);                 // residual = rhs - L(phi). 
}

void EBHelmholtzOp::preCond(LevelData<EBCellFAB>& a_corr, const LevelData<EBCellFAB>& a_residual) {
  EBLevelDataOps::assign(a_corr, a_residual);
  EBLevelDataOps::scale(a_corr,  m_relCoef);

  this->relax(a_corr, a_residual, 40);
}

void EBHelmholtzOp::create(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) {
  EBCellFactory fact(m_eblg.getEBISL());
  a_lhs.define(m_eblg.getDBL(), a_rhs.nComp(), a_rhs.ghostVect(), fact);
}

void EBHelmholtzOp::assign(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) {
  a_rhs.copyTo(a_lhs);
}

Real EBHelmholtzOp::dotProduct(const LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs) {
  ProblemDomain domain;
  Real volume;

  return EBLevelDataOps::kappaDotProduct(volume, a_lhs, a_rhs, EBLEVELDATAOPS_ALLVOFS, domain);
}

void EBHelmholtzOp::incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale) {
  EBLevelDataOps::incr(a_lhs, a_rhs, a_scale);
}

void EBHelmholtzOp::axby(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_x, const LevelData<EBCellFAB>& a_y, const Real a_a, const Real a_b) {
  EBLevelDataOps::axby(a_lhs, a_x, a_y, a_a, a_b);
}

void EBHelmholtzOp::scale(LevelData<EBCellFAB>& a_lhs, const Real& a_scale) {
  EBLevelDataOps::scale(a_lhs, a_scale);
}

Real EBHelmholtzOp::norm(const LevelData<EBCellFAB>& a_rhs, const int a_order) {
  Real maxNorm = EBAMRPoissonOp::staticMaxNorm(a_rhs, m_eblg);

#ifdef CH_MPI
  Real tmp = 1;
  MPI_Allreduce(&maxNorm, &tmp, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
  maxNorm = tmp;
#endif

  return maxNorm;
}

void EBHelmholtzOp::setToZero(LevelData<EBCellFAB>& a_lhs) {
  EBLevelDataOps::setToZero(a_lhs);
}

void EBHelmholtzOp::createCoarser(LevelData<EBCellFAB>& a_coarse, const LevelData<EBCellFAB>& a_fine, bool a_ghosted) {
  EBCellFactory factCoar(m_eblgCoarMG.getEBISL());
  a_coarse.define(m_eblgCoarMG.getDBL(), a_fine.nComp(), a_fine.ghostVect(), factCoar);
}

void EBHelmholtzOp::createCoarsened(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int& a_refRat) {
  CH_assert(m_hasCoar);

  EBCellFactory factCoFi(m_eblgCoFi.getEBISL());
  a_lhs.define(m_eblgCoFi.getDBL(), a_rhs.nComp(), a_rhs.ghostVect(), factCoFi);
}

void EBHelmholtzOp::restrictResidual(LevelData<EBCellFAB>& a_resCoar, LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_rhs) {

  // Compute the residual on this level first. Make a temporary for that.
  LevelData<EBCellFAB> res;
  this->create(res, a_phi);
  this->setToZero(res);
  this->residual(res, a_phi, a_rhs, true);

  m_ebAverageMG.average(a_resCoar, res, m_interval);
}

void EBHelmholtzOp::prolongIncrement(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_correctCoarse) {
  m_ebInterpMG.pwcInterp(a_phi, a_correctCoarse, m_interval);
}

int EBHelmholtzOp::refToCoarser() {
  return m_refToCoar;
}

void EBHelmholtzOp::AMROperator(LevelData<EBCellFAB>&              a_Lphi,
				const LevelData<EBCellFAB>&        a_phiFine,
				const LevelData<EBCellFAB>&        a_phi,
				const LevelData<EBCellFAB>&        a_phiCoar,
				const bool                         a_homogeneousPhysBC,
				AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp){
  this->applyOp(a_Lphi, a_phi, &a_phiCoar, a_homogeneousPhysBC, false);

  if(m_hasFine){
    this->reflux(a_Lphi, a_phiFine, a_phi, *a_finerOp);
  }
}

void EBHelmholtzOp::AMROperatorNF(LevelData<EBCellFAB>&       a_Lphi,
				  const LevelData<EBCellFAB>& a_phi,
				  const LevelData<EBCellFAB>& a_phiCoar,
				  bool                        a_homogeneousPhysBC) {
  this->applyOp(a_Lphi, a_phi, &a_phiCoar, a_homogeneousPhysBC, false);
}

void EBHelmholtzOp::AMROperatorNC(LevelData<EBCellFAB>&              a_Lphi,
				  const LevelData<EBCellFAB>&        a_phiFine,
				  const LevelData<EBCellFAB>&        a_phi,
				  bool                               a_homogeneousPhysBC,
				  AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp) {
  LevelData<EBCellFAB> phiCoar; // Should be safe on the bottom AMR level because only multigrid levels exist below. 
  this->AMROperator(a_Lphi, a_phiFine, a_phi, phiCoar, a_homogeneousPhysBC, a_finerOp); 
}

void EBHelmholtzOp::AMRResidual(LevelData<EBCellFAB>&              a_residual,
				const LevelData<EBCellFAB>&        a_phiFine,
				const LevelData<EBCellFAB>&        a_phi,
				const LevelData<EBCellFAB>&        a_phiCoar,
				const LevelData<EBCellFAB>&        a_rhs,
				bool                               a_homogeneousPhysBC,
				AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp) {
  this->AMROperator(a_residual, a_phiFine, a_phi, a_phiCoar, a_homogeneousPhysBC, a_finerOp); // a_residual =  L(phi) 
  this->axby(a_residual, a_rhs, a_residual, 1., -1.);                                         // a_residual = rhs - -L(phi)
}

void EBHelmholtzOp::AMRResidualNF(LevelData<EBCellFAB>&              a_residual,
				  const LevelData<EBCellFAB>&        a_phi,
				  const LevelData<EBCellFAB>&        a_phiCoar,
				  const LevelData<EBCellFAB>&        a_rhs,
				  bool                               a_homogeneousPhysBC) {

  // Simple, because we don't need to reflux. 
  this->AMROperatorNF(a_residual, a_phi, a_phiCoar, a_homogeneousPhysBC); // a_residual = L(phi)
  this->axby(a_residual, a_rhs, a_residual, 1., -1.);                     // a_residual = rhs - L(phi)
}

void EBHelmholtzOp::AMRResidualNC(LevelData<EBCellFAB>&              a_residual,
				  const LevelData<EBCellFAB>&        a_phiFine,
				  const LevelData<EBCellFAB>&        a_phi,
				  const LevelData<EBCellFAB>&        a_rhs,
				  bool                               a_homogeneousPhysBC,
				  AMRLevelOp<LevelData<EBCellFAB> >* a_finerOp) {
  this->AMROperatorNC(a_residual, a_phiFine, a_phi, a_homogeneousPhysBC, a_finerOp); // a_residual = L(phi)
  this->axby(a_residual, a_rhs, a_residual, 1., -1.);                                // a_residual = rhs - L(phi)
}

void EBHelmholtzOp::AMRRestrict(LevelData<EBCellFAB>&       a_residualCoarse,
				const LevelData<EBCellFAB>& a_residual,
				const LevelData<EBCellFAB>& a_correction,
				const LevelData<EBCellFAB>& a_coarseCorrection,
				bool                        a_skip_res) {

  LevelData<EBCellFAB> resThisLevel;
  this->create(resThisLevel, a_residual);
  this->setToZero(resThisLevel);

  const bool homogeneousPhysBC = true;
  const bool homogeneousCFBC   = false;

  // We should average a_residual - L(correction, coarCorrection).
  this->applyOp(resThisLevel, a_correction, &a_coarseCorrection, homogeneousPhysBC, homogeneousCFBC);
  this->incr(resThisLevel, a_residual, -1.0);
  this->scale(resThisLevel, -1.0);

  m_ebAverage.average(a_residualCoarse, resThisLevel, m_interval);
}

void EBHelmholtzOp::AMRProlong(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_coarseCorrection) {
  m_ebInterp.pwcInterp(a_correction, a_coarseCorrection, m_interval);
}

void EBHelmholtzOp::AMRUpdateResidual(LevelData<EBCellFAB>&       a_residual,
				      const LevelData<EBCellFAB>& a_correction,
				      const LevelData<EBCellFAB>& a_coarseCorrection) {

  const bool homogeneousPhysBC = true;
  const bool homogeneousCFBC   = false;
  
  LevelData<EBCellFAB> lcorr;
  this->create(lcorr, a_correction);
  this->applyOp(lcorr, a_correction, &a_coarseCorrection, homogeneousPhysBC, homogeneousCFBC); // lcorr = L(phi, phiCoar)
  this->incr(a_residual, lcorr, -1.0);                                                         // a_residual = a_residual - L(phi)
}

void EBHelmholtzOp::applyOp(LevelData<EBCellFAB>& a_Lphi, const LevelData<EBCellFAB>& a_phi, bool a_homogeneousPhysBC) {
  this->applyOp(a_Lphi, a_phi, nullptr, a_homogeneousPhysBC, true);
}

void EBHelmholtzOp::applyOp(LevelData<EBCellFAB>&             a_Lphi,
			    const LevelData<EBCellFAB>&       a_phi,
			    const LevelData<EBCellFAB>* const a_phiCoar,
			    const bool                        a_homogeneousPhysBC,
			    const bool                        a_homogeneousCFBC){

  LevelData<EBCellFAB>& phi = (LevelData<EBCellFAB>& ) a_phi;

  if(m_hasCoar && !m_turnOffBCs){
    this->interpolateCF(phi, *a_phiCoar, a_homogeneousCFBC);
  }
  
  phi.exchange();
  
  const DisjointBoxLayout& dbl = a_Lphi.disjointBoxLayout();
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box cellBox = dbl[dit()];
    
    // Now do the terms with beta*phi
    this->applyOpRegular(  a_Lphi[dit()], phi[dit()], cellBox, dit(), a_homogeneousPhysBC);
    this->applyOpIrregular(a_Lphi[dit()], phi[dit()], cellBox, dit(), a_homogeneousPhysBC); // Overwrites result in irregular cells. 
  }
}

void EBHelmholtzOp::applyOpRegular(EBCellFAB& a_Lphi, EBCellFAB& a_phi, const Box& a_cellBox, const DataIndex& a_dit, const bool a_homogeneousPhysBC){
  
  // Alpha term
  a_Lphi.setVal(0.0);
  a_Lphi += (*m_Acoef)[a_dit];
  a_Lphi *= m_alpha;

  // Add domain flux to a_Lphi. 
  this->applyDomainFlux(a_Lphi, a_phi, a_cellBox, a_dit, a_homogeneousPhysBC);

  // It's an older code, sir, but it checks out.
  BaseFab<Real>& Lphi      = a_Lphi.getSingleValuedFAB();
  const BaseFab<Real>& phi = a_phi.getSingleValuedFAB(); 
  BaseFab<Real> dummy(Box(IntVect::Zero, IntVect::Zero), 1);
  BaseFab<Real>* bcoef[3];

  for (int dir = 0; dir < 3; dir++){
    if(dir >= SpaceDim){
      bcoef[dir] = &dummy;
    }
    else{
      bcoef[dir] = &(*m_Bcoef)[a_dit][dir].getSingleValuedFAB();
    }
  }

  // Fill the ghost cells outside the domain so that we don't get an extra flux from the domain side. 
  FORT_LAPLACIANINPLACE(CHF_FRA1(Lphi, m_comp),
			CHF_CONST_FRA1(phi, m_comp),
			CHF_CONST_FRA1((*bcoef[0]), m_comp),
			CHF_CONST_FRA1((*bcoef[1]), m_comp),
			CHF_CONST_FRA1((*bcoef[2]), m_comp), // Not used for 2D. 
			CHF_CONST_REAL(m_dx),
			CHF_CONST_REAL(m_beta),
			CHF_BOX(a_cellBox));
}

void EBHelmholtzOp::applyDomainFlux(EBCellFAB& a_Lphi, EBCellFAB& a_phi, const Box& a_cellBox, const DataIndex& a_dit, const bool a_homogeneousPhysBC){
  // TLDR: While we're developing this operator I'm resuing the old BC classes. This routine is weird because the domain bc needs a cell-centered flux
  //       box, so that's what all the shifting is about.
  //
  //       I am leaving both a_Lphi and a_phi as options here because we may want to fill a_Lphi directly so as to avoid the current division-by-B in
  //       the Fortran kernels. In this case one can zero the domain flux by calling suppressDomainFlux...

  BaseFab<Real>& phiFAB = a_phi.getSingleValuedFAB();
  
  for (int dir = 0; dir < SpaceDim; dir++){

    Box loBox;
    Box hiBox;
    int hasLo;
    int hasHi;
    EBArith::loHi(loBox, hasLo, hiBox, hasHi, m_eblg.getDomain(), a_cellBox, dir);
    
    if(hasLo == 1){
      int side = -1;
      Box lbox = loBox;
      lbox.shift(dir, -1);
      FArrayBox faceFlux(loBox, m_nComp);
      m_domainBc->setCoef(m_eblg, m_beta, m_Bcoef);
      m_domainBc->getFaceFlux(faceFlux, a_phi.getSingleValuedFAB(), m_probLo, m_vecDx, dir, Side::Lo, a_dit, 0.0, a_homogeneousPhysBC);

      BaseFab<Real>& bco = (*m_Bcoef)[a_dit][dir].getSingleValuedFAB();
      bco.shiftHalf(dir, 1);
      FORT_HELMHOLTZAPPLYDOMAINFLUX(CHF_FRA1(phiFAB, m_comp),
      				    CHF_CONST_FRA1(faceFlux, m_comp),
      				    CHF_CONST_FRA1(bco,m_comp),
      				    CHF_CONST_REAL(m_dx),
				    CHF_CONST_INT(side),
				    CHF_CONST_INT(dir),
				    CHF_BOX(lbox));
      bco.shiftHalf(dir, -1);
    }

    if(hasHi == 1){
      int side = 1;
      Box hbox = hiBox;
      hbox.shift(dir, 1);
      FArrayBox faceFlux(hiBox, m_nComp);
      m_domainBc->setCoef(m_eblg, m_beta, m_Bcoef);
      m_domainBc->getFaceFlux(faceFlux, a_phi.getSingleValuedFAB(), m_probLo, m_vecDx, dir, Side::Hi, a_dit, 0.0, a_homogeneousPhysBC);

      BaseFab<Real>& bco     = (*m_Bcoef)[a_dit][dir].getSingleValuedFAB();
      bco.shiftHalf(dir, -1);
      FORT_HELMHOLTZAPPLYDOMAINFLUX(CHF_FRA1(phiFAB, m_comp),
      				    CHF_CONST_FRA1(faceFlux, m_comp),
      				    CHF_CONST_FRA1(bco,m_comp),
      				    CHF_CONST_REAL(m_dx),
				    CHF_CONST_INT(side),
				    CHF_CONST_INT(dir),
				    CHF_BOX(hbox));
      bco.shiftHalf(dir, 1);
    }
  }

  //  this->suppressDomainFlux(a_phi, a_cellBox, a_dit);
}

void EBHelmholtzOp::suppressDomainFlux(EBCellFAB& a_phi, const Box& a_cellBox, const DataIndex& a_dit){
  const Box domainBox = m_eblg.getDomain().domainBox();

  BaseFab<Real>& regFab = a_phi.getSingleValuedFAB();
  
  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const Box ghostBox = adjCellBox(domainBox, dir, sit(), 1); // Strip of cells outside domain.
      const int s        = sign(flip(sit())); // Make +1 for Side::Lo and -1 otherwise.
      
      // Phi defined over a_cellBox - get the ghost in a_phi that lie outside the domain and set them to the value in the cell 
      // immediately on the inside. 
      Box grownBox = a_cellBox;
      grownBox.growDir(dir, sit(), 1);
      grownBox &= ghostBox;

      // Set the ghost cells to the same value as the next "inside" cell. This will suppress the gradient on the domain faces. 
      for (BoxIterator bit(grownBox); bit.ok(); ++bit){
	const IntVect iv = bit();
	const IntVect ivNeigh = iv + s*BASISV(dir);

	regFab(iv, m_comp) = regFab(ivNeigh, m_comp);
      }
    }
  }
}

void EBHelmholtzOp::applyOpIrregular(EBCellFAB& a_Lphi, const EBCellFAB& a_phi, const Box& a_cellBox, const DataIndex& a_dit, const bool a_homogeneousPhysBC){
  m_opEBStencil[a_dit]->apply(a_Lphi, a_phi, m_alphaDiagWeight[a_dit], m_alpha, m_beta, false); // This includes stencils for interpolation to faces. 

  if(!a_homogeneousPhysBC){
    const Real factor = m_beta/m_dx; // Beta not handled inside m_ebBc but we will fix this (later). 
    m_ebBc->applyEBFlux(a_Lphi, a_phi, m_vofIterIrreg[a_dit], (*m_eblg.getCFIVS()), a_dit, m_probLo, m_vecDx, factor, a_homogeneousPhysBC, 0.0);
  }

  // Do irregular faces on domain sides. m_domainBc should give the centroid-centered flux so we don't do interpolations here. 
  for (int dir = 0; dir < SpaceDim; dir++){
    Real flux = 0.0;

    // Lo side.
    VoFIterator& vofitLo = m_vofIterDomLo[dir][a_dit];    
    for (vofitLo.reset(); vofitLo.ok(); ++vofitLo){
      const VolIndex& vof = vofitLo();
      m_domainBc->getFaceFlux(flux, vof, m_comp, a_phi, m_probLo, m_vecDx, dir, Side::Lo, a_dit, 0.0, a_homogeneousPhysBC);

      a_Lphi(vof, m_comp) -= flux*m_beta/m_dx;
    }

    // Hi side. 
    VoFIterator& vofitHi = m_vofIterDomHi[dir][a_dit];
    for (vofitHi.reset(); vofitHi.ok(); ++vofitHi){
      const VolIndex& vof = vofitHi();

      m_domainBc->getFaceFlux(flux, vof, m_comp, a_phi, m_probLo, m_vecDx, dir, Side::Hi, a_dit, 0.0, a_homogeneousPhysBC);

      a_Lphi(vof, m_comp) += flux*m_beta/m_dx;
    }
  }
}

void EBHelmholtzOp::divideByIdentityCoef(LevelData<EBCellFAB>& a_rhs) {
  for (DataIterator dit(a_rhs.disjointBoxLayout()); dit.ok(); ++dit){
    a_rhs[dit()] /= (*m_Acoef)[dit()];
  }
}

void EBHelmholtzOp::applyOpNoBoundary(LevelData<EBCellFAB>& a_Lphi, const LevelData<EBCellFAB>& a_phi) {
  m_turnOffBCs = true;
  this->applyOp(a_Lphi, a_phi, true);
  m_turnOffBCs = false;
}

void EBHelmholtzOp::fillGrad(const LevelData<EBCellFAB>& a_phi){
  MayDay::Abort("EBHelmholtzOp::fillGrad - not implemented (yet)");
}

void EBHelmholtzOp::getFlux(EBFluxFAB&                  a_flux,
			    const LevelData<EBCellFAB>& a_data,
			    const Box&                  a_grid,
			    const DataIndex&            a_dit,
			    Real                        a_scale) {
  MayDay::Abort("EBHelmholtzOp::getFlux - not implemented (yet)");
}

void EBHelmholtzOp::homogeneousCFInterp(LevelData<EBCellFAB>& a_phi){
  if(m_hasCoar) m_interpolator->coarseFineInterpH(a_phi, m_interval);
}

void EBHelmholtzOp::inhomogeneousCFInterp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar){
  if(m_hasCoar) m_interpolator->coarseFineInterp(a_phiFine, a_phiCoar, m_interval);
}

void EBHelmholtzOp::interpolateCF(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar, const bool a_homogeneousCFBC){
  if(m_hasCoar){
    if(a_homogeneousCFBC){
      this->homogeneousCFInterp(a_phiFine);
    }
    else{
      this->inhomogeneousCFInterp(a_phiFine, a_phiCoar);
    }
  }
}

void EBHelmholtzOp::relax(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, int a_iterations){
  switch(m_relaxationMethod){
  case RelaxationMethod::Jacobi:
    this->relaxJacobi(a_correction, a_residual, a_iterations);
    break;
  case RelaxationMethod::GSRB:
    this->relaxGauSai(a_correction, a_residual, a_iterations);
    break;
  case RelaxationMethod::GSRBFast:
    this->relaxGSRBFast(a_correction, a_residual, a_iterations);
    break;
  default:
    MayDay::Abort("EBHelmholtzOp::relax - bogus relaxation method requested");
  };
}

void EBHelmholtzOp::relaxJacobi(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations){
  LevelData<EBCellFAB> Lcorr;
  this->create(Lcorr, a_residual);

  for (int iter = 0; iter < a_iterations; iter++){
    this->homogeneousCFInterp(a_correction);
    this->applyOp(Lcorr, a_correction, true);

    for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit){
      Lcorr[dit()]        -= a_residual[dit()];
      Lcorr[dit()]        *= m_relCoef[dit()];
      
      a_correction[dit()] -= Lcorr[dit()];
    }
  }
}

void EBHelmholtzOp::relaxGauSai(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations){
  LevelData<EBCellFAB> Lcorr;
  this->create(Lcorr, a_residual);

  for (int iter = 0; iter < a_iterations; iter++){
    for (int icolor = 0; icolor < m_colors.size(); icolor++){
      this->homogeneousCFInterp(a_correction);
      this->applyOp(Lcorr, a_correction, true);
      this->gsrbColor(a_correction, Lcorr, a_residual, m_colors[icolor]);
    }
  }
}

void EBHelmholtzOp::gsrbColor(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_Lphi, const LevelData<EBCellFAB>& a_rhs, const IntVect& a_color){
  const DisjointBoxLayout& dbl = a_phi.disjointBoxLayout();
  
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box box = dbl[dit()];

    BaseFab<Real>& phiReg        = a_phi    [dit()].getSingleValuedFAB();
    const BaseFab<Real>& LphiReg = a_Lphi   [dit()].getSingleValuedFAB();
    const BaseFab<Real>& rhsReg  = a_rhs    [dit()].getSingleValuedFAB();
    const BaseFab<Real>& relReg  = m_relCoef[dit()].getSingleValuedFAB();

    // Regular cells (well, plus whatever is not multi-valued)
    IntVect loIV = box.smallEnd();
    IntVect hiIV = box.bigEnd();
    for (int dir = 0; dir < SpaceDim; dir++){
      if(loIV[dir] % 2 != a_color[dir]) loIV[dir]++;
    }

    if(loIV <= hiIV){
      const Box colorBox(loIV, hiIV);

      // Note - replace this with our own kernel so we can develop code without coupling too hard to Chombo. 
      FORT_GSRBEBCO(CHF_FRA1(phiReg, m_comp),
		    CHF_CONST_FRA1(LphiReg, m_comp),
		    CHF_CONST_FRA1(rhsReg, m_comp),
		    CHF_CONST_FRA1(relReg, m_comp),
		    CHF_BOX(colorBox));
    }
    

    // Fortran took care of the irregular cells but multi-cells still remain
    VoFIterator& vofit = m_vofIterMulti[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const IntVect&  iv  = vof.gridIndex();

      bool doThisCell = true;
      for(int dir = 0; dir < SpaceDim; dir++){
	if(iv[dir]%2 != a_color[dir]) doThisCell = false;
      }

      if(doThisCell){
	const Real lambda = m_relCoef[dit()](vof, m_comp);
	const Real Lphi   = a_Lphi[dit()](vof, m_comp);
	const Real rhs    = a_rhs[dit()](vof, m_comp);
	const Real resid  = rhs - Lphi;

	a_phi[dit()](vof, m_comp) += lambda*resid;
      }
    }
  }
}

void EBHelmholtzOp::relaxGSRBFast(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations){
  MayDay::Abort("EBHelmholtzOp::relaxGauSai - not implemented");
}

void EBHelmholtzOp::computeAlphaWeight(){
  for (DataIterator dit(m_eblg.getDBL().dataIterator()); dit.ok(); ++dit){
    VoFIterator& vofit = m_vofIterIrreg[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      const Real volFrac = m_eblg.getEBISL()[dit()].volFrac(vof);
      const Real Aco     = (*m_Acoef)[dit()](vof, m_comp);

      m_alphaDiagWeight[dit()](vof, m_comp) = volFrac*Aco;
    }
  }
}

void EBHelmholtzOp::computeRelaxationCoefficient(){
  for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit){
    const Box cellBox = m_eblg.getDBL()[dit()];

    // Set relaxation coefficient = aco*alpha
    m_relCoef[dit()].setVal(0.0);
    m_relCoef[dit()] += (*m_Acoef)[dit()];
    m_relCoef[dit()] *= m_alpha;

    // Add in the diagonal term for the variable-coefficient Laplacian operator
    BaseFab<Real>& regRel = m_relCoef[dit()].getSingleValuedFAB();
    for (int dir = 0; dir < SpaceDim; dir++){

      // This adds -beta*(bcoef(loFace) + bcoef(hiFace))/dx^2 to the relaxation term.
      BaseFab<Real>& regBcoDir = (*m_Bcoef)[dit()][dir].getSingleValuedFAB();
      FORT_ADDBCOTERMTOINVRELCOEF(CHF_FRA1(regRel, m_comp),
				  CHF_CONST_FRA1(regBcoDir, m_comp),
				  CHF_CONST_REAL(m_beta),
				  CHF_CONST_REAL(m_dx),
				  CHF_CONST_INT(dir),
				  CHF_BOX(cellBox));
    }

    // Invert the relaxation coefficient (in the irregular cells)
    FORT_INVERTRELAXATIONCOEFFICIENT(CHF_FRA1(regRel, m_comp),
				     CHF_BOX(cellBox));

    // Do the same for the irregular cells.
    VoFIterator& vofit = m_vofIterIrreg[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      const Real alphaWeight = m_alpha * m_alphaDiagWeight[dit()](vof, m_comp);
      const Real  betaWeight = m_beta  * m_betaDiagWeight [dit()](vof, m_comp);

      m_relCoef[dit()](vof, m_comp) = 1./(alphaWeight + betaWeight);
    }
  }

  // Add safety factors for relaxations 
  switch(m_relaxationMethod){
  case RelaxationMethod::Jacobi: 
    this->scale(m_relCoef, 0.5); 
    break;
  default:
    break;
  }
}

VoFStencil EBHelmholtzOp::getFaceCenterFluxStencil(const FaceIndex& a_face, const DataIndex& a_dit) const {
  VoFStencil fluxStencil;
  
  if(!a_face.isBoundary()){ // BC handles the boundary fluxes. 
    fluxStencil.add(a_face.getVoF(Side::Hi),  1.0/m_dx);
    fluxStencil.add(a_face.getVoF(Side::Lo), -1.0/m_dx);
    fluxStencil *= (*m_Bcoef)[a_dit][a_face.direction()](a_face, m_comp);
  }

  return fluxStencil;
}

VoFStencil EBHelmholtzOp::getFaceCentroidFluxStencil(const FaceIndex& a_face, const DataIndex& a_dit) const {
  VoFStencil fluxStencil;

  const FaceStencil& interpolationStencil = m_interpStencil[a_face.direction()][a_dit](a_face, m_comp);

  for (int i = 0; i < interpolationStencil.size(); i++){
    const FaceIndex& iface = interpolationStencil.face(i);
    const Real& iweight    = interpolationStencil.weight(i);

    VoFStencil fluxCenterStencil = this->getFaceCenterFluxStencil(iface, a_dit);
    fluxCenterStencil *= iweight;

    fluxStencil += fluxCenterStencil;
  }

  return fluxStencil;
}

VoFStencil EBHelmholtzOp::getDivFStencil(const VolIndex& a_vof, const DataIndex& a_dit) const {
  VoFStencil divStencil;

  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];

  for (int dir = 0; dir < SpaceDim; dir++){
    for (SideIterator sit; sit.ok(); ++sit){
      const int isign = sign(sit());
      
      Vector<FaceIndex> faces = ebisbox.getFaces(a_vof, dir, sit());
      for (auto f : faces.stdVector()){
	VoFStencil centroidFluxStencil = this->getFaceCentroidFluxStencil(f, a_dit); // This includes multiplication by bco. 
	centroidFluxStencil *= sign(sit())*ebisbox.areaFrac(f)/m_dx;

	divStencil += centroidFluxStencil;
      }
    }
  }

  return divStencil;
}

void EBHelmholtzOp::getFaceCentroidFlux(EBFluxFAB&       a_fluxCentroid,
					const EBCellFAB& a_phi,
					const Box&       a_cellBox,
					const DataIndex& a_dit){
  for (int dir = 0; dir < SpaceDim; dir++){
    EBFaceFAB& fluxCentroid = a_fluxCentroid[dir];
    this->getFaceCentroidFlux(fluxCentroid, a_phi, a_cellBox, a_dit, dir);
  }
}

void EBHelmholtzOp::getFaceCentroidFlux(EBFaceFAB&       a_fluxCentroid,
					const EBCellFAB& a_phi,
					const Box&       a_cellBox,
					const DataIndex& a_dit,
					const int        a_dir){

  
  this->getFaceCenteredFlux(a_fluxCentroid, a_phi, a_cellBox, a_dit, a_dir);
  //  this->interpolateFluxes(  a_fluxCentroid,        a_cellBox, a_dit, a_dir);
}

void EBHelmholtzOp::getFaceCenteredFlux(EBFaceFAB&       a_fluxCenter,
					const EBCellFAB& a_phi,
					const Box&       a_cellBox,
					const DataIndex& a_dit,
					const int        a_dir){

  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const EBGraph& ebgraph = ebisbox.getEBGraph();
    
  // Make a face box which has no boundary faces.
  // R.M comment: I don't think this is necessary....
  Box faceBox = a_cellBox;
  faceBox.grow(a_dir, 1);
  faceBox &= m_eblg.getDomain();
  faceBox.grow(a_dir, -1);
  faceBox.surroundingNodes(a_dir);

  BaseFab<Real>& regFlux       = a_fluxCenter.getSingleValuedFAB();
  const BaseFab<Real>& regPhi  = a_phi.getSingleValuedFAB();
  const BaseFab<Real>& regBco  = (*m_Bcoef)[a_dit][a_dir].getSingleValuedFAB();

  FORT_GETINTERIORREGFLUX(CHF_FRA1(regFlux, m_comp),
			   CHF_CONST_FRA1(regPhi, m_comp),
			   CHF_CONST_FRA1(regBco, m_comp),
			   CHF_CONST_INT(a_dir),
			   CHF_CONST_REAL(m_dx),
			   CHF_BOX(faceBox));
  
  // Do the irregular cells
  const IntVectSet irregIVS = ebisbox.getIrregIVS(a_cellBox);
  for (FaceIterator faceIt(irregIVS, ebgraph, a_dir, FaceStop::SurroundingNoBoundary); faceIt.ok(); ++faceIt){
    const FaceIndex& face = faceIt();

    const Real phiHi = a_phi(face.getVoF(Side::Hi), m_comp);
    const Real phiLo = a_phi(face.getVoF(Side::Lo), m_comp);
    const Real bco   = (*m_Bcoef)[a_dit][a_dir](face, m_comp);

    a_fluxCenter(face, m_comp) = bco*(phiHi-phiLo)/m_dx;
  }

  a_fluxCenter *= m_beta;
}

void EBHelmholtzOp::interpolateFluxes(EBFaceFAB& a_flux, const Box& a_cellBox, const DataIndex& a_dit, const int a_dir){
  const BaseIFFAB<FaceStencil>& faceStencils = m_interpStencil[a_dir][a_dit];

  BaseIFFAB<Real>& interpol = m_interpolant[a_dir][a_dit];
  interpol.setVal(0.0);

  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const EBGraph& ebgraph = ebisbox.getEBGraph();

  const IntVectSet ivs = ebisbox.getIrregIVS(a_cellBox);

  // Interpolate flux to centroid.
  for (FaceIterator faceIt(ivs, ebgraph, a_dir, FaceStop::SurroundingNoBoundary); faceIt.ok(); ++faceIt){
    const FaceIndex& face      = faceIt();
    const FaceStencil& stencil = faceStencils(face, m_comp);


    Real centroidFlux = 0.0;
    for (int i = 0; i < stencil.size(); i++){
      centroidFlux += stencil.weight(i)*a_flux(stencil.face(i), m_comp);
    }

    interpol(face, m_comp) = centroidFlux;
  }


  // Copy the centroid flux back to the data holder.
  for (FaceIterator faceIt(ivs, ebgraph, a_dir, FaceStop::SurroundingNoBoundary); faceIt.ok(); ++faceIt){
    const FaceIndex& face = faceIt();
    a_flux(face, m_comp) = interpol(face, m_comp);
  }
}

void EBHelmholtzOp::incrementFRCoar(const LevelData<EBCellFAB>& a_phi){
  CH_assert(m_hasFine);

  const Real scale = 1.0;

  for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit){
    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flux   = m_flux[dit()][dir];

      const Box cellBox = m_eblg.getDBL()[dit()];
      for (SideIterator sit; sit.ok(); ++sit){

	// Computes fluxes for all faces oriented in +/- dir, but which are not boundary faces. Restrict the
	// calculation to stripBox. Must do this for all faces because the CF interface does not align with a coarse box. 
	this->getFaceCentroidFlux(flux, a_phi[dit()], cellBox, dit(), dir);

	// Increment flux register, recall that beta and the b-coefficient are included in getFaceCentroidFlux. 
	m_fluxReg->incrementCoarseBoth(flux, scale, dit(), m_interval, dir, sit());
      }
    }
  }
}

void EBHelmholtzOp::incrementFRFine(const LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phi, AMRLevelOp<LevelData<EBCellFAB> >& a_finerOp){
  CH_assert(m_hasFine);
  
  // Doing the nasty here -- phiFine needs new ghost cells.
  LevelData<EBCellFAB>& phiFine = (LevelData<EBCellFAB>&) a_phiFine;
  EBHelmholtzOp& finerOp        = (EBHelmholtzOp&)       (a_finerOp);
  finerOp.inhomogeneousCFInterp(phiFine, a_phi);
  phiFine.exchange();

  LevelData<EBFluxFAB>& fineFlux = finerOp.getFlux();
  const Real scale = 1.0; 

  // Compute
  for (DataIterator dit(m_eblgFine.getDBL()); dit.ok(); ++dit){
    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flux   = fineFlux[dit()][dir];
      
      const Box cellBox = m_eblgFine.getDBL()[dit()];
      for (SideIterator sit; sit.ok(); ++sit){

	// The CF interface always ends at the boundary of a fine-level grid box, so we can run the computation for
	// the subset of cells that align with the CF. 
	Box stripBox = adjCellBox(cellBox, dir, sit(), 1);
	stripBox.shift(dir, -sign(sit()));

	// Computes fluxes for all faces oriented in +/- dir, but which are not boundary faces. 
	finerOp.getFaceCentroidFlux(flux, phiFine[dit()], stripBox, dit(), dir);

	// Increment flux register, recall that beta and the b-coefficient are included in getFaceCentroidFlux. 
	m_fluxReg->incrementFineBoth(flux, scale, dit(), m_interval, dir, sit());
      }
    }
  }
}

void EBHelmholtzOp::reflux(LevelData<EBCellFAB>&              a_Lphi,
			   const LevelData<EBCellFAB>&        a_phiFine,
			   const LevelData<EBCellFAB>&        a_phi,
			   AMRLevelOp<LevelData<EBCellFAB> >& a_finerOp) {

  m_fluxReg->setToZero();
  
  this->incrementFRCoar(a_phi);
  this->incrementFRFine(a_phiFine, a_phi, a_finerOp);

  m_fluxReg->reflux(a_Lphi, m_interval, 1./m_dx);
}

#include <CD_NamespaceFooter.H>
