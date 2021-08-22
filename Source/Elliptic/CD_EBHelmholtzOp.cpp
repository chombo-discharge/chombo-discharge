/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzOp.cpp
  @brief  Implementation of CD_EBHelmholtzOp.H
  @author Robert Marskar
  @todo   See if we can incorporate EBBCs into the centroid stencils!
  @todo   The domain BC classes are wrong when using centroid discretization, and they need to be fixed. 
*/

// Chombo includes
#include <EBLevelDataOps.H>
#include <EBCellFactory.H>

// Our includes
#include <CD_EBHelmholtzOp.H>
#include <CD_LeastSquares.H>
#include <CD_EBHelmholtzOpF_F.H>
#include <CD_NamespaceHeader.H>

constexpr int EBHelmholtzOp::m_nComp;
constexpr int EBHelmholtzOp::m_comp;

EBHelmholtzOp::EBHelmholtzOp(const Location::Cell                               a_dataLocation,
			     const EBLevelGrid&                                 a_eblgFine,
			     const EBLevelGrid&                                 a_eblg,
			     const EBLevelGrid&                                 a_eblgCoFi,
			     const EBLevelGrid&                                 a_eblgCoar,
			     const EBLevelGrid&                                 a_eblgCoarMG,
			     const RefCountedPtr<EBMultigridInterpolator>&      a_interpolator,
			     const RefCountedPtr<EBFluxRegister>&               a_fluxReg,
			     const RefCountedPtr<EbCoarAve>&                    a_coarAve,			       
			     const RefCountedPtr<EBHelmholtzDomainBC>&          a_domainBc,
			     const RefCountedPtr<EBHelmholtzEBBC>&              a_ebBc,
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
			     const Smoother&                                    a_smoother) :
  LevelTGAHelmOp<LevelData<EBCellFAB>, EBFluxFAB>(false), // Time-independent
  m_dataLocation(a_dataLocation),
  m_eblgFine(),
  m_eblg(a_eblg),
  m_eblgCoFi(),
  m_eblgCoar(),
  m_eblgCoarMG(),
  m_interpolator(a_interpolator),
  m_coarAve(a_coarAve),
  m_fluxReg(a_fluxReg),
  m_domainBc(a_domainBc),
  m_ebBc(a_ebBc),
  m_probLo(a_probLo),
  m_dx(a_dx),
  m_refToCoar(a_hasCoar ? a_refToCoar : 1),
  m_refToFine(a_hasFine ? a_refToFine : 1),
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
  m_smoother(a_smoother) {

  // Default settings. Always solve for comp = 0. If you want something different, copy your
  // input two different data holders before you use AMRMultiGrid. 
  m_turnOffBCs = false;
  m_interval   = Interval(m_comp, m_comp);

  if(m_hasFine){
    m_eblgFine = a_eblgFine;
  }

  // Issue a warning if there are
  if(m_dataLocation == Location::Cell::Centroid && m_hasCoar){
    const int numFill = m_interpolator->getGhostCF();
    if(numFill < 3) MayDay::Warning("EBHelmholtzOp::EBHelmholtzOp - interpolator is not filling enough ghost cells. Your discretization may be incorrect!");
  }

  if(m_hasCoar){
    m_eblgCoFi = a_eblgCoFi;
    m_eblgCoar = a_eblgCoar;

    m_ebInterp.define(m_eblg.getDBL(),   m_eblgCoar.getDBL(), m_eblg.getEBISL(), m_eblgCoar.getEBISL(), m_eblgCoar.getDomain(),
    		      m_refToCoar, m_nComp, m_eblg.getEBIS(), m_ghostPhi);

    m_ebAverage.define(m_eblg.getDBL(),   m_eblgCoFi.getDBL(), m_eblg.getEBISL(), m_eblgCoar.getEBISL(), m_eblgCoar.getDomain(),
		       m_refToCoar, m_nComp, m_eblg.getEBIS(), m_ghostPhi);
  }

  if(m_hasMGObjects){
    constexpr int mgRef = 2;
    
    m_eblgCoarMG = a_eblgCoarMG;

    m_ebInterpMG.define(m_eblg.getDBL(), m_eblgCoarMG.getDBL(), m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(), m_eblgCoarMG.getDomain(),
    			mgRef, m_nComp, m_eblg.getEBIS(), m_ghostPhi);

    m_ebAverageMG.define(m_eblg.getDBL(), m_eblgCoarMG.getDBL(), m_eblg.getEBISL(), m_eblgCoarMG.getEBISL(), m_eblgCoarMG.getDomain(),
			 mgRef, m_nComp, m_eblg.getEBIS(), m_ghostPhi);
  }

  // Define data holders and stencils
  this->defineStencils();
  }

EBHelmholtzOp::~EBHelmholtzOp(){

}

void EBHelmholtzOp::turnOffBCs(){
  m_turnOffBCs = true;
}

void EBHelmholtzOp::turnOnBCs(){
  m_turnOffBCs = false;
}

LevelData<EBFluxFAB>& EBHelmholtzOp::getFlux(){
  return m_flux;
}

LevelData<EBCellFAB>& EBHelmholtzOp::getRelaxationWeights(){
  return m_relCoef;
}

void EBHelmholtzOp::defineStencils(){
  // Basic defines. 
  EBCellFactory cellFact(m_eblg.getEBISL());
  EBFluxFactory fluxFact(m_eblg.getEBISL());
  
  m_relCoef.define(m_eblg.getDBL(), m_nComp, IntVect::Zero, cellFact); // Don't need ghost cells for this one. 
  m_flux.define(   m_eblg.getDBL(), m_nComp, IntVect::Unit, fluxFact); // One ghost cell needed for face center to face centroid interpolation. 

  m_vofIterIrreg.   define(m_eblg.getDBL());
  m_vofIterMulti.   define(m_eblg.getDBL());
  m_vofIterStenc.   define(m_eblg.getDBL()); 
  m_alphaDiagWeight.define(m_eblg.getDBL());
  m_betaDiagWeight. define(m_eblg.getDBL());
  m_opStencils.     define(m_eblg.getDBL());

  for (int dir = 0; dir < SpaceDim; dir++){
    m_vofIterDomLo       [dir].define(m_eblg.getDBL());
    m_vofIterDomHi       [dir].define(m_eblg.getDBL());
    m_centroidFluxStencil[dir].define(m_eblg.getDBL());
  }

  // Get the "colors" for multi-colored relaxation. 
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

  // Define BC objects. Can't do this in the factory because the BC objects will need the b-coefficient,
  // but the factories won't know about that.
  const int ghostCF = m_hasCoar ? m_interpolator->getGhostCF() : 99;
  m_domainBc->define(m_dataLocation, m_eblg, m_Bcoef,      m_probLo, m_dx);
  m_ebBc    ->define(m_dataLocation, m_eblg, m_BcoefIrreg, m_probLo, m_dx, ghostCF);

  // This contains the part of the eb flux that contains interior cells. 
  const LayoutData<BaseIVFAB<VoFStencil> >& ebFluxStencil = m_ebBc->getKappaDivFStencils();

  // Define everything
  for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit){
    const Box cellBox      = m_eblg.getDBL()[dit()];
    const EBISBox& ebisbox = m_eblg.getEBISL()[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    
    const IntVectSet irregIVS = ebisbox.getIrregIVS  (cellBox);
    const IntVectSet multiIVS = ebisbox.getMultiCells(cellBox);


    // Define the cells where we explicitly store stencils. If we use cell-centered data we only
    // need explicit stencils for kappa*L(phi) on the cut-cells. If this is a centroid-based discretization
    // we also need those stencils for cells that share a grid face with a cut-cell. Those cells will have at least one
    // grid face where we can't use centered differencing. 
    IntVectSet stencIVS = irregIVS;
    if(m_dataLocation == Location::Cell::Centroid){
      for (VoFIterator vofit(irregIVS, ebgraph); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();

	for (int dir = 0; dir < SpaceDim; dir++){
	  for (SideIterator sit; sit.ok(); ++sit){
	    Vector<VolIndex> neighborVoFs = ebisbox.getVoFs(vof, dir, sit(), 1);

	    for (const auto& curVoF : neighborVoFs.stdVector()){
	      if(ebisbox.isIrregular(curVoF.gridIndex())){
		stencIVS |= curVoF.gridIndex();
	      }
	    }
	  }
	}
      }
      
      stencIVS &= cellBox;
    }

    // Define iterators. These iterators run over the irregular cells, multi-valued cells only, and cells where we have explicit
    // stencils for kappa*L(phi). The domain iterators are iterators for cut-cells that neighbor a domain edge (face)
    m_vofIterIrreg[dit()].define(irregIVS, ebgraph);
    m_vofIterMulti[dit()].define(multiIVS, ebgraph);
    m_vofIterStenc[dit()].define(stencIVS, ebgraph);

    for (int dir = 0; dir < SpaceDim; dir++){
      const IntVectSet loIrreg = irregIVS & m_sideBox.at(std::make_pair(dir, Side::Lo));
      const IntVectSet hiIrreg = irregIVS & m_sideBox.at(std::make_pair(dir, Side::Hi));

      m_vofIterDomLo[dir][dit()].define(loIrreg, ebgraph);
      m_vofIterDomHi[dir][dit()].define(hiIrreg, ebgraph);
    }

    // Define data holders. 
    m_opStencils     [dit()].define(stencIVS, ebgraph, m_nComp);
    m_alphaDiagWeight[dit()].define(stencIVS, ebgraph, m_nComp);
    m_betaDiagWeight [dit()].define(stencIVS, ebgraph, m_nComp);

    // The below code may seem intimidating at first. What happens is that we explicitly store stencils for all cells that is either a cut-cell
    // or shares a face with a cut-cell. Now, we have to compute stencils explicitly for this subset of cells, and we need representations both
    // of centroid fluxes, i.e. b*grad(phi) (because of refluxing), and also kappa*div(F). The latter is obviously found by summing the finite
    // volume contributions.
    //
    //                    EB
    //                   /
    //    |------------|/ ----------|------------|
    //    |A           /           B|           C|
    //    |           /|            |            |
    //    |          / |            |     x      |
    //    |         /  |       x    |            |
    //    |        /  x|            |            |    
    //    |-------/----|------------|------------|
    //    |D     /     |           E|           F|
    //    |     /      |            |            |
    //    |    /       |     x      |     x      |
    //    |   /    x   |            |            |
    //    |  /         |            |            |
    //    |-/----------|------------|------------|
    //     /
    //    EB
    //
    // The sketch above represents the principle behind the discretization. All the faces A-B, B-C, D-E, A-D, B-E require
    // explicit stencils for the fluxes through the faces. This is true even though E is a regular cell and B-C is a regular face. 
    //
    // The code computes the required stencils using (essentially) the following procedure:
    //
    //   1. Iterate through all faces in cut-cells, and for a regular cell which shares at least one edge (face in 3D) with a cut-cell.
    //   2. Compute an approximation to the flux on the face centroid
    //   3. Store this flux stencil explicitly (it might be needed for refluxing operations).
    //   4. The face obviously connects two vofs, and if one of those vofs is a part of the subset of cells discussed above, add the flux
    //      contribution to that vof. For example, we need to store a flux stencil for 
    //   5. Irregular cells have an extra face, the EB face. Add the flux contribution from that face into storage for kappa*div(F).
    //
    // Referring to the sketch above, we iterate through cells A, B, D, E and compute the face centroid fluxes for all these cells. This includes
    // e.g. the face connecting E and F. However, since F only has regular faces we don't store the stencil explicitly for that cell.
    //
    BaseIVFAB<VoFStencil>& opStencil = m_opStencils  [dit()];
    VoFIterator& vofitStenc          = m_vofIterStenc[dit()];
    VoFIterator& vofitIrreg          = m_vofIterIrreg[dit()];

    for (vofitStenc.reset(); vofitStenc.ok(); ++vofitStenc){
      opStencil(vofitStenc(), m_comp).clear(); 
    }
    
    for (int dir = 0; dir < SpaceDim; dir++){
      m_centroidFluxStencil[dir][dit()].define(stencIVS, ebgraph, dir, m_nComp);
      
      BaseIFFAB<VoFStencil>& fluxStencils = m_centroidFluxStencil[dir][dit()];

      // 1. 
      for (FaceIterator faceIt(stencIVS, ebgraph, dir, FaceStop::SurroundingNoBoundary); faceIt.ok(); ++faceIt){ 
	const FaceIndex& face = faceIt();

	if(!face.isBoundary()){
	  const VolIndex vofLo = face.getVoF(Side::Lo);
	  const VolIndex vofHi = face.getVoF(Side::Hi);

	  // 2.
	  const VoFStencil fluxSten = this->getFaceCentroidFluxStencil(face, dit()); 

	  // 3.
	  fluxStencils(face, m_comp) = fluxSten; 

	  VoFStencil loKappaDivFSten = fluxSten;
	  VoFStencil hiKappaDivFSten = fluxSten;

	  loKappaDivFSten *=  ebisbox.areaFrac(face)/m_dx; // Sign explanation. For vofLo, faceIt() is the face on the "high" side and
	  hiKappaDivFSten *= -ebisbox.areaFrac(face)/m_dx; // vice versa for vofHi.

	  // 4. Note that we prune storage here. 
	  if(stencIVS.contains(vofLo.gridIndex())) opStencil(vofLo, m_comp) += loKappaDivFSten;
	  if(stencIVS.contains(vofHi.gridIndex())) opStencil(vofHi, m_comp) += hiKappaDivFSten;
	}
      }
    }

    // 5. Add contributions to the operator from the EB faces. 
    for (vofitIrreg.reset(); vofitIrreg.ok(); ++vofitIrreg){
      const VolIndex& vof = vofitIrreg();
      opStencil(vof, m_comp) += ebFluxStencil[dit()](vof, m_comp);
    }

    // Compute relaxation factor. Adjust the weight with domain boundary faces. 
    for (vofitStenc.reset(); vofitStenc.ok(); ++vofitStenc){
      const VolIndex& vof = vofitStenc();
      const IntVect iv    = vof.gridIndex();

      VoFStencil& curStencil = opStencil(vof, m_comp);

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
  }
    
  // Compute the alpha-weight and relaxation coefficient. 
  this->computeAlphaWeight();
  this->computeRelaxationCoefficient();
}

void EBHelmholtzOp::setAlphaAndBeta(const Real& a_alpha, const Real& a_beta) {
  m_alpha = a_alpha;
  m_beta  = a_beta;

  // When we change alpha and beta we need to recompute relaxation coefficients...
  this->computeAlphaWeight(); 
  this->computeRelaxationCoefficient();
}

void EBHelmholtzOp::residual(LevelData<EBCellFAB>& a_residual, const LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_rhs, bool a_homogeneousPhysBC) {
  this->applyOp(a_residual, a_phi, nullptr, a_homogeneousPhysBC, true); // residual = L(phi)
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
  Real globalNorm = 0.0;
  Real localNorm  = 0.0;
  
  for (DataIterator dit = a_rhs.dataIterator(); dit.ok(); ++dit){
    const EBCellFAB& rhs       = a_rhs[dit()];
    const FArrayBox& regRhs    = rhs.getFArrayBox();
    const EBISBox& ebisbox     = rhs.getEBISBox();
    const EBGraph& ebgraph     = ebisbox.getEBGraph();
    const Box box              = a_rhs.disjointBoxLayout()[dit()];
    const IntVectSet& irregIVS = ebisbox.getIrregIVS(box);

    // Can replace with Fortran kernel if we need to.
    for (BoxIterator bit(box); bit.ok(); ++bit){
      const IntVect iv = bit();
      if(ebisbox.isRegular(iv)) localNorm = std::max(localNorm, std::abs(regRhs(iv, m_comp)));
    }

    for (VoFIterator vofit(irregIVS, ebgraph); vofit.ok(); ++vofit){
      localNorm = std::max(localNorm, std::abs(rhs(vofit(), m_comp)));
    }
  }

#ifdef CH_MPI
  MPI_Allreduce(&localNorm, &globalNorm, 1, MPI_CH_REAL, MPI_MAX, Chombo_MPI::comm);
#else
  globalNorm = localNorm;
#endif
  
  return globalNorm;
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
  if(m_hasFine){
    EBHelmholtzOp* fineOp = (EBHelmholtzOp*) a_finerOp;
    fineOp->coarsen((LevelData<EBCellFAB>&) a_phi, a_phiFine);
  }
  
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

  constexpr bool homogeneousPhysBC = true;
  constexpr bool homogeneousCFBC   = false;

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

  constexpr bool homogeneousPhysBC = true;
  constexpr bool homogeneousCFBC   = false;
  
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
    this->interpolateCF(phi, a_phiCoar, a_homogeneousCFBC);
  }
  
  phi.exchange();
  
  const DisjointBoxLayout& dbl = a_Lphi.disjointBoxLayout();
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box cellBox = dbl[dit()];

    this->applyOp(a_Lphi[dit()], phi[dit()], cellBox, dit(), a_homogeneousPhysBC);
  }
}

void EBHelmholtzOp::applyOp(EBCellFAB& a_Lphi, EBCellFAB& a_phi, const Box& a_cellBox, const DataIndex& a_dit, const bool a_homogeneousPhysBC){
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];

  if(!ebisbox.isAllCovered()){
    this->applyOpRegular  (a_Lphi, a_phi, a_cellBox, a_dit, a_homogeneousPhysBC);
    this->applyOpIrregular(a_Lphi, a_phi, a_cellBox, a_dit, a_homogeneousPhysBC);
  }
}

void EBHelmholtzOp::applyOpRegular(EBCellFAB& a_Lphi, EBCellFAB& a_phi, const Box& a_cellBox, const DataIndex& a_dit, const bool a_homogeneousPhysBC){
  // Fill a_phi such that centered differences pushes in the domain flux. 
  this->applyDomainFlux(a_phi, a_cellBox, a_dit, a_homogeneousPhysBC);

  BaseFab<Real>& Lphi      = a_Lphi.getSingleValuedFAB();
  const BaseFab<Real>& phi = a_phi.getSingleValuedFAB();
  const BaseFab<Real>& aco = (*m_Acoef)[a_dit].getSingleValuedFAB();

  // It's an older code, sir, but it checks out.
  const BaseFab<Real> dummy(Box(IntVect::Zero, IntVect::Zero), 1);
  const BaseFab<Real>* bcoef[3];

  for (int dir = 0; dir < 3; dir++){
    if(dir >= SpaceDim){
      bcoef[dir] = &dummy;
    }
    else{
      bcoef[dir] = &(*m_Bcoef)[a_dit][dir].getSingleValuedFAB();
    }
  }

  // Fill the ghost cells outside the domain so that we don't get an extra flux from the domain side. 
  FORT_HELMHOLTZINPLACE(CHF_FRA1(Lphi, m_comp),
			CHF_CONST_FRA1(phi, m_comp),
			CHF_CONST_FRA1(aco, m_comp),
			CHF_CONST_FRA1((*bcoef[0]), m_comp),
			CHF_CONST_FRA1((*bcoef[1]), m_comp),
			CHF_CONST_FRA1((*bcoef[2]), m_comp), // Not used for 2D. 
			CHF_CONST_REAL(m_dx),
			CHF_CONST_REAL(m_alpha),
			CHF_CONST_REAL(m_beta),
			CHF_BOX(a_cellBox));
}

void EBHelmholtzOp::applyDomainFlux(EBCellFAB& a_phi, const Box& a_cellBox, const DataIndex& a_dit, const bool a_homogeneousPhysBC){
  // TLDR: We compute the flux on the domain edges and store it in a cell-centered box. We then monkey with the ghost cells
  //       so that centered differences on the edge cells inject said flux.
  
  for (int dir = 0; dir < SpaceDim; dir++){

    BaseFab<Real>& phiFAB = a_phi.getSingleValuedFAB();

    Box loBox;
    Box hiBox;
    int hasLo;
    int hasHi;
    EBArith::loHi(loBox, hasLo, hiBox, hasHi, m_eblg.getDomain(), a_cellBox, dir);
    
    if(hasLo == 1){
      const int side = -1;

      // Get domain flux. 
      FArrayBox faceFlux(loBox, m_nComp);
      m_domainBc->getFaceFlux(faceFlux, a_phi.getSingleValuedFAB(), dir, Side::Lo, a_dit, a_homogeneousPhysBC);

      // loBox is cell-centered interior region. 
      Box ghostBox = loBox;    
      ghostBox.shift(dir, -1);
      
      BaseFab<Real>& bco = (*m_Bcoef)[a_dit][dir].getSingleValuedFAB();
      FORT_HELMHOLTZAPPLYDOMAINFLUX(CHF_FRA1(phiFAB, m_comp),
      				    CHF_CONST_FRA1(faceFlux, m_comp),
      				    CHF_CONST_FRA1(bco,m_comp),
      				    CHF_CONST_REAL(m_dx),
				    CHF_CONST_INT(side),
				    CHF_CONST_INT(dir),
				    CHF_BOX(ghostBox));
    }

    if(hasHi == 1){
      const int side = 1;

      // Get domain flux. 
      FArrayBox faceFlux(hiBox, m_nComp);
      m_domainBc->getFaceFlux(faceFlux, a_phi.getSingleValuedFAB(), dir, Side::Hi, a_dit, a_homogeneousPhysBC);

      // hiBox is cell-centered interior region. We will monkey with the ghost cells so that centered differences in this
      // region injects the domain flux. 
      Box ghostBox = hiBox;
      ghostBox.shift(dir, 1);

      const BaseFab<Real>& bco = (*m_Bcoef)[a_dit][dir].getSingleValuedFAB();
      FORT_HELMHOLTZAPPLYDOMAINFLUX(CHF_FRA1(phiFAB, m_comp),
      				    CHF_CONST_FRA1(faceFlux, m_comp),
      				    CHF_CONST_FRA1(bco,m_comp),
      				    CHF_CONST_REAL(m_dx),
				    CHF_CONST_INT(side),
				    CHF_CONST_INT(dir),
				    CHF_BOX(ghostBox));
    }
  }
}

void EBHelmholtzOp::applyOpIrregular(EBCellFAB& a_Lphi, const EBCellFAB& a_phi, const Box& a_cellBox, const DataIndex& a_dit, const bool a_homogeneousPhysBC){
  // Apply the operator in all cells where we needed an explicit stencil. 
  VoFIterator& vofit = m_vofIterStenc[a_dit];
  for (vofit.reset(); vofit.ok(); ++vofit){
    const VolIndex& vof     = vofit();
    const VoFStencil& stenc = m_opStencils[a_dit](vof, m_comp);                // = Finite volume stencil representation of Int[B*grad(phi)]dA
    const Real& alphaDiag   = m_alpha * m_alphaDiagWeight[a_dit](vof, m_comp); // = kappa * alpha * aco (m_alphaDiagWeight holds kappa* aco)
    
    a_Lphi(vof, m_comp) = alphaDiag * a_phi(vof, m_comp);
    
    for (int i = 0; i < stenc.size(); i++){
      const VolIndex& ivof = stenc.vof(i);
      const Real& iweight  = stenc.weight(i);

      a_Lphi(vof, m_comp) += m_beta * iweight * a_phi(ivof, m_comp); // Note that bco is a part of the stencil weight. 
    }
  }
  
  m_ebBc->applyEBFlux(m_vofIterIrreg[a_dit], a_Lphi, a_phi, a_dit, m_beta, a_homogeneousPhysBC);

  // Do irregular faces on domain sides. m_domainBc should give the centroid-centered flux so we don't do interpolations here. 
  for (int dir = 0; dir < SpaceDim; dir++){
    // Lo side.
    VoFIterator& vofitLo = m_vofIterDomLo[dir][a_dit];    
    for (vofitLo.reset(); vofitLo.ok(); ++vofitLo){
      const VolIndex& vof = vofitLo();
      
      const Real flux = m_domainBc->getFaceFlux(vof, a_phi, dir, Side::Lo, a_dit, a_homogeneousPhysBC);

      a_Lphi(vof, m_comp) -= flux*m_beta/m_dx;
    }

    // Hi side. 
    VoFIterator& vofitHi = m_vofIterDomHi[dir][a_dit];
    for (vofitHi.reset(); vofitHi.ok(); ++vofitHi){
      const VolIndex& vof = vofitHi();

      const Real flux = m_domainBc->getFaceFlux(vof, a_phi, dir, Side::Hi, a_dit, a_homogeneousPhysBC);      

      a_Lphi(vof, m_comp) += flux*m_beta/m_dx;
    }
  }
}

void EBHelmholtzOp::diagonalScale(LevelData<EBCellFAB> & a_rhs, bool a_kappaWeighted){

  // Scale by volume fraction if asked. 
  if(a_kappaWeighted) EBLevelDataOps::kappaWeight(a_rhs);

  // Scale by a-coefficient and alpha, too.
  for (DataIterator dit(a_rhs.dataIterator()); dit.ok(); ++dit){
    a_rhs[dit()] *= (*m_Acoef)[dit()];
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
  MayDay::Warning("EBHelmholtzOp::fillGrad - not implemented (yet)");
}

void EBHelmholtzOp::getFlux(EBFluxFAB&                  a_flux,
			    const LevelData<EBCellFAB>& a_data,
			    const Box&                  a_grid,
			    const DataIndex&            a_dit,
			    Real                        a_scale) {
  MayDay::Warning("EBHelmholtzOp::getFlux - not implemented (yet)");
}

void EBHelmholtzOp::homogeneousCFInterp(LevelData<EBCellFAB>& a_phi){
  if(m_hasCoar) m_interpolator->coarseFineInterpH(a_phi, m_interval);
}

void EBHelmholtzOp::inhomogeneousCFInterp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar){
  if(m_hasCoar) m_interpolator->coarseFineInterp(a_phiFine, a_phiCoar, m_interval);
}

void EBHelmholtzOp::interpolateCF(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>* a_phiCoar, const bool a_homogeneousCFBC){
  if(m_hasCoar){
    if(a_homogeneousCFBC){
      this->homogeneousCFInterp(a_phiFine);
    }
    else{
      if(a_phiCoar == nullptr) MayDay::Error("EBHelmholtzOp::interpolateCF -- calling inhomogeneousCFInterp with nullptr coarse is an error.");
      this->inhomogeneousCFInterp(a_phiFine, *a_phiCoar);
    }
  }
}

void EBHelmholtzOp::relax(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, int a_iterations){
  switch(m_smoother){
  case Smoother::NoRelax: // Don't know why you would do this. 
    break;
  case Smoother::PointJacobi:
    this->relaxPointJacobi(a_correction, a_residual, a_iterations);
    break;
  case Smoother::GauSaiRedBlack:
    this->relaxGSRedBlack(a_correction, a_residual, a_iterations);
    break;
  case Smoother::GauSaiMultiColor:
    this->relaxGSMultiColor(a_correction, a_residual, a_iterations);
    break;
  default:
    MayDay::Error("EBHelmholtzOp::relax - bogus relaxation method requested");
  };
}

void EBHelmholtzOp::relaxPointJacobi(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations){
  LevelData<EBCellFAB> Lcorr;
  this->create(Lcorr, a_residual);

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  for (int iter = 0; iter < a_iterations; iter++){
    a_correction.exchange();

    this->homogeneousCFInterp(a_correction);

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      const Box cellBox = dbl[dit()];
      this->pointJacobiKernel(Lcorr[dit()], a_correction[dit()], a_residual[dit()], cellBox, dit());
    }
  }
}

void EBHelmholtzOp::pointJacobiKernel(EBCellFAB& a_Lcorr, EBCellFAB& a_correction, const EBCellFAB& a_residual, const Box& a_cellBox, const DataIndex& a_dit) {
  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];

  if(!ebisbox.isAllCovered()){
    this->applyOp(a_Lcorr, a_correction, a_cellBox, a_dit, true);
  
    a_Lcorr      -= a_residual;
    a_Lcorr      *= m_relCoef[a_dit];
    a_Lcorr      *= 0.5;
    a_correction -= a_Lcorr;
  }
}

void EBHelmholtzOp::relaxGSRedBlack(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations){
  LevelData<EBCellFAB> Lcorr;
  this->create(Lcorr, a_residual);

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  for (int iter = 0; iter < a_iterations; iter++){
    for (int redBlack = 0; redBlack <= 1; redBlack++){
      a_correction.exchange();
      
      this->homogeneousCFInterp(a_correction);
      
      for (DataIterator dit(dbl); dit.ok(); ++dit){
	const Box cellBox = dbl[dit()];
	
	this->gauSaiRedBlackKernel(Lcorr[dit()], a_correction[dit()], a_residual[dit()], cellBox, dit(), redBlack);
      }
    }
  }
}

void EBHelmholtzOp::gauSaiRedBlackKernel(EBCellFAB& a_Lcorr, EBCellFAB& a_corr, const EBCellFAB& a_resid, const Box& a_cellBox, const DataIndex& a_dit, const int& a_redBlack){
  const EBISBox& ebisbox   = m_eblg.getEBISL()[a_dit];
  const EBCellFAB& relCoef = m_relCoef[a_dit];
  
  if(!ebisbox.isAllCovered()){
    this->applyOp(a_Lcorr, a_corr, a_cellBox, a_dit, true);

    BaseFab<Real>& phiReg        = a_corr .getSingleValuedFAB();
    const BaseFab<Real>& LphiReg = a_Lcorr.getSingleValuedFAB();
    const BaseFab<Real>& rhsReg  = a_resid.getSingleValuedFAB();
    const BaseFab<Real>& relReg  = relCoef.getSingleValuedFAB();

    FORT_HELMHOLTZGAUSAIREDBLACK(CHF_FRA1(phiReg, m_comp),
				 CHF_CONST_FRA1(LphiReg, m_comp),
				 CHF_CONST_FRA1(rhsReg, m_comp),
				 CHF_CONST_FRA1(relReg, m_comp),
				 CHF_CONST_INT(a_redBlack),
				 CHF_BOX(a_cellBox));

    // Fortran took care of the irregular cells but multi-cells still remain
    VoFIterator& vofit = m_vofIterMulti[a_dit];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const IntVect&  iv  = vof.gridIndex();


      const bool doThisCell = std::abs((iv.sum() + a_redBlack)%2) == 0;

      if(doThisCell){
	const Real lambda = relCoef(vof, m_comp);
	const Real Lcorr  = a_Lcorr(vof, m_comp);
	const Real rhs    = a_resid(vof, m_comp);
	
	const Real resid  = rhs - Lcorr;

	a_corr(vof, m_comp) += lambda*resid;
      }
    }
  }
}

void EBHelmholtzOp::relaxGSMultiColor(LevelData<EBCellFAB>& a_correction, const LevelData<EBCellFAB>& a_residual, const int a_iterations){
  LevelData<EBCellFAB> Lcorr;
  this->create(Lcorr, a_residual);

  const DisjointBoxLayout& dbl = m_eblg.getDBL();

  for (int iter = 0; iter < a_iterations; iter++){
    for (int icolor = 0; icolor < m_colors.size(); icolor++){
      a_correction.exchange();
      
      this->homogeneousCFInterp(a_correction);
      
      for (DataIterator dit(dbl); dit.ok(); ++dit){
	const Box cellBox = dbl[dit()];

	this->gauSaiMultiColorKernel(Lcorr[dit()], a_correction[dit()], a_residual[dit()], cellBox, dit(), m_colors[icolor]);
      }
    }
  }
}

void EBHelmholtzOp::gauSaiMultiColorKernel(EBCellFAB&       a_Lcorr,
					   EBCellFAB&       a_corr,
					   const EBCellFAB& a_resid,
					   const Box&       a_cellBox,
					   const DataIndex& a_dit,
					   const IntVect&   a_color){
  const EBISBox& ebisbox   = m_eblg.getEBISL()[a_dit];
  const EBCellFAB& relCoef = m_relCoef[a_dit];
  
  if(!ebisbox.isAllCovered()){
    this->applyOp(a_Lcorr, a_corr, a_cellBox, a_dit, true);

    BaseFab<Real>& phiReg        = a_corr. getSingleValuedFAB();
    const BaseFab<Real>& LphiReg = a_Lcorr.getSingleValuedFAB();
    const BaseFab<Real>& rhsReg  = a_resid.getSingleValuedFAB();
    const BaseFab<Real>& relReg  = relCoef.getSingleValuedFAB();

    // Regular cells (well, plus whatever is not multi-valued)
    IntVect loIV = a_cellBox.smallEnd();
    IntVect hiIV = a_cellBox.bigEnd();
    for (int dir = 0; dir < SpaceDim; dir++){
      if(loIV[dir] % 2 != a_color[dir]) loIV[dir]++;
    }

    if(loIV <= hiIV){
      const Box colorBox(loIV, hiIV);

      FORT_HELMHOLTZGAUSAICOLOR(CHF_FRA1(phiReg, m_comp),
				CHF_CONST_FRA1(LphiReg, m_comp),
				CHF_CONST_FRA1(rhsReg, m_comp),
				CHF_CONST_FRA1(relReg, m_comp),
				CHF_BOX(colorBox));
    }
    

    // Fortran took care of the irregular cells but multi-cells still remain
    VoFIterator& vofit = m_vofIterMulti[a_dit];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const IntVect&  iv  = vof.gridIndex();

      bool doThisCell = true;
      for(int dir = 0; dir < SpaceDim; dir++){
	if(iv[dir]%2 != a_color[dir]) doThisCell = false;
      }

      if(doThisCell){
	const Real lambda = relCoef(vof, m_comp);
	const Real Lcorr  = a_Lcorr(vof, m_comp);
	const Real rhs    = a_resid(vof, m_comp);
	const Real resid  = rhs - Lcorr;

	a_corr(vof, m_comp) += lambda*resid;
      }
    }
  }
}

void EBHelmholtzOp::computeAlphaWeight(){
  for (DataIterator dit(m_eblg.getDBL().dataIterator()); dit.ok(); ++dit){
    VoFIterator& vofit = m_vofIterStenc[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      const Real volFrac = m_eblg.getEBISL()[dit()].volFrac(vof);
      const Real Aco     = (*m_Acoef)[dit()](vof, m_comp);

      m_alphaDiagWeight[dit()](vof, m_comp) = volFrac * Aco;
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
    VoFIterator& vofit = m_vofIterStenc[dit()];
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      const Real alphaWeight = m_alpha * m_alphaDiagWeight[dit()](vof, m_comp);
      const Real  betaWeight = m_beta  * m_betaDiagWeight [dit()](vof, m_comp);

      m_relCoef[dit()](vof, m_comp) = 1./(alphaWeight + betaWeight);
    }
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

  if(!a_face.isBoundary()){ // Domain BC classes handle domain faces. 
    const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
    
    const VolIndex loVoF   = a_face.getVoF(Side::Lo);
    const VolIndex hiVoF   = a_face.getVoF(Side::Hi);
    
    const bool irregFace   = ebisbox.isIrregular(loVoF.gridIndex()) || ebisbox.isIrregular(hiVoF.gridIndex());

    // Centered differencing for regular faces. 
    if(!irregFace){
      fluxStencil = this->getFaceCenterFluxStencil(a_face, a_dit); 
    }
    else {
      if(m_dataLocation == Location::Cell::Center){ // Irregular face, but using cell-centered data. Get face-centered stencils and interpolate them to centroids
	const FaceStencil interpolationStencil = EBArith::getInterpStencil(a_face, IntVectSet(), ebisbox, m_eblg.getDomain());

	for (int i = 0; i < interpolationStencil.size(); i++){
	  const FaceIndex& iface = interpolationStencil.face(i);
	  const Real& iweight    = interpolationStencil.weight(i);

	  // Get the face-centered stencil.
	  VoFStencil fluxCenterStencil = this->getFaceCenterFluxStencil(iface, a_dit);
	
	  fluxCenterStencil *= iweight;

	  fluxStencil += fluxCenterStencil;
	}
      }
      else{ // Irregular face, but with centroid-centered data. In this case we reconstruct the gradient using a least squares reconstruction of the solution.
	constexpr int stencilWeight = 2;
	constexpr int stencilRadius = 2;
	constexpr int stencilOrder  = 2;
	
	const VoFStencil gradSten = LeastSquares::getGradSten(a_face,
							      LeastSquares::FaceLocation::Centroid,
							      LeastSquares::CellLocation::Centroid,
							      ebisbox,
							      m_dx,
							      stencilRadius,
							      stencilWeight,
							      stencilOrder,
							      IntVectSet());

	if(gradSten.size() > 0) {
	  fluxStencil  = LeastSquares::projectGradSten(gradSten, BASISREALV(a_face.direction()));
	  fluxStencil *= (*m_Bcoef)[a_dit][a_face.direction()](a_face, m_comp);
	}
	else{
	  MayDay::Warning("EBHelmholtz::getFaceCentroidFluxStencil -- could not find centroid stencil. Maybe your multigrid runs too deep?");
	}
      }
    }
  }
  
  return fluxStencil;
}

void EBHelmholtzOp::getFaceCentroidFlux(EBFaceFAB&       a_fluxCentroid,
					const EBCellFAB& a_phi,
					const Box&       a_cellBox,
					const DataIndex& a_dit,
					const int        a_dir){
  this->computeFaceCenteredFlux(a_fluxCentroid, a_phi, a_cellBox, a_dit, a_dir);
  this->computeFaceCentroidFlux(a_fluxCentroid, a_phi, a_cellBox, a_dit, a_dir);

  a_fluxCentroid *= m_beta;
}

void EBHelmholtzOp::computeFaceCenteredFlux(EBFaceFAB&       a_fluxCenter,
					    const EBCellFAB& a_phi,
					    const Box&       a_cellBox,
					    const DataIndex& a_dit,
					    const int        a_dir){

  const EBISBox& ebisbox = m_eblg.getEBISL()[a_dit];
  const EBGraph& ebgraph = ebisbox.getEBGraph();
    
  // Since a_fluxCenter is used as a basis for interpolation near the cut-cells, we need to fill fluxes on faces in direction that extends one cell
  // outside of the cellbox in directions that are tangential to a_dir. We do not want to include domain faces here since that
  // is taken care of elsewhere.

  // Discard domain faces
  Box compCellBox = a_cellBox;
  compCellBox.grow(a_dir, 1);
  compCellBox &= m_eblg.getDomain();
  compCellBox.grow(a_dir, -1);

  // cellbox -> facebox
  Box faceBox = compCellBox;
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
}

void EBHelmholtzOp::computeFaceCentroidFlux(EBFaceFAB&       a_flux,
					    const EBCellFAB& a_phi,
					    const Box&       a_cellBox,
					    const DataIndex& a_dit,
					    const int        a_dir){

  const BaseIFFAB<VoFStencil>& fluxStencils = m_centroidFluxStencil[a_dir][a_dit];
  const EBGraph& ebgraph                    = fluxStencils.getEBGraph();
  IntVectSet ivs                            = fluxStencils.getIVS();

  ivs &= a_cellBox;

  for (FaceIterator faceIt(ivs, ebgraph, a_dir, FaceStop::SurroundingNoBoundary); faceIt.ok(); ++faceIt){
    const FaceIndex& face  = faceIt();
    const VoFStencil& sten = fluxStencils(face, m_comp);

    a_flux(face, m_comp) = 0.0;
    for (int i = 0; i < sten.size(); i++){
      const VolIndex& ivof = sten.vof(i);
      const Real& iweight  = sten.weight(i);

      a_flux(face, m_comp) += iweight * a_phi(ivof, m_comp);
    }
  }
}

void EBHelmholtzOp::incrementFRCoar(const LevelData<EBCellFAB>& a_phi){
  CH_assert(m_hasFine);

  const Real scale = 1.0;

  LevelData<EBCellFAB>& phiCoar = (LevelData<EBCellFAB>&) a_phi;
  phiCoar.exchange();
  
  for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit){
    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flux   = m_flux[dit()][dir];

      const Box cellBox = m_eblg.getDBL()[dit()];
      for (SideIterator sit; sit.ok(); ++sit){

	// Computes fluxes for all faces oriented in +/- dir, but which are not boundary faces. Must do this for all faces because
	// the CF interface does not align with a coarse box. 
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
	Box stripBox = adjCellBox(cellBox, dir, sit(), 1); // This box is outside the cellBox patch
	stripBox.shift(dir, -sign(sit()));                 // This box is inside the cellBox patch

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

void EBHelmholtzOp::coarsen(LevelData<EBCellFAB>& a_phi, const LevelData<EBCellFAB>& a_phiFine) {
  m_coarAve->average(a_phi, a_phiFine, m_interval);
}

#include <CD_NamespaceFooter.H>
