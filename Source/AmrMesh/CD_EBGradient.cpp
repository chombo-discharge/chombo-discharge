/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBGradient.cpp
  @brief  Implementation of CD_EBGradient.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>
#include <EBArith.H>
#include <EBCellFactory.H>
#include <EBFluxFactory.H>
#include <BaseIVFactory.H>
#include <CH_Timer.H>

// Our includes
#include <CD_Timer.H>
#include <CD_EBGradient.H>
#include <CD_EBGradientF_F.H>
#include <CD_LeastSquares.H>
#include <CD_VofUtils.H>
#include <CD_LoadBalancing.H>
#include <CD_NamespaceHeader.H>

constexpr int EBGradient::m_comp;
constexpr int EBGradient::m_nComp;

EBGradient::EBGradient(const EBLevelGrid& a_eblg,
		       const EBLevelGrid& a_eblgFine,
		       const Real         a_dx,
		       const int          a_refRat,
		       const int          a_order,
		       const int          a_weighting){
  CH_TIME("EBGradient::EBGradient");

  CH_assert(a_order    >  0);
  CH_assert(a_weight   >= 0);  
  CH_assert(a_refRat%2 == 0);

  m_eblg         = a_eblg;
  m_eblgFine     = a_eblgFine;
  m_dataLocation = Location::Cell::Center;
  m_dx           = a_dx;
  m_refRat       = a_refRat;
  m_order        = a_order;
  m_weighting    = a_weighting;
  m_hasEBCF      = false;  

  if(a_eblgFine.isDefined()){
    m_hasFine  = true;
    m_dxFine   = m_dx/m_refRat;
  }
  else{
    m_hasFine  = false;
    m_dxFine   = 1;
  }

  ParmParse pp("EBGradient");
  bool profile     = false;
  bool forceNoEBCF = false;
  pp.query("profile",       profile);
  pp.query("force_no_ebcf", forceNoEBCF);    

  Timer timer("EBGradient::EBGradient");

  // Define the level stencils. These are regular finite-difference stencils. 
  timer.startEvent("Define level stencils");
  this->defineLevelStencils();
  timer.stopEvent("Define level stencils");

  // If there's finer level there might also be an EBCF crossing. Those can be tricky to deal with because
  // the stencil on this level might reach underneath the finer level and obtain bogus data. This loop defines
  // a few masks for figuring out which cells lie on the CF region, and which cells can't be used for finite
  // differencing. Once those masks have been defined, we iterate through them and see if there are cells
  // where we can't find good stencils. If we do find such cells, we trigger m_hasEBCF which will define
  // a new set of grids where we run a least squares procedure. THOSE stencils are "dual-level" stencils
  // that only reach into valid regions (i.e., not covered by a finer grid level). 
  if(m_hasFine && !forceNoEBCF) {

    // Masks with transient lifetimes. 
    BoxLayoutData<FArrayBox> coarMaskCF;
    BoxLayoutData<FArrayBox> coarMaskInvalid;

    BoxLayoutData<FArrayBox> bufferCoarMaskCF;
    BoxLayoutData<FArrayBox> bufferCoarMaskInvalid;        

    // Define masks on the input grids. 
    timer.startEvent("EBCF masks");
    this->defineMasks(coarMaskCF, coarMaskInvalid);
    timer.stopEvent("EBCF masks");

    // Define iterators for the input grid AND define the simplified buffer grids. 
    timer.startEvent("EBCF Iterators");
    this->defineIteratorsEBCF(coarMaskCF, coarMaskInvalid);
    timer.stopEvent("EBCF Iterators");

    if(m_hasEBCF){
      timer.startEvent("EBCF buffers");
      this->defineBuffers(bufferCoarMaskCF, bufferCoarMaskInvalid, coarMaskCF, coarMaskInvalid);
      timer.stopEvent("EBCF buffers");            

      timer.startEvent("EBCF stencils");
      this->defineStencilsEBCF(bufferCoarMaskInvalid);
      timer.stopEvent("EBCF stencils");
    }
  }

  if(profile){
    timer.eventReport(pout(), false);
  }
}

EBGradient::~EBGradient(){
  CH_TIME("EBGradient::~EBGradient");  
}

void EBGradient::computeLevelGradient(LevelData<EBCellFAB>&       a_gradient,
				      const LevelData<EBCellFAB>& a_phi) const {
  CH_TIME("EBGradient::computeLevelGradient");

  CH_assert(a_gradient.nComp() == SpaceDim);
  CH_assert(a_phi.     nComp() == 1       );

  // TLDR: This routine computes the level gradient, i.e. using finite difference stencils isolated to this level.
  
  const DisjointBoxLayout& dbl   = m_eblg.getDBL();
  const EBISLayout&        ebisl = m_eblg.getEBISL();

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    EBCellFAB&       grad    = a_gradient[dit()];
    const EBCellFAB& phi     = a_phi     [dit()];
    const EBISBox&   ebisBox = ebisl     [dit()];        
    const Box        cellBox = dbl       [dit()];

    const bool isAllCovered = ebisBox.isAllCovered();

    if(!isAllCovered){
      // Compute regular cells
      BaseFab<Real>&       gradFAB = grad.getSingleValuedFAB();
      const BaseFab<Real>& phiFAB  = phi. getSingleValuedFAB();
      FORT_GRADIENT(CHF_FRA(gradFAB),
		    CHF_CONST_FRA1(phiFAB, m_comp),
		    CHF_CONST_REAL(m_dx),
		    CHF_BOX(cellBox));

      // Now do the boundary stencils.
      const BaseIVFAB<VoFStencil>& bndryStencils = m_levelStencils[dit()];

      VoFIterator& vofit = m_levelIterator[dit()];
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex&   vof  = vofit();
	const VoFStencil& sten = bndryStencils(vof, m_comp);

	for (int dir = 0; dir < SpaceDim; dir++){
	  grad(vof, dir) = 0.0;
	}

	for (int i = 0; i < sten.size(); i++){
	  const VolIndex& ivof    = sten.vof(i);
	  const Real&     iweight = sten.weight(i);
	  const int&      ivar    = sten.variable(i); // Note: For the gradient stencil the component is the direction. 

	  grad(vof, ivar) += phi(ivof, m_comp)*iweight;
	}
      }
    }
    else{
      grad.setVal(0.0);
    }

    // Covered data is bogus. 
    for (int dir = 0; dir < SpaceDim; dir++){
      a_gradient[dit()].setCoveredCellVal(0.0, dir);
    }    
  }
}

void EBGradient::computeNormalDerivative(LevelData<EBFluxFAB>& a_gradient,
					 const LevelData<EBCellFAB>& a_phi) const {
  CH_TIME("EBGradient::computeNormalDerivative");

  CH_assert(a_gradient.nComp() == SpaceDim);
  CH_assert(a_phi.     nComp() == 1       );


  // TLDR: This routine computes the level gradient, i.e. using finite difference stencils isolated to this level. It only does this for interior
  // faces, and it only computes the component of the gradient that is normal to the face. 
  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const EBISLayout&        ebisl  = m_eblg.getEBISL();
  const ProblemDomain&     domain = m_eblg.getDomain();

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const EBCellFAB& phi     = a_phi[dit()];
    const EBISBox&   ebisBox = ebisl[dit()];
    const Box        cellBox = dbl  [dit()];

    const EBGraph&   ebgraph = ebisBox.getEBGraph();

    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& grad = a_gradient[dit()][dir];
      
      BaseFab<Real>&       regGradient = grad.getSingleValuedFAB();
      const BaseFab<Real>& regPhi      = phi. getSingleValuedFAB();

      // Do interior faces. This computes the gradient using a standard two-point stencil reaching across the face. 
      Box interiorFaces = cellBox;
      interiorFaces.grow(dir, 1);
      interiorFaces &= domain;
      interiorFaces.grow(dir, -1);
      interiorFaces.surroundingNodes(dir);

      FORT_FACENORMALDERIVATIVE(CHF_FRA1(regGradient, dir),
				CHF_CONST_FRA1(regPhi, m_comp),
				CHF_CONST_INT(dir),
				CHF_CONST_REAL(m_dx),
				CHF_BOX(interiorFaces));

      // Do the same for all interior faces.
      Box interiorCells = cellBox;
      interiorCells.grow(dir, 1);
      interiorCells &= domain;
      interiorCells.grow(dir, -1);

      const IntVectSet irregIVS = ebisBox.getIrregIVS(interiorCells);
      
      for (FaceIterator faceIt(irregIVS, ebgraph, dir, FaceStop::SurroundingNoBoundary); faceIt.ok(); ++faceIt){
	const FaceIndex& face = faceIt();

	const VolIndex& vofHi = face.getVoF(Side::Hi);
	const VolIndex& vofLo = face.getVoF(Side::Lo);	

	grad(face, dir) = (phi(vofHi, m_comp) - phi(vofLo, m_comp))/m_dx;
      }
    }
  }
}

void EBGradient::computeAMRGradient(LevelData<EBCellFAB>&       a_gradient,
				    const LevelData<EBCellFAB>& a_phi,
				    const LevelData<EBCellFAB>& a_phiFine) const {
  CH_TIME("EBGradient::computeAMRGradient(no finer)");

  // TLDR: This routine computes the two-level gradient. It first computes the level gradient and then iterates through
  //       all cells that require a "two-level" view of the gradient. The two-level stencils are applied to buffer data
  //       holders. Once the gradient has been computed (in the buffer data holder), it is copied into a_gradient. 

  this->computeLevelGradient(a_gradient, a_phi);

  // Do corrections near the EBCF. 
  if(m_hasFine && m_hasEBCF){

    // Copy input data to buffers. 
    const Interval interv(m_comp, m_comp);    
    a_phi.    copyTo(interv, m_bufferCoar, interv);
    a_phiFine.copyTo(interv, m_bufferFine, interv);

    // Go through data buffers and compute the gradient (on the buffer grids). 
    for (DataIterator dit(m_bufferDblCoar); dit.ok(); ++dit){
      VoFIterator&     vofit    = m_bufferIterator[dit()];
      BaseIVFAB<Real>& gradCoar = m_bufferGradient[dit()];
      const EBCellFAB& phiCoar  = m_bufferCoar    [dit()];
      const EBCellFAB& phiFine  = m_bufferFine    [dit()];

      // Go through all cells that have a two-level stencil. 
      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();

	const VoFStencil& coarSten = m_bufferStencilsCoar[dit()](vof, m_comp);
	const VoFStencil& fineSten = m_bufferStencilsFine[dit()](vof, m_comp);	

	for (int dir = 0; dir < SpaceDim; dir++){
	  gradCoar(vof, dir) = 0.0;
	}

	// Apply the coarse stencil. 
	for (int i = 0; i < coarSten.size(); i++){
	  const VolIndex& ivof    = coarSten.vof(i);
	  const Real&     iweight = coarSten.weight(i);
	  const int&      ivar    = coarSten.variable(i);

	  gradCoar(vof, ivar) += iweight * phiCoar(ivof, m_comp);
	}

	// Apply the fine stencil. 
	for (int i = 0; i < fineSten.size(); i++){
	  const VolIndex& ivof    = fineSten.vof(i);
	  const Real&     iweight = fineSten.weight(i);
	  const int&      ivar    = fineSten.variable(i);

	  gradCoar(vof, ivar) += iweight * phiFine(ivof, m_comp);
	}	
      }
    }

    // Now copy the result from the buffer grids to the "real" grids. Here, m_bufferGradient is a BaseIVFAB<Real> (for memory reasons)
    // so we can't copy directly to a_gradient. After that we just iterate through those cells and copy the result to a_gradient. 
    m_bufferGradient.copyTo(m_ebcfGradient);

    for (DataIterator dit(m_eblg.getDBL()); dit.ok(); ++dit){
      EBCellFAB&             gradient     = a_gradient    [dit()];
      const BaseIVFAB<Real>& ebcfGradient = m_ebcfGradient[dit()];

      VoFIterator& vofit = m_ebcfIterator[dit()];

      for (vofit.reset(); vofit.ok(); ++vofit){
	const VolIndex& vof = vofit();

	for (int dir = 0; dir < SpaceDim; dir++){
	  gradient(vof, dir) = ebcfGradient(vof, dir);
	}
      }
    }
  }
}

void EBGradient::defineLevelStencils(){
  CH_TIME("EBGradient::defineLevelStencils");

  // TLDR: This just defines finite-difference stencils on this level. We need explicit stencils at the domain boundary AND at
  //       the EB. There are other cases, too, but those are handled by defineStencilsEBCF. 

  const DisjointBoxLayout& dbl    = m_eblg.getDBL();
  const EBISLayout&        ebisl  = m_eblg.getEBISL();
  const ProblemDomain&     domain = m_eblg.getDomain();

  m_levelStencils.define(dbl);
  m_levelIterator.define(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box      cellBox  = dbl  [dit()];
    const EBISBox& ebisbox  = ebisl[dit()];
    const EBGraph& ebgraph  = ebisbox.getEBGraph();


    // Define the stencils for computing the 
    IntVectSet bndryIVS = ebisbox.getIrregIVS(cellBox);
    for (int dir = 0; dir < SpaceDim; dir++){
      Box loBox;
      Box hiBox;
      int hasLo;
      int hasHi;
	
      EBArith::loHi(loBox, hasLo, hiBox, hasHi, domain, cellBox, dir);
	
      if(hasLo){
	bndryIVS |= IntVectSet(loBox);
      }
      if(hasHi){
	bndryIVS |= IntVectSet(hiBox);
      }
    }

    VoFIterator&           bndryIterator = m_levelIterator[dit()];        
    BaseIVFAB<VoFStencil>& bndryStencils = m_levelStencils[dit()];

    bndryIterator.define(bndryIVS, ebgraph);
    bndryStencils.define(bndryIVS, ebgraph, m_nComp);

    for (bndryIterator.reset(); bndryIterator.ok(); ++bndryIterator){
      const VolIndex& vof = bndryIterator();
      VoFStencil& stencil = bndryStencils(vof, m_comp);
      
      stencil.clear();

      // Get derivative in each direction. Note that EBArith might decide to drop order. 
      for (int dir = 0; dir < SpaceDim; dir++){
	VoFStencil derivDirStencil;
	EBArith::getFirstDerivStencilWidthOne(derivDirStencil, vof, ebisbox, dir, m_dx, nullptr, dir);
	stencil += derivDirStencil;	
      }
    }
  }
}

void EBGradient::defineMasks(BoxLayoutData<FArrayBox>& a_coarMaskCF, BoxLayoutData<FArrayBox>& a_coarMaskInvalid){
  CH_TIME("EBGradient::defineMasks");

  CH_assert(m_hasFine);

  // TLDR: This routine defines the masks in a_coarMaskCF and a_coarMaskInvalid. These masks live on the coarse level
  //       and they hold a value of 1 in coarse cells that abut the refinement boundary (a_coarMaskCF) and a value of 1
  //       in cells that are covered by a finer level (a_coarMaskInvalid). These are just constructed by creating some
  //       buffer storage on the (coarsened) fine grid with a value of 1 in the appropriate regions, and then copying
  //       the result to the coarse grid. 

  constexpr Real     zero   = 0.0;
  constexpr Real     one    = 1.0;    
  const     Interval interv = Interval(m_comp, m_comp);

  // Handle to computational grids. 
  const DisjointBoxLayout& dbl        = m_eblg.    getDBL();
  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
  
  const EBISLayout&        ebisl      = m_eblg.    getEBISL();
  const EBISLayout&        ebislFine  = m_eblgFine.getEBISL();
  
  const ProblemDomain&     domain     = m_eblg.    getDomain();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();

  // Create some temporary storage. We coarsen the fine grid and create a BoxLayoutData<FArrayBox> object where each
  // Box is grown by one. We do this because we want to put some data in the ghost cells outside the valid region
  // of the fine level. 
  DisjointBoxLayout dblCoFi;
  coarsen(dblCoFi, dblFine, m_refRat);

  Vector<int> coFiRanks      = dblCoFi.procIDs();  
  Vector<Box> coFiBoxes      = dblCoFi.boxArray();
  Vector<Box> grownCoFiBoxes = dblCoFi.boxArray();
  for(int i = 0; i < coFiBoxes.size(); i++){
    grownCoFiBoxes[i].grow(1);
  }

  BoxLayout coFiLayout     (     coFiBoxes, coFiRanks);
  BoxLayout grownCoFiLayout(grownCoFiBoxes, coFiRanks);
  
  BoxLayoutData<FArrayBox> coFiMaskCF     (grownCoFiLayout, m_nComp); 
  BoxLayoutData<FArrayBox> coFiMaskInvalid(     coFiLayout, m_nComp); 

  // Build the mask regions. 
  for (DataIterator dit(coFiLayout); dit.ok(); ++dit){
    const Box cellBox  =      coFiLayout[dit()];
    const Box grownBox = grownCoFiLayout[dit()];

    // Set the "coarse-fine" mask. It is set to 1 in cells that are not valid cells. We do this by setting the FArrayBox data to one
    // in the entire patch, and then we iterate through the parts of the patch that overlap with valid grid boxes. That part of the
    // data is set to zero. We end up with a bunch of grid patches which hold a value of 1 in ghost cells that do not overlap with
    // the valid region of the grid. 
    coFiMaskCF[dit()].setVal(one);
    coFiMaskCF[dit()].setVal(zero, cellBox, m_comp, m_nComp);
    for (LayoutIterator lit = dblCoFi.layoutIterator(); lit.ok(); ++lit){
      const Box neighborBox = dblCoFi[lit()];
      const Box overlapBox  = neighborBox & grownBox;

      if(!overlapBox.isEmpty()){
	coFiMaskCF[dit()].setVal(zero, overlapBox, m_comp, m_nComp);
      }
    }

    // Set the invalid cell mask on the coarsened fine grids. We set it to 1 everywhere on the fine grid and 0 elsewhere. When we add this to the coarse
    // level we will add 1 into every cell that is covered by a finer level. 
    coFiMaskInvalid[dit()].setVal(one, cellBox, m_comp, m_nComp);
  }

  // Define the input masks. After this, coarMaskCF will have a value of 1 in all cells that abut the fine level. Likewise, coarMaskInvalid
  // will have a value of 1 in all cells that are covered by a finer level.
  Vector<int>      coarRanks = dbl.procIDs();
  Vector<Box>      coarBoxes = dbl.boxArray();
  Vector<Box> grownCoarBoxes = dbl.boxArray();

  for (int i = 0; i < grownCoarBoxes.size(); i++){
    grownCoarBoxes[i].grow(1);
  }

  BoxLayout      coarLayout(     coarBoxes, coarRanks);
  BoxLayout grownCoarLayout(grownCoarBoxes, coarRanks);
  
  a_coarMaskCF     .define(     coarLayout, m_nComp); 
  a_coarMaskInvalid.define(grownCoarLayout, m_nComp);

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    a_coarMaskCF     [dit()].setVal(zero);
    a_coarMaskInvalid[dit()].setVal(zero);
  }

  coFiMaskCF.     addTo(interv, a_coarMaskCF,      interv, ebisl.getDomain()); // Contains > 0 in coarse cells that abut the refinement bounary. 
  coFiMaskInvalid.addTo(interv, a_coarMaskInvalid, interv, ebisl.getDomain()); // Contains > 0 in coarse cells that are covered by a finer level.
}

bool EBGradient::isFiniteDifferenceStencilValid(const IntVect&   a_ivCoar,
						const EBISBox&   a_ebisBox,
						const FArrayBox& a_invalidRegion){
  CH_TIME("EBGradient::isFiniteDifferenceStencilValid");

  // TLDR: Routine which checks if a finite difference stencil is "valid" in the sense that
  //       it does not reach into cut-cells in invalid regions. 

  constexpr Real zero = 0.0;

  bool validStencil = true;

  // Get all VoFs in the input cell
  Vector<VolIndex> coarVoFs = a_ebisBox.getVoFs(a_ivCoar);

  for (const auto& coarVoF : coarVoFs.stdVector()){

    // Build a finite difference stencil. 
    VoFStencil gradSten;
    for (int dir = 0; dir < SpaceDim; dir++){
      VoFStencil derivSten;	      
      EBArith::getFirstDerivStencilWidthOne(derivSten, coarVoF, a_ebisBox, dir, 1.0, nullptr, dir);
      gradSten += derivSten;	      
    }

    // Check if the stencil is OK, i.e. that it does not reach into cut-cells. 
    for (int i = 0; i < gradSten.size(); i++){
      const VolIndex& ivof = gradSten.vof(i);
      const IntVect   iv   = ivof.gridIndex();
      
      const bool isCoveredByFinerCell = a_invalidRegion(iv, m_comp) > zero;
      const bool isIrregularCell      = a_ebisBox.isIrregular(iv);
      
      if(isCoveredByFinerCell && isIrregularCell){
	validStencil = false;
      }
    }
  }

  return validStencil;
}

void EBGradient::defineIteratorsEBCF(const BoxLayoutData<FArrayBox>& a_coarMaskCF, const BoxLayoutData<FArrayBox>& a_coarMaskInvalid){
  CH_TIME("EBGradient::defineIteratorsEBCF");

  // TLDR: This defines which cells need explicit two-level stencils. We move along the refinement boundary on the coarse
  //       level and check if we can use finite-differencing. If not, we flag the cell and the box. We then generate the
  //       buffer grids, which consist of all boxes that contain at least one cell that requires a least squares stencil. 
  
  constexpr Real zero   = 0.0;
  constexpr Real one    = 1.0;      

  const DisjointBoxLayout& dbl        = m_eblg.    getDBL();
  const DisjointBoxLayout& dblFine    = m_eblgFine.getDBL();
  
  const EBISLayout&        ebisl      = m_eblg.    getEBISL();
  const EBISLayout&        ebislFine  = m_eblgFine.getEBISL();
  
  const ProblemDomain&     domain     = m_eblg.    getDomain();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();

  // Define iterators and stencils for EBCF
  m_ebcfIterator.define(dbl);

  Vector<int> ebcfRanks;
  Vector<Box> ebcfBoxes;

  LayoutData<IntVectSet> sets(dbl);

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box      cellBox = dbl  [dit()];
    const EBISBox& ebisBox = ebisl[dit()];
    const EBGraph& ebgraph = ebisBox.getEBGraph();

    const bool isAllRegular = ebisBox.isAllRegular();
    const bool isAllCovered = ebisBox.isAllCovered();
    const bool isIrregular  = !isAllRegular && !isAllCovered;

    const FArrayBox& coarseFineRegion = a_coarMaskCF     [dit()];
    const FArrayBox& invalidRegion    = a_coarMaskInvalid[dit()];

    bool addedBox = false;

    // Determine cells where we need to drop order.
    IntVectSet& ebcfIVS = sets[dit()];
    if(isIrregular){

      // Iterate through the coarse-fine region and check if the finite difference stencil reaches into a cut-cell.
      for (BoxIterator bit(cellBox); bit.ok(); ++bit){
	const IntVect ivCoar = bit();

	if(coarseFineRegion(ivCoar, m_comp) > zero){
	  const bool hasStencil = this->isFiniteDifferenceStencilValid(ivCoar, ebisBox, invalidRegion);

	  // Ok, found a stencil that won't work. 
	  if(!hasStencil){
	    ebcfIVS |= ivCoar;

	    if(addedBox == false){ // One of the cells needs a better stencil so we allocate data in this box. 
	      ebcfBoxes.push_back(cellBox);
	      addedBox = true; 
	    }
	  }
	}
      }

    }

    // Define the iterator.
    m_ebcfIterator[dit()].define(ebcfIVS, ebgraph);
  }

  // Define storage for holding gradient
  m_ebcfGradient.define(dbl, SpaceDim, IntVect::Zero, BaseIVFactory<Real>(ebisl, sets));

  // Make the grids that contain irregular boxes that require better stencils. We will use this as a buffered layout.
  LoadBalancing::gatherBoxes(ebcfBoxes);
  LoadBalancing::makeBalance(ebcfRanks, ebcfBoxes);

  // Fill the buffer layouts. Note that we need one ghost cell on the coarse level but m_refRat ghost cells on the fine level. 
  if(ebcfBoxes.size() > 0){
    m_hasEBCF = true;
    const EBIndexSpace* const ebisPtr = m_eblg.getEBIS();
    m_bufferDblCoar.define(ebcfBoxes, ebcfRanks);
    refine(m_bufferDblFine, m_bufferDblCoar, m_refRat);

    ebisPtr->fillEBISLayout(m_bufferEBISLCoar, m_bufferDblCoar, domain,  1);
    ebisPtr->fillEBISLayout(m_bufferEBISLFine, m_bufferDblFine, domainFine, m_refRat);

    m_bufferEBISLCoar.setMaxRefinementRatio(m_refRat, ebisPtr);
  }
  else{
    m_hasEBCF = false;
  }
}

void EBGradient::defineBuffers(BoxLayoutData<FArrayBox>&       a_bufferCoarMaskCF,
			       BoxLayoutData<FArrayBox>&       a_bufferCoarMaskInvalid,
			       const BoxLayoutData<FArrayBox>& a_coarMaskCF,
			       const BoxLayoutData<FArrayBox>& a_coarMaskInvalid){
  CH_TIME("EBGradient::defineBuffers");

  // TLDR: This routine defines the required buffers for holding data for dual level stencils. We use the input masks for determining
  //       which cells need a least squares stencil. Much of this code simply consists of creating masks that are viewable from the
  //       buffer grids (rather than m_eblg). 

  constexpr Real zero   = 0.0;
  constexpr Real one    = 1.0;
  
  const Interval interv = Interval(m_comp, m_comp);

  const DisjointBoxLayout& dbl        = m_bufferDblCoar;
  const DisjointBoxLayout& dblFine    = m_bufferDblFine;
  
  const EBISLayout&        ebisl      = m_bufferEBISLCoar;
  const EBISLayout&        ebislFine  = m_bufferEBISLFine;
  
  const ProblemDomain&     domain     = m_eblg.    getDomain();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();

  // We need buffers on this layout
  Vector<int>      coarRanks = dbl.procIDs();
  Vector<Box>      coarBoxes = dbl.boxArray();
  Vector<Box> grownFineBoxes = dbl.boxArray();  
  Vector<Box> grownCoarBoxes = dbl.boxArray();
  for (auto& bx : grownCoarBoxes.stdVector()){
    bx.grow(1);
  }

  for (auto& bx : grownFineBoxes.stdVector()){
    bx.grow(1);
    bx.refine(m_refRat);
  }  

  BoxLayout      coarLayout(     coarBoxes, coarRanks);
  BoxLayout grownCoarLayout(grownCoarBoxes, coarRanks);
  BoxLayout grownFineLayout(grownFineBoxes, coarRanks);  

  a_bufferCoarMaskCF.     define(     coarLayout, m_nComp);
  a_bufferCoarMaskInvalid.define(grownCoarLayout, m_nComp);

  m_bufferCoar.define(grownCoarLayout, m_nComp, EBCellFactory(ebisl    ));
  m_bufferFine.define(grownFineLayout, m_nComp, EBCellFactory(ebislFine));    

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    a_bufferCoarMaskCF     [dit()].setVal(zero);
    a_bufferCoarMaskInvalid[dit()].setVal(zero);

    m_bufferCoar[dit()].setVal(zero);
    m_bufferFine[dit()].setVal(zero);
  }

  // Transfer the masks to the reduced grids.
  a_coarMaskCF.     addTo(interv, a_bufferCoarMaskCF,      interv, domain);
  a_coarMaskInvalid.addTo(interv, a_bufferCoarMaskInvalid, interv, domain);

  // Define the necessary data on the buffer grids
  m_bufferIterator.    define(dbl);  
  m_bufferStencilsFine.define(dbl);  
  m_bufferStencilsCoar.define(dbl);

  // Go through the "buffer" grids and define all cells where we need better stencils.
  LayoutData<IntVectSet> sets(dbl);
  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box      cellBox = dbl  [dit()];
    const EBISBox& ebisBox = ebisl[dit()];
    const EBGraph& ebgraph = ebisBox.getEBGraph();

    IntVectSet& ebcfIVS = sets[dit()];

    const FArrayBox& coarMaskCF      = a_bufferCoarMaskCF     [dit()];
    const FArrayBox& coarMaskInvalid = a_bufferCoarMaskInvalid[dit()];    

    for (BoxIterator bit(cellBox); bit.ok(); ++bit){
      const IntVect ivCoar = bit();

      if(coarMaskCF(ivCoar, m_comp) > zero){ // Cell lies on the CF bounary
    	const bool validStencil = this->isFiniteDifferenceStencilValid(ivCoar, ebisBox, coarMaskInvalid);

    	if(!validStencil){
    	  ebcfIVS |= ivCoar;
    	}
      }
    }

    m_bufferIterator    [dit()].define(ebcfIVS, ebgraph          );
    m_bufferStencilsFine[dit()].define(ebcfIVS, ebgraph, m_nComp );
    m_bufferStencilsCoar[dit()].define(ebcfIVS, ebgraph, m_nComp );
  }

  // Define storage for gradient
  m_bufferGradient.define(dbl, SpaceDim, IntVect::Zero, BaseIVFactory<Real>(ebisl, sets));
}

void EBGradient::defineStencilsEBCF(const BoxLayoutData<FArrayBox>& a_bufferCoarMaskInvalid){
  CH_TIME("EBGradient::defineStencilsEBCF");

  // TLDR: This code iterates through the buffer grids and defines least squares stencils in the
  //       various cells that need it. The input mask is used to determine which cells we can use
  //       on the fine/coarse levels. 

  const DisjointBoxLayout& dbl        = m_bufferDblCoar;
  const DisjointBoxLayout& dblFine    = m_bufferDblFine;
  
  const EBISLayout&        ebisl      = m_bufferEBISLCoar;
  const EBISLayout&        ebislFine  = m_bufferEBISLFine;
  
  const ProblemDomain&     domain     = m_eblg.    getDomain();
  const ProblemDomain&     domainFine = m_eblgFine.getDomain();

  for (DataIterator dit(dbl); dit.ok(); ++dit){
    const Box      cellBox     = dbl      [dit()];
    const Box      cellBoxFine = dblFine  [dit()];
    
    const EBISBox& ebisBox     = ebisl    [dit()];
    const EBISBox& ebisBoxFine = ebislFine[dit()];
    
    const EBGraph& ebgraph     = ebisBox.    getEBGraph();
    const EBGraph& ebgraphFine = ebisBoxFine.getEBGraph();

    const Box grownCellBox     = grow(cellBox, 1);
    const Box grownCellBoxFine = refine(grownCellBox, m_refRat);

    const FArrayBox& invalidRegion = a_bufferCoarMaskInvalid[dit()];

    // Make the invalid region mask into DenseIntVectSets because that is what
    // LeastSquares mostly want. 
    DenseIntVectSet validRegionCoar(grownCellBox,     true);
    DenseIntVectSet validRegionFine(grownCellBoxFine, false);

    for (BoxIterator bit(grownCellBox); bit.ok(); ++bit){
      const IntVect ivCoar = bit();

      if(invalidRegion(ivCoar, m_comp)){
	validRegionCoar -= ivCoar;

	Box bx(ivCoar, ivCoar);
	bx.refine(m_refRat);
	validRegionFine |= bx;
      }
    }

    // Define the iterator and stencils for places where we need to drop order.
    VoFIterator&           vofit         = m_bufferIterator    [dit()];
    BaseIVFAB<VoFStencil>& coarStencils  = m_bufferStencilsCoar[dit()];
    BaseIVFAB<VoFStencil>& fineStencils  = m_bufferStencilsFine[dit()];

    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      VoFStencil& coarStencil = coarStencils(vof, m_comp);
      VoFStencil& fineStencil = fineStencils(vof, m_comp);
      
      coarStencil.clear();
      fineStencil.clear();

      int  order        = m_order;
      bool foundStencil = false;

      while(order > 0 && !foundStencil){
	foundStencil = this->getLeastSquaresStencil(coarStencil,
						    fineStencil,
						    vof,
						    m_dataLocation,
						    ebisl,
						    ebislFine,
						    dit(),
						    validRegionCoar,
						    validRegionFine,
						    m_dx,
						    m_dxFine,
						    order,
						    m_weighting);
	order--;
      }

      // As a last effort, get a finite difference stencil after all. 
      if(!foundStencil){
	coarStencil.clear();
	fineStencil.clear();

	this->getFiniteDifferenceStencil(coarStencil, vof, ebisBox, DenseIntVectSet(cellBox, true), m_dx);

	MayDay::Warning("CD_EBGradient::defineStencilsEBCF -- could not find stencil!");
      }
    }
  }
}

bool EBGradient::getFiniteDifferenceStencil(VoFStencil&            a_stencil,
					    const VolIndex&        a_vof,
					    const EBISBox&         a_ebisBox,
					    const DenseIntVectSet& a_validRegion,
					    const Real             a_dx){
  CH_TIME("EBGradient::getFiniteDifferenceStencil");

  // TLDR: This routine computes a finite difference stencil. If it does not have any valid cells to use (defined by a_validRegion)
  //       it will not return a stencil. 

  a_stencil.clear();

  bool foundStencil = false;

  const IntVect iv = a_vof.gridIndex();
  
  for (int dir = 0; dir < SpaceDim; dir++){
    VoFStencil derivStencil;
	
    // Check for cells to the low/high side of this one.
    const IntVect ivLo = iv - BASISV(dir);
    const IntVect ivHi = iv + BASISV(dir);

    const bool isLoCellCovered    = a_validRegion[ivLo];
    const bool isHiCellCovered    = a_validRegion[ivHi];

    const bool isLoCellIrregular  = a_ebisBox.isIrregular(ivLo);
    const bool isHiCellIrregular  = a_ebisBox.isIrregular(ivHi);

    const Vector<VolIndex> loVoFs = a_ebisBox.getVoFs(a_vof, dir, Side::Lo, 1);
    const Vector<VolIndex> hiVoFs = a_ebisBox.getVoFs(a_vof, dir, Side::Hi, 1);

    const int numLoVoFs = loVoFs.size();
    const int numHiVoFs = hiVoFs.size();

    const bool useLoCell = (numLoVoFs > 0) && !(isLoCellIrregular && isLoCellCovered);
    const bool useHiCell = (numHiVoFs > 0) && !(isHiCellIrregular && isHiCellCovered);			

    // Use centered differencing in dir-direction if we can. Otherwise, drop order to forward/backward differences. 
    if(useLoCell && useHiCell){
      for (int i = 0; i < numHiVoFs; i++){
	derivStencil.add(hiVoFs[i], 1./numHiVoFs, dir);
      }

      for (int i = 0; i < numLoVoFs; i++){
	derivStencil.add(loVoFs[i], -1./numLoVoFs, dir);
      }

      derivStencil *= 1./(2.0*a_dx);
    }
    else if(useLoCell){ 
      derivStencil.add(a_vof, 1.0, dir);

      for (int i = 0; i < numLoVoFs; i++){
	derivStencil.add(loVoFs[i], -1./numLoVoFs, dir);
      }

      derivStencil *= 1./a_dx;
    }
    else if(useHiCell){ 
      derivStencil.add(a_vof, -1.0, dir);

      for (int i = 0; i < numHiVoFs; i++){
	derivStencil.add(hiVoFs[i], 1./numHiVoFs, dir);
      }

      derivStencil *= 1./a_dx;
    }

    if(useLoCell || useHiCell){
      foundStencil = true;
      
      a_stencil += derivStencil;
    }
    else{
      foundStencil = false;
    }
  }

  return foundStencil;
}

bool EBGradient::getLeastSquaresStencil(VoFStencil&            a_stencilCoar,
					VoFStencil&            a_stencilFine,
					const VolIndex&        a_vofCoar,
					const CellLocation&    a_dataLocation,
					const EBISLayout&      a_ebislCoar,
					const EBISLayout&      a_ebislFine,
					const DataIndex&       a_dit,					
					const DenseIntVectSet& a_validCellsCoar,
					const DenseIntVectSet& a_validCellsFine,					
					const Real&            a_dxCoar,
					const Real&            a_dxFine,					
					const int&             a_order,
					const int&             a_weight) {
  CH_TIME("EBGradient::getLeastSquaresStencil");

  // TLDR: This routine computes a two-level least squares stencil. We include cells that are within a range of 1 on the coarse level, but only
  //       those in a monotone path from the input vof. The fine cells are defined as the cells that are available through a coarsening of the
  //       coarse cells. Since the resulting system might sometimes be too big, we trim the system down to a specified size. Once the least squares
  //       system has been "solved", i.e. a minimization of the truncation order for the various expansions has been achieved, we formulate the stencil
  //       using the output of LeastSquares::computeDualLevelStencils. Note that the solution on the input vof is known, so it is pruned from the
  //       system of equations. 

  bool foundStencil;

  a_stencilCoar.clear();
  a_stencilFine.clear();

  const EBISBox& ebisBoxCoar = a_ebislCoar[a_dit];
  const EBISBox& ebisBoxFine = a_ebislFine[a_dit];  
    
  // Get all coarse cells within radius 1.
  const int coarRadius = 1;
  const int fineRadius = m_refRat;

  // Get coar vofs in range. The fine VoFs are defined as the VoFs that are available through refinement of the coarse vofs. 
  Vector<VolIndex> coarVoFs = VofUtils::getVofsInRadius(a_vofCoar, ebisBoxCoar, coarRadius, VofUtils::Connectivity::MonotonePath, false);
  Vector<VolIndex> fineVoFs;
  for (int i = 0; i < coarVoFs.size(); i++){
    fineVoFs.append(a_ebislCoar.refine(coarVoFs[i], m_refRat, a_dit));
  }

  // Only unique, in case something went wrong. 
  //  VofUtils::onlyUnique(coarVoFs);
  VofUtils::onlyUnique(fineVoFs);  

  // Only include fine VoFs that lie with a_validCellsFine.
  VofUtils::includeCells(coarVoFs, a_validCellsCoar);  
  VofUtils::includeCells(fineVoFs, a_validCellsFine);

  // See if we have enough cells to solve the system of equations. The "-1" is because we interpolate to the cell center/centroid, but we already know
  // phi at this point. 
  const int numEquations = fineVoFs.size() + coarVoFs.size();
  const int numUnknowns  = LeastSquares::getTaylorExpansionSize(a_order) - 1;

  if(numEquations > numUnknowns){

    // In many cases we will have WAY too many equations for the specified order. This is particularly true in 3D
    // because the number of coar vofs included in a radius r from the ghost vof can be 3^3 = 27. In addition, if we
    // use a refinement factor of 4 and just half of those cells are covered by finer cells, we end up with a system
    // size of = 13 + 14*4^3 = 97. That's way more than we need for order 2 which requires 10 equations. Since the SVD
    // decomposition scales as O(n^3), we trim the system size to bring the cost down, discarding the cells that are
    // furthest away. 
    std::vector<VolIndex>& fineVofsTrimmedSize = fineVoFs.stdVector();
    std::vector<VolIndex>& coarVofsTrimmedSize = coarVoFs.stdVector();

    // Coordinates of the vof that we will interpolate to (excluding lower-left corner because of the subtraction
    // in the comparators). 
    const RealVect x0 = Location::position(a_dataLocation, a_vofCoar, ebisBoxCoar, a_dxCoar);

    // For sorting fine vofs, based on distance to the ghost vof. Shortest distance goes first. 
    auto comparatorFine = [&loc     = a_dataLocation,
			   &p       = x0,
			   &ebisbox = ebisBoxFine,
			   &dx      = a_dxFine](const VolIndex& v1, const VolIndex& v2) -> bool {
      const RealVect d1 = Location::position(loc, v1, ebisbox, dx) - p;
      const RealVect d2 = Location::position(loc, v2, ebisbox, dx) - p;

      const Real l1 = d1.vectorLength();
      const Real l2 = d2.vectorLength();
      
      return l1 < l2;
    };

    // For sorting coar vofs, based on distance to the ghost vof. Shortest distance goes first. 
    auto comparatorCoar = [&loc     = a_dataLocation,
			   &p       = x0,
			   &ebisbox = ebisBoxCoar,
			   &dx      = a_dxCoar](const VolIndex& v1, const VolIndex& v2) -> bool {
      const RealVect d1 = Location::position(loc, v1, ebisbox, dx) - p;
      const RealVect d2 = Location::position(loc, v2, ebisbox, dx) - p;

      const Real l1 = d1.vectorLength();
      const Real l2 = d2.vectorLength();
      
      return l1 < l2;
    };    

    // Sort and trim the system size. 
    std::sort(fineVofsTrimmedSize.begin(), fineVofsTrimmedSize.end(), comparatorFine);
    std::sort(coarVofsTrimmedSize.begin(), coarVofsTrimmedSize.end(), comparatorCoar);

    const int curFineSize = fineVofsTrimmedSize.size();
    const int curCoarSize = coarVofsTrimmedSize.size();
    
    fineVofsTrimmedSize.resize(std::min(2*numUnknowns, curFineSize));
    coarVofsTrimmedSize.resize(std::min(2*numUnknowns, curCoarSize));
    
    // Build the displacement vectors
    Vector<RealVect> fineDisplacements;
    Vector<RealVect> coarDisplacements;
    
    for (const auto& coarVoF : coarVoFs.stdVector()){
      coarDisplacements.push_back(LeastSquares::displacement(a_dataLocation, a_dataLocation, a_vofCoar, coarVoF, ebisBoxCoar, a_dxCoar));
    }

    for (const auto& fineVoF : fineVoFs.stdVector()){
      fineDisplacements.push_back(LeastSquares::displacement(a_dataLocation, a_dataLocation, a_vofCoar, fineVoF, ebisBoxCoar, ebisBoxFine, a_dxCoar, a_dxFine));
    }    

    // Now solve the least squares system for the gradient at the input VoF. After this we have found an expression which minimizes the
    // truncation error in the neighborhood of the VoF. The result is is a two-level stencil consisting of cells on both the fine and
    // the coarse levels.
    // These are the unknown terms that we want stencils for.
    IntVectSet derivTerms;
    IntVectSet knownTerms;

    for (int dir = 0; dir < SpaceDim; dir++){
      derivTerms |= BASISV(dir);
    }
    knownTerms |= IntVect::Zero;
    
    const std::map<IntVect, std::pair<VoFStencil, VoFStencil> > stencils = LeastSquares::computeDualLevelStencils<float>(derivTerms,
															 knownTerms,
															 fineVoFs,
															 coarVoFs,
															 fineDisplacements,
															 coarDisplacements,
															 a_weight,
															 a_order);
    // LeastSquares returns a map over all derivatives (unknowns) in the Taylor series. These are stored
    // as IntVects so that IntVect(1,1) = d^2/(dxdy) and so on. We are after IntVect(1,0) = d/dx, IntVect(0,1) = d/dy. We fetch
    // those and place them in a_stencilFine and a_stencilCoar. We encode the direction in the stencil variable (in a_stencilFine) and
    // a_stencilCoar. 
    for (int dir = 0; dir < SpaceDim; dir++){
      const IntVect deriv = BASISV(dir);

      const VoFStencil& fineDerivStencil = stencils.at(deriv).first;
      const VoFStencil& coarDerivStencil = stencils.at(deriv).second;

      for (int i = 0; i < fineDerivStencil.size(); i++){
	const VolIndex& vof    = fineDerivStencil.vof   (i);
	const Real&     weight = fineDerivStencil.weight(i);

	a_stencilFine.add(vof, weight, dir);
      }

      for (int i = 0; i < coarDerivStencil.size(); i++){
	const VolIndex& vof    = coarDerivStencil.vof   (i);
	const Real&     weight = coarDerivStencil.weight(i);

	a_stencilCoar.add(vof, weight, dir);
      }
    }


    // We have modified the right-hand side of the least squares system by pruning a_vofCoar from the system of equations (it's value is known). So,
    // our least squares solution is actually something like
    //
    //    [df/dx]     [A11 A12 ... ] [f(x0) - f(x)]
    //    [df/dy]  =  [A21 A22 ... ] [f(x1) - f(x)]
    //    [df/dz]     [A31 A32 ... ] [f(x2) - f(x)]
    // 
    // where f(x) = phi(a_vofCoar). LeastSquares::computeGradSten does not case about the modified right-hand side and the stencils it returns
    // only account for f(x0), f(x1) etc. We need to add in the contribution from a_vofCoar, which is fortunately just the sum of the weights for 
    // each derivative.     
    for (int dir = 0; dir < SpaceDim; dir++){
	Real coarVofWeight = 0.0;
	
	for (int i = 0; i < a_stencilFine.size(); i++){
	  const int  curVar    = a_stencilFine.variable(i);
	  const Real curWeight = a_stencilFine.weight(i);
	  
	  if(curVar == dir){
	    coarVofWeight += curWeight;
	  }
	}

	for (int i = 0; i < a_stencilCoar.size(); i++){
	  const int  curVar    = a_stencilCoar.variable(i);
	  const Real curWeight = a_stencilCoar.weight(i);
	  
	  if(curVar == dir){
	    coarVofWeight += curWeight;
	  }
	}	
	
	a_stencilCoar.add(a_vofCoar, -1.0*coarVofWeight, dir);
      }          

    foundStencil = true;    
  }
  else{
    foundStencil = false;
  }

  return foundStencil;
}

#include <CD_NamespaceFooter.H>
