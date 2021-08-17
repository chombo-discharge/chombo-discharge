/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrMuscl.cpp
  @brief  Implementation of CD_CdrMuscl.H
  @author Robert Marskar
  @todo   computeSlopes is missing proper treatment of domain boundaries!
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CdrMuscl.H>
#include <CD_CdrMusclF_F.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

CdrMuscl::CdrMuscl(){
  CH_TIME("CdrMuscl::CdrMuscl()");
  
  // Class and instantiatio name
  m_className   = "CdrMuscl";
  m_name        = "CdrMuscl";

  m_limitSlopes = true;
}

CdrMuscl::~CdrMuscl(){
  CH_TIME("CdrMuscl::~CdrMuscl()");
}

void CdrMuscl::parseOptions(){
  CH_TIME("CdrGodunov::parseOptions()");
  if(m_verbosity > 5){
    pout() << m_name + "::parseOptions()" << endl;
  }
  
  this->parseRngSeed();                // Parses RNG seed
  this->parsePlotMode();               // Parses plot mode
  this->parseDomainBc();               // Parses domain BC options
  this->parseSlopeLimiter();           // Parses slope limiter settings
  this->parsePlotVariables();          // Parses plot variables
  this->parseMultigridSettings();      // Parses multigrid settings. 
  this->parseDivergenceComputation();  // Non-conservative divergence blending
}

void CdrMuscl::parseRuntimeOptions(){
  CH_TIME("CdrMuscl::parseRuntimeOptions()");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRuntimeOptions()" << endl;
  }

  this->parsePlotMode();               // Parses plot mode
  this->parseDomainBc();               // Parses domain BC options
  this->parseSlopeLimiter();           // Parses slope limiter settings
  this->parsePlotVariables();          // Parses plot variables
  this->parseMultigridSettings();      // Parses multigrid settings. 
  this->parseDivergenceComputation();  // Non-conservative divergence blending. 
}

void CdrMuscl::parseSlopeLimiter(){
  CH_TIME("CdrMuscl::parseSlopeLimiter()");
  if(m_verbosity > 5){
    pout() << m_name + "::parseSlopeLimiter()" << endl;
  }
  
  ParmParse pp(m_className.c_str());
  
  pp.get("limit_slopes", m_limitSlopes);
}

void CdrMuscl::advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_cellPhi, const Real a_extrapDt){
  CH_TIME("CdrMuscl::advectToFaces(EBAMRFluxData, EBAMRCellData, Real)");
  if(m_verbosity > 5){
    pout() << m_name + "::advectToFaces(EBAMRFluxData, EBAMRCellData, Real)" << endl;
  }

  CH_assert(a_facePhi[0]->nComp() == 1);
  CH_assert(a_cellPhi[0]->nComp() == 1);  

  DataOps::setValue(a_facePhi, 0.0);

  // Ghost cells need to be interpolated. We make a copy of a_cellPhi which we use for that. 
  EBAMRCellData phi; 
  m_amr->allocate(phi, m_realm, m_phase, m_nComp);
  DataOps::copy(phi, a_cellPhi);

  m_amr->averageDown(phi, m_realm, m_phase);
  m_amr->interpGhost(phi, m_realm, m_phase);

  for (int lvl = 0; lvl <= m_amr->getFinestLevel(); lvl++){
    const DisjointBoxLayout& dbl    = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout&        ebisl  = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const ProblemDomain&     domain = m_amr->getDomains()[lvl];

    for (DataIterator dit(dbl); dit.ok(); ++dit){
      EBFluxFAB&       facePhi = (*a_facePhi     [lvl])[dit()];
      const EBCellFAB& cellPhi = (*phi           [lvl])[dit()];
      const EBFluxFAB& velo    = (*m_faceVelocity[lvl])[dit()];
      
      const Box        cellBox = dbl  [dit()];      
      const EBISBox&   ebisbox = ebisl[dit()];

      // Limit slopes and solve Riemann problem (this is just the upwind state). 
      Box grownBox = cellBox;
      grownBox.grow(1);
      EBCellFAB deltaC(ebisbox, grownBox, SpaceDim); // Cell-centered slopes
      deltaC.setVal(0.0);

      if(m_limitSlopes){
	this->computeSlopes(deltaC, cellPhi, cellBox, domain, lvl, dit());
      }
      
      this->upwind(facePhi, deltaC, cellPhi, velo, domain, cellBox, lvl, dit());
    }
  }
}

void CdrMuscl::computeSlopes(EBCellFAB&           a_deltaC,
			     const EBCellFAB&     a_cellPhi,
			     const Box&           a_cellBox,
			     const ProblemDomain& a_domain,
			     const int            a_level,
			     const DataIndex&     a_dit){
  CH_TIME("CdrMuscl::computeSlopes(EBCellFAB, EBCellFAB, Box, ProblemDomain, int, DataIndex)");
  if(m_verbosity > 5){
    pout() << m_name + "::computeSlopes(EBCellFAB, EBCellFAB, Box, ProblemDomain, int, DataIndex)" << endl;
  }

  CH_assert(a_deltaC. nComp() == SpaceDim);
  CH_assert(a_cellPhi.nComp() == 1);

  const Box& domainBox      = a_domain.domainBox();
  const EBISBox& ebisbox     = a_cellPhi.getEBISBox();

  // Compute slopes in regular cells. 
  BaseFab<Real>&       slopesReg = a_deltaC.getSingleValuedFAB();
  const BaseFab<Real>& phiReg    = a_cellPhi.getSingleValuedFAB();

  // Compute slopes in regular grid cells. 
  for (int dir = 0; dir < SpaceDim; dir++){
    Box compBox = a_cellBox;
    compBox.grow(1);
    FORT_MUSCL_SLOPES(CHF_FRA1(slopesReg, dir),
		      CHF_CONST_FRA1(phiReg, m_comp),
		      CHF_CONST_INT(dir),
		      CHF_BOX(compBox));
  }

  // Same kernel, but in irregular cells. 
  for (int dir = 0; dir < SpaceDim; dir++){

    VoFIterator& vofit = (*m_amr->getVofIterator(m_realm, m_phase)[a_level])[a_dit];
    
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const IntVect iv    = vof.gridIndex();

      const bool onLoSide = (iv[dir] == domainBox.smallEnd(dir));
      const bool onHiSide = (iv[dir] == domainBox.bigEnd  (dir));
      
      const bool hasFacesLeft = (ebisbox.numFaces(vof, dir, Side::Lo) == 1) && !onLoSide;
      const bool hasFacesRigh = (ebisbox.numFaces(vof, dir, Side::Hi) == 1) && !onHiSide;

      VolIndex vofLeft;
      VolIndex vofRigh;
      
      Real dwc   = 0.0;
      Real dwl   = 0.0;
      Real dwr   = 0.0;
      Real dwmin = 0.0;

      Real phiLeft = 0.0;
      Real phiRigh = 0.0;
      Real phiCent = a_cellPhi(vof, m_comp);

      // Compute left and right slope
      if(hasFacesLeft){
	Vector<FaceIndex> facesLeft = ebisbox.getFaces(vof, dir, Side::Lo);
	vofLeft = facesLeft[0].getVoF(Side::Lo);
	phiLeft = a_cellPhi(vofLeft, m_comp);
	dwl     = phiCent - phiLeft;
      }
      if(hasFacesRigh){
	Vector<FaceIndex> facesRigh = ebisbox.getFaces(vof, dir, Side::Hi);
	vofRigh = facesRigh[0].getVoF(Side::Hi);
	phiRigh = a_cellPhi(vofRigh, m_comp);
	dwr     = phiRigh - phiCent;
      }

      if(hasFacesLeft && hasFacesRigh){
	dwc = 0.5*(dwl + dwr);
      }
      else if(!hasFacesLeft && !hasFacesRigh){
	dwl = 0.0;
	dwc = 0.0;
	dwr = 0.0;
      }
      else if(!hasFacesLeft && hasFacesRigh){
	dwl = dwr;
	dwc = dwr;
      }
      else if(hasFacesLeft && !hasFacesRigh){
	dwr = dwl;
	dwc = dwl;
      }
      else{
	MayDay::Error("CdrMuscl::computeSlopes - missed a case");
      }

      // Limit the slopes. 
      dwmin = 2.0*std::min(std::abs(dwl), std::abs(dwr));
      if(dwl*dwr < 0.0){
	dwc = 0.0;
      }
      else{
	if(dwc < 0.0){
	  dwc = -std::min(dwmin, std::abs(dwc));
	}
	else if (dwc > 0.0){
	  dwc = std::min(dwmin, std::abs(dwc));
	}
	else{
	  dwc = 0.0;
	}
      }

      if(std::isnan(dwc)){
	MayDay::Error("CdrMuscl::computeSlopes - dwc != dwc.");
      }

      //      a_deltaC(vof, m_comp) = dwc;
      a_deltaC(vof, m_comp) = 0.0;
    }
  }
}

void CdrMuscl::upwind(EBFluxFAB&           a_facePhi,
		      const EBCellFAB&     a_slopes,
		      const EBCellFAB&     a_cellPhi,
		      const EBFluxFAB&     a_faceVel,
		      const ProblemDomain& a_domain,
		      const Box&           a_cellBox,
		      const int&           a_level,
		      const DataIndex&     a_dit){
  CH_TIME("CdrMuscl::upwind(EBFluxFAB, EBCellFAB, EBCellFAB, EBFluxFAB, ProblemDomain, Box, int, DataIndex)");
  if(m_verbosity > 99){
    pout() << m_name + "::upwind(EBFluxFAB, EBCellFAB, EBCellFAB, EBFluxFAB, ProblemDomain, Box, int, DataIndex)" << endl;
  }

  CH_assert(a_facePhi.nComp() == 1);
  CH_assert(a_slopes. nComp() == SpaceDim);
  CH_assert(a_cellPhi.nComp() == 1);
  CH_assert(a_faceVel.nComp() == 1);  

  for (int dir = 0; dir < SpaceDim; dir++){
    EBFaceFAB& facePhi       = a_facePhi[dir];
    const EBFaceFAB& faceVel = a_faceVel[dir];
    
    BaseFab<Real>&      regFacePhi = a_facePhi[dir].getSingleValuedFAB();
    const BaseFab<Real>& regSlopes = a_slopes.      getSingleValuedFAB();
    const BaseFab<Real>& regStates = a_cellPhi.     getSingleValuedFAB();
    const BaseFab<Real>& regVelo   = a_faceVel[dir].getSingleValuedFAB();
    
    // Fortran kernel. Computes the extrapolation of cell-centered states to the face centers and selects
    // the upwind side. 
    Box facebox = a_cellBox;
    facebox.surroundingNodes();
    regFacePhi.setVal(0.0);
    
    FORT_MUSCL_UPWIND(CHF_FRA1(regFacePhi,   m_comp),
		      CHF_CONST_FRA1(regSlopes, dir),
		      CHF_CONST_FRA1(regStates, m_comp),
		      CHF_CONST_FRA1(regVelo,   m_comp),
		      CHF_CONST_INT(dir),
		      CHF_BOX(facebox));

    // Compute extrapolated primitives in for irregular faces and solve the Riemann problem (i.e., just upwind). Just
    // like the Fortran kernel, but for irregular grid faces.
    const EBISBox& ebisbox = a_cellPhi.getEBISBox();
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    
    const Vector<FaceIndex> irregFaces = ebgraph.getIrregFaces(a_cellBox, dir);
    
    for (int i = 0; i < irregFaces.size(); i++){
      const FaceIndex& face = irregFaces[i];
      
      if(!face.isBoundary()){
	const VolIndex& vofLeft = face.getVoF(Side::Lo);
	const VolIndex& vofRigh = face.getVoF(Side::Hi);

	const Real velo     = faceVel(face, m_comp);
	const Real primLeft = a_cellPhi(vofLeft, m_comp) + 0.5*a_slopes(vofLeft, dir);
	const Real primRigh = a_cellPhi(vofRigh, m_comp) - 0.5*a_slopes(vofRigh, dir);

	if(velo > 0.0){
	  facePhi(face, m_comp) = primLeft;
	}
	else if (velo < 0.0){
	  facePhi(face, m_comp) = primRigh;
	}
	else{
	  facePhi(face, m_comp) = 0.0;
	}
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
