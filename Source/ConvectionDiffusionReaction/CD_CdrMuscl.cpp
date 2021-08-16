/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CdrMuscl.cpp
  @brief  Implementation of CD_CdrMuscl.H
  @author Robert Marskar
*/

// Chombo includes
#include <ParmParse.H>

// Our includes
#include <CD_CdrMuscl.H>
#include <CD_CdrMusclF_F.H>
#include <CD_DataOps.H>
#include <CD_NamespaceHeader.H>

CdrMuscl::CdrMuscl(){
  m_className = "CdrMuscl";
  m_name      = "CdrMuscl";
}

CdrMuscl::~CdrMuscl(){

}

void CdrMuscl::parseOptions(){
  parseRngSeed();                // Parses RNG seed
  parsePlotMode();               // Parses plot mode
  parseDomainBc();               // Parses domain BC options
  parseSlopeLimiter();           // Parses slope limiter settings
  parsePlotVariables();          // Parses plot variables
  parseMultigridSettings();      // Parses solver parameters for geometric multigrid
  parseDivergenceComputation();  // Nonlinear divergence blending
}

void CdrMuscl::parseRuntimeOptions(){
  CH_TIME("CdrMuscl::parseRuntimeOptions");
  if(m_verbosity > 5){
    pout() << m_name + "::parseRuntimeOptions" << endl;
  }

  parsePlotMode();               // Parses plot mode
  parseDomainBc();               // Parses domain BC options
  parseSlopeLimiter();           // Parses slope limiter settings
  parsePlotVariables();          // Parses plot variables
  parseMultigridSettings();      // Parses solver parameters for geometric multigrid
  parseDivergenceComputation();  // Nonlinear divergence blending
}

void CdrMuscl::parseSlopeLimiter(){
  ParmParse pp(m_className.c_str());

  std::string str;
  pp.get("limit_slopes", str);
  m_slopelim = (str == "true") ? true : false;
}

void CdrMuscl::advance_advect(EBAMRCellData& a_phi, const Real a_dt){
  CH_TIME("CdrMuscl::advance_advect");
  if(m_verbosity > 5){
    pout() << m_name + "::advance_advect" << endl;
  }

  MayDay::Abort("CdrMuscl::advance_advect - not implemented (yet)");
}

void CdrMuscl::advectToFaces(EBAMRFluxData& a_facePhi, const EBAMRCellData& a_phi, const Real a_extrapDt){
  CH_TIME("CdrMuscl::advectToFaces");
  if(m_verbosity > 5){
    pout() << m_name + "::advectToFaces" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;
  const int finest_level = m_amr->getFinestLevel();


  // Need to overwrite ghost cells, so I do that here.
  EBAMRCellData copy_state; 
  m_amr->allocate(copy_state, m_realm, m_phase, ncomp);
  DataOps::setValue(copy_state, 0.0);
  DataOps::incr(copy_state, a_phi, 1.0);

  m_amr->averageDown(copy_state,     m_realm, m_phase);
  m_amr->interpGhostPwl(copy_state, m_realm, m_phase);

  DataOps::setValue(a_facePhi, 0.0);

  for (int lvl = 0; lvl <= finest_level; lvl++){
    const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[lvl];
    const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[lvl];
    const ProblemDomain& domain  = m_amr->getDomains()[lvl];

    for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
      EBFluxFAB& face_state  = (*a_facePhi[lvl])[dit()];
      const Box box          = dbl.get(dit());
      const EBCellFAB& state = (*copy_state[lvl])[dit()];
      const EBFluxFAB& velo  = (*m_faceVelocity[lvl])[dit()];
      const EBISBox& ebisbox = state.getEBISBox();

      // Limit slopes and solve Riemann problem
      Box grown_box = box;
      grown_box.grow(1);
      EBCellFAB deltaC(ebisbox, grown_box, SpaceDim); // Cell-centered slopes
      deltaC.setVal(0.0);
      if(m_slopelim){
	this->computeSlopes(deltaC, state, box, domain, lvl, dit());
      }
      this->upwind(face_state, deltaC, state, velo, domain, box, lvl, dit());
    }

    this->computeBoundaryOutflow(*a_facePhi[lvl], lvl);
  }
}

void CdrMuscl::allocateInternals(){
  CdrTGA::allocateInternals();

  if(m_isDiffusive){
    this->setupDiffusionSolver();
  }
}

void CdrMuscl::computeSlopes(EBCellFAB&           a_deltaC,
			       const EBCellFAB&     a_phi,
			       const Box&           a_box,
			       const ProblemDomain& a_domain,
			       const int            a_level,
			       const DataIndex&     a_dit){
  CH_TIME("CdrMuscl::computeSlopes");
  if(m_verbosity > 99){
    pout() << m_name + "::computeSlopes" << endl;
  }

  const int comp         = 0;
  const int ncomp        = 1;

  // Regular cells
  BaseFab<Real>& regDeltaC      = a_deltaC.getSingleValuedFAB();
  const BaseFab<Real>& regState = a_phi.getSingleValuedFAB();


  for (int dir = 0; dir < SpaceDim; dir++){
    Box cellbox = a_box;
    cellbox.grow(1);
    FORT_MUSCL_SLOPES(CHF_FRA1(regDeltaC, dir),
		      CHF_CONST_FRA1(regState, comp),
		      CHF_CONST_INT(dir),
		      CHF_BOX(cellbox));
  }

  // Irregular cells
  const Box& domain_box      = a_domain.domainBox();
  const EBISBox& ebisbox     = a_phi.getEBISBox();
  const IntVectSet irreg_ivs = ebisbox.getIrregIVS(a_box);
  const EBGraph& ebgraph     = ebisbox.getEBGraph();
  VoFIterator vofit(irreg_ivs, ebgraph);
  
  for (int dir = 0; dir < SpaceDim; dir++){
    for (vofit.reset(); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();
      const IntVect iv    = vof.gridIndex();

      const bool on_lo_side = (iv[dir] == domain_box.smallEnd(dir));
      const bool on_hi_side = (iv[dir] == domain_box.bigEnd(dir));
      
      const bool has_faces_left = (ebisbox.numFaces(vof, dir, Side::Lo) == 1) && !on_lo_side;
      const bool has_faces_righ = (ebisbox.numFaces(vof, dir, Side::Hi) == 1) && !on_hi_side;

      VolIndex vof_left;
      VolIndex vof_righ;
      
      Real dwc   = 0.0;
      Real dwl   = 0.0;
      Real dwr   = 0.0;
      Real dwmin = 0.0;

      Real state_left = 0.0;
      Real state_righ = 0.0;
      Real state_cent = a_phi(vof, comp);

      // Compute left and right slope
      if(has_faces_left){
	Vector<FaceIndex> faces_left = ebisbox.getFaces(vof, dir, Side::Lo);
	vof_left   = faces_left[0].getVoF(Side::Lo);
	state_left = a_phi(vof_left, comp);
	dwl        = state_cent - state_left;
      }
      if(has_faces_righ){
	Vector<FaceIndex> faces_righ = ebisbox.getFaces(vof, dir, Side::Hi);
	vof_righ   = faces_righ[0].getVoF(Side::Hi);
	state_righ = a_phi(vof_righ, comp);
	dwr        = state_righ - state_cent;
      }

      if(has_faces_left && has_faces_righ){
	dwc = 0.5*(dwl + dwr);
      }
      else if(!has_faces_left && !has_faces_righ){
	dwl = 0.0;
	dwc = 0.0;
	dwr = 0.0;
      }
      else if(!has_faces_left && has_faces_righ){
	dwl = dwr;
	dwc = dwr;
      }
      else if(has_faces_left && !has_faces_righ){
	dwr = dwl;
	dwc = dwl;
      }
      else{
	MayDay::Abort("CdrMuscl::computeSlopes - missed a case");
      }


      // Limit slopes
      dwmin = 2.0*Min(Abs(dwl), Abs(dwr));
      if(dwl*dwr < 0.0){
	dwc = 0.0;
      }
      else{
	if(dwc < 0.0){
	  dwc = -Min(dwmin, Abs(dwc));
	}
	else if (dwc > 0.0){
	  dwc = Min(dwmin, Abs(dwc));
	}
	else{
	  dwc = 0.0;
	}
      }

      if(dwc != dwc){
	MayDay::Abort("CdrMuscl::computeSlopes - dwc != dwc.");
      }

      a_deltaC(vof, comp) = dwc;
    }
  }
}

void CdrMuscl::upwind(EBFluxFAB&           a_facePhi,
		       const EBCellFAB&     a_slopes,
		       const EBCellFAB&     a_phi,
		       const EBFluxFAB&     a_velo,
		       const ProblemDomain& a_domain,
		       const Box&           a_box,
		       const int&           a_level,
		       const DataIndex&     a_dit){
  CH_TIME("CdrMuscl::upwind");
  if(m_verbosity > 99){
    pout() << m_name + "::upwind" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  const EBISBox& ebisbox = a_phi.getEBISBox();
  const EBGraph& ebgraph = ebisbox.getEBGraph();
  const IntVectSet irreg = ebisbox.getIrregIVS(a_box);
  
  VoFIterator vofit(irreg, ebgraph);

  for (int dir = 0; dir < SpaceDim; dir++){
    BaseFab<Real>& regFaceStates   = a_facePhi[dir].getSingleValuedFAB();
    const BaseFab<Real>& regSlopes = a_slopes.getSingleValuedFAB();
    const BaseFab<Real>& regStates = a_phi.getSingleValuedFAB();
    const BaseFab<Real>& regVelo   = a_velo[dir].getSingleValuedFAB();
    const Box& face_box            = a_facePhi.getRegion();

    Box facebox = a_box;
    facebox.surroundingNodes();
    // Regular cells
    regFaceStates.setVal(0.0);
    FORT_MUSCL_UPWIND(CHF_FRA1(regFaceStates,   comp),
		      CHF_CONST_FRA1(regSlopes, dir),
		      CHF_CONST_FRA1(regStates, comp),
		      CHF_CONST_FRA1(regVelo,   comp),
		      CHF_CONST_INT(dir),
		      CHF_BOX(facebox));

    // Irregular cells
    EBFaceFAB& face_states     = a_facePhi[dir];
    const EBFaceFAB& face_velo = a_velo[dir];
    
    // The box that was sent in was cell-centered. We need to grow it by 1 in order to get 
    // all the cells around this box
    Box box = a_box;
    box.grow(1, dir);
    Vector<FaceIndex> multi_faces = ebgraph.getIrregFaces(box, dir);
    for (int i = 0; i < multi_faces.size(); i++){
      const FaceIndex& face = multi_faces[i];
      if(!face.isBoundary()){
	const VolIndex& vof_left = face.getVoF(Side::Lo);
	const VolIndex& vof_righ = face.getVoF(Side::Hi);

	const Real velo      = a_velo[dir](face, comp);
	const Real prim_left = a_phi(vof_left, comp) + 0.5*a_slopes(vof_left, dir);
	const Real prim_righ = a_phi(vof_righ, comp) - 0.5*a_slopes(vof_righ, dir);

	if(velo > 0.0){
	  face_states(face, comp) = prim_left;
	}
	else if (velo < 0.0){
	  face_states(face, comp) = prim_righ;
	}
	else{
	  face_states(face, comp) = 0.0;
	}
      }
    }
  }
}

void CdrMuscl::computeBoundaryOutflow(LevelData<EBFluxFAB>&       a_flux,
				      const int                   a_lvl){
  CH_TIME("CdrTGA::compute_sg_flux(lvl)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_sg_flux(lvl)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout& dbl = m_amr->getGrids(m_realm)[a_lvl];
  const ProblemDomain& domain  = m_amr->getDomains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->getEBISLayout(m_realm, m_phase)[a_lvl];
  const Real dx                = m_amr->getDx()[a_lvl];
  const Real zero              = 0.0;

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();

    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flux = a_flux[dit()][dir];

      // Get boundary domains
      Box lobox, hibox;
      int has_lo, has_hi;
      EBArith::loHi(lobox, has_lo, hibox, has_hi, domain, box, dir);
      const FaceStop::WhichFaces stopcrit = FaceStop::AllBoundaryOnly; // Only outflow

      // Lo faces
      if(has_lo == 1){
	const IntVectSet ivs(lobox);
	for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face = faceit();

	  const VolIndex vof_hi   = face.getVoF(Side::Hi);
	  const Vector<VolIndex> vof_hihi = ebisbox.getVoFs(vof_hi, dir, Side::Hi, 1);

	  Vector<FaceIndex> faces = ebisbox.getFaces(vof_hi, dir, Side::Hi);

	  Real f1 = 0.;
	  for (int i = 0; i < faces.size(); i++){
	    const FaceIndex& iface = faces[i];
	    f1 += flux(iface, comp)*ebisbox.areaFrac(iface);
	  }

	  Real f2 = 0.;
	  for (int ivof = 0; ivof < vof_hihi.size(); ivof++){
	    Vector<FaceIndex> faces = ebisbox.getFaces(vof_hihi[ivof], dir, Side::Hi);
	    for (int i = 0; i < faces.size(); i++){
	      const FaceIndex& iface = faces[i];
	      f2 += flux(iface, comp)*ebisbox.areaFrac(iface);
	    }
	  }

	  // Extrapolate flux
	  const Real extrap = 2.0*f1 - f2;

	  flux(face, comp) = Min(extrap, zero); // No inflow
	  flux(face, comp) = 0.0;
	}
      }

      // Hi faces
      if(has_hi == 1){
	const IntVectSet ivs(hibox);
	for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	  const FaceIndex& face   = faceit();
	  const VolIndex vof_lo   = face.getVoF(Side::Lo);
	  const Vector<VolIndex> vof_lolo = ebisbox.getVoFs(vof_lo, dir, Side::Lo, 1);

	  Vector<FaceIndex> faces = ebisbox.getFaces(vof_lo, dir, Side::Lo);

	  Real f1 = 0.;
	  for (int i = 0; i < faces.size(); i++){
	    const FaceIndex& iface = faces[i];
	    f1 += flux(iface, comp)*ebisbox.areaFrac(iface);
	  }

	  Real f2 = 0.;
	  for (int ivof = 0; ivof < vof_lolo.size(); ivof++){
	    Vector<FaceIndex> faces = ebisbox.getFaces(vof_lolo[ivof], dir, Side::Lo);
	    for (int i = 0; i < faces.size(); i++){
	      const FaceIndex& iface = faces[i];
	      f2 += flux(iface, comp)*ebisbox.areaFrac(iface);
	    }
	  }

	  // Extrapolate flux
	  const Real extrap = 2.0*f1 - f2;

	  flux(face, comp) = Max(extrap, zero); // No inflow
	  flux(face, comp) = 0.0;
	}
      }
    }
  }
}

#include <CD_NamespaceFooter.H>
