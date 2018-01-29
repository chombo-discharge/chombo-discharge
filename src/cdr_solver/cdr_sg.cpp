/*!
  @file   cdr_sg.cpp
  @brief  Implementation of cdr_SG.H
  @author Robert Marskar
  @date   Jan. 2018
*/

#include "cdr_sg.H"
#include "data_ops.H"

#include <EBArith.H>

cdr_sg::cdr_sg() : cdr_solver() {

  m_name = "cdr_sg";

  this->set_mass_redist(false);
}

cdr_sg::~cdr_sg(){
  
}

int cdr_sg::query_ghost() const {
  CH_TIME("cdr_sg::query_ghost");
  if(m_verbosity > 5){
    pout() << m_name + "::query_ghost" << endl;
  }
  return 3;
}

void cdr_sg::compute_divJ(EBAMRCellData& a_divJ, const EBAMRCellData& a_state, const Real a_extrap_dt){
  CH_TIME("cdr_sg::compute_divJ");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_divJ" << endl;
  }

  const int comp       = 0;
  const int ncomp      = 1;
  const int redist_rad = m_amr->get_redist_rad();

  EBAMRFluxData flux;
  EBAMRIVData   div_nc;
  EBAMRIVData   mass_diff;
  EBAMRCellData weights;

  m_amr->allocate(flux,      m_phase, ncomp);
  m_amr->allocate(div_nc,    m_phase, ncomp);
  m_amr->allocate(mass_diff, m_phase, ncomp);
  m_amr->allocate(weights,   m_phase, ncomp, 2*redist_rad);

  data_ops::set_value(a_divJ,    0.0); 
  data_ops::set_value(flux,      0.0);
  data_ops::set_value(div_nc,    0.0);
  data_ops::set_value(mass_diff, 0.0);
  data_ops::set_value(weights,   0.0);

  this->average_velo_to_faces(m_velo_face, m_velo_cell);       // Average cell-centered velocities to faces
  this->compute_sg_flux(flux, a_state, m_velo_face, m_diffco); // Compute Scharfetter-Gummel flux
  this->conservative_divergence(a_divJ, flux);                 // Compute conservative divergence
  this->nonconservative_divergence(div_nc, a_divJ, flux);      // Compute non-conservative divergence. Last argument is a dummy.
  this->hybrid_divergence(a_divJ, mass_diff, div_nc);          // Make divJ = hybrid divergence. Compute mass diff.
  this->increment_flux_register(flux);                         // Increment flux registers
  this->increment_redist(mass_diff);                           // Increment redistribution objects

  // Mass weights
  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    data_ops::incr(*weights[lvl], *a_state[lvl], 1.0);
  }

  const bool ebcf = m_amr->get_ebcf();
  if(ebcf){ 
    this->coarse_fine_increment(mass_diff);     // Increment the coarse-fine redistribution objects
    this->increment_redist_flux();              // Increment flux registers with the redistribution stuff
    this->coarse_fine_redistribution(a_divJ);   // Redistribute
  }
  else{
    this->hyperbolic_redistribution(a_divJ, mass_diff, weights);  // Redistribute mass into hybrid divergence
    this->reflux(a_divJ);                                         // Reflux at coarse-fine interfaces
  }

  for (int lvl = 0; lvl <= m_amr->get_finest_level(); lvl++){
    a_divJ[lvl]->exchange();
  }
}

void cdr_sg::compute_sg_flux(EBAMRFluxData&       a_flux,
			     const EBAMRCellData& a_state,
			     const EBAMRFluxData& a_velo,
			     const EBAMRFluxData& a_diffco){
  CH_TIME("cdr_tga::compute_sg_flux");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_sg_flux" << endl;
  }

  const int finest_level = m_amr->get_finest_level();
  for (int lvl = 0; lvl <= finest_level; lvl++){
    this->compute_sg_flux(*a_flux[lvl], *a_state[lvl], *a_velo[lvl], *a_diffco[lvl], lvl);
    this->compute_bndry_outflow(*a_flux[lvl], lvl);
    a_flux[lvl]->exchange();
  }

  m_amr->average_down(a_flux, m_phase);
}

void cdr_sg::compute_sg_flux(LevelData<EBFluxFAB>&       a_flux,
			     const LevelData<EBCellFAB>& a_state,
			     const LevelData<EBFluxFAB>& a_velo,
			     const LevelData<EBFluxFAB>& a_diffco,
			     const int                   a_lvl){
  CH_TIME("cdr_tga::compute_sg_flux(lvl)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_sg_flux(lvl)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];
  const Real dx                = m_amr->get_dx()[a_lvl];

  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box box          = dbl.get(dit());
    const EBCellFAB& state = a_state[dit()];
    const EBISBox& ebisbox = ebisl[dit()];
    const EBGraph& ebgraph = ebisbox.getEBGraph();
    const IntVectSet ivs(box);

    for (int dir = 0; dir < SpaceDim; dir++){
      EBFaceFAB& flux         = a_flux[dit()][dir];
      const EBFaceFAB& velo   = a_velo[dit()][dir];
      const EBFaceFAB& diffco = a_diffco[dit()][dir];

      const FaceStop::WhichFaces stopcrit = FaceStop::SurroundingNoBoundary; // Boundary fluxes from elsewhere. 
      for (FaceIterator faceit(ivs, ebgraph, dir, stopcrit); faceit.ok(); ++faceit){
	const FaceIndex& face = faceit();
	const VolIndex vof_lo = face.getVoF(Side::Lo);
	const VolIndex vof_hi = face.getVoF(Side::Hi);
	  
	flux(face, comp) = 0.;

	const Real n_lo     = state(vof_lo, comp);
	const Real n_hi     = state(vof_hi, comp);
	const Real v        = velo(face,   comp);

	if(m_diffusive){
	  const Real d        = diffco(face, comp);
	  const Real exp_plus = exp( v*dx/d);
	  const Real exp_minu = exp(-v*dx/d);

	  flux(face, comp) = v*(n_lo*exp_plus - n_hi*exp_minu)/(exp_plus - exp_minu);
	}
	else{
	  if(v > 0.){
	    flux(face, comp) = n_lo*v;
	  }
	  else if(v < 0.){
	    flux(face, comp) = n_hi*v;
	  }
	}
      }
    }
  }
}

void cdr_sg::compute_bndry_outflow(LevelData<EBFluxFAB>&       a_flux,
				   const int                   a_lvl){
  CH_TIME("cdr_tga::compute_sg_flux(lvl)");
  if(m_verbosity > 5){
    pout() << m_name + "::compute_sg_flux(lvl)" << endl;
  }

  const int comp  = 0;
  const int ncomp = 1;

  const DisjointBoxLayout& dbl = m_amr->get_grids()[a_lvl];
  const ProblemDomain& domain  = m_amr->get_domains()[a_lvl];
  const EBISLayout& ebisl      = m_amr->get_ebisl(m_phase)[a_lvl];
  const Real dx                = m_amr->get_dx()[a_lvl];

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
      const FaceStop::WhichFaces stopcrit = FaceStop::AllBoundaryOnly; // Boundary fluxes from elsewhere.

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

	  flux(face, comp) = Min(extrap, 0.); // No inflow
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

	  flux(face, comp) = Max(extrap, 0.); // No inflow
	}
      }
    }
  }
}
