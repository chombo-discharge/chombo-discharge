/*!
  @file mfdirichletconductivityebbc.cpp
  @brief Implementation of mfdiricheltconductivityebbc.H
  @author Robert Marskar
  @date Dec. 2017
*/

#include "mfdirichletconductivityebbc.H"

bool mfdirichletconductivityebbc::s_quadrant_based = true;
int  mfdirichletconductivityebbc::s_lsq_radius     = 2;

mfdirichletconductivityebbc::mfdirichletconductivityebbc(const ProblemDomain& a_domain,
							 const EBISLayout&    a_ebisl,
							 const RealVect&      a_dx,
							 const IntVect*       a_phig,
							 const IntVect*       a_rhsg) : DirichletConductivityEBBC(a_domain,
														  a_ebisl,
														  a_dx,
														  a_phig,
														  a_rhsg){

  m_domain = a_domain;
  m_ebisl  = a_ebisl;
  m_dx     = a_dx;

  this->setOrder(1);
}

mfdirichletconductivityebbc::~mfdirichletconductivityebbc(){

}

#if 1
LayoutData<BaseIVFAB<VoFStencil> >* mfdirichletconductivityebbc::getFluxStencil(int ivar){
  return &m_irreg_stencils;
}
#endif

void mfdirichletconductivityebbc::setOrder(int a_order){
  m_order = a_order;
}

void mfdirichletconductivityebbc::define_ivs(const MFLevelGrid& a_mflg){
    
  m_ivs.define(a_mflg.get_grids());

  const DisjointBoxLayout& dbl = a_mflg.get_grids();

  int num = 0;
  m_ivs.define(dbl);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    m_ivs[dit()] = a_mflg.interface_region(dbl.get(dit()), dit());
    num += m_ivs[dit()].numPts();
  }
 
  m_definedmf = true;
}

void mfdirichletconductivityebbc::define(const LayoutData<IntVectSet>& a_cfivs, const Real& a_factor){

  if(!m_definedmf){
    MayDay::Error("mfdirichletconductivityebbc::define - must set multifluid cells first!");
  }
  else if(!m_coefSet){
    MayDay::Error("mfdirichletconductivityebbc::define - must call setCoef BEFORE calling define");
  }
  
#if 0 // Fall back to dirichletconductivityebbc
  DirichletConductivityEBBC::define(a_cfivs, a_factor);
#else // This is where the new code goes.
  const int comp      = 0;
  const int num_comps = 1;

  const DisjointBoxLayout& dbl = m_ebisl.getDisjointLayout();

  m_irreg_weights.define(dbl);
  m_matching_weights.define(dbl);
  m_irreg_stencils.define(dbl);
  m_matching_stencils.define(dbl);

  // Create stencils on irregular cells that are not part of m_ivs
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    const Box& box          = dbl[dit()];
    const EBISBox& ebisbox  = m_ebisl[dit()];
    const EBGraph& ebgraph  = ebisbox.getEBGraph();
    const IntVectSet& cfivs = a_cfivs[dit()];

    IntVectSet irreg_ivs;  // All irregular cells
    IntVectSet diri_ivs;   // Pure dirichlet cells (i.e. non-matched irregular cells)
    IntVectSet match_ivs;  // Quasi-dirichlet cells (i.e. matched irregular cells)
    
    irreg_ivs  = ebisbox.getIrregIVS(box);
    irreg_ivs |= ebisbox.getMultiCells(box);

    diri_ivs  = irreg_ivs - m_ivs[dit()];
    match_ivs = m_ivs[dit()];

    m_irreg_weights[dit()].define(irreg_ivs,  ebgraph, num_comps);
    m_irreg_stencils[dit()].define(irreg_ivs, ebgraph, num_comps);

    m_matching_weights[dit()].define(match_ivs,  ebgraph, num_comps);
    m_matching_stencils[dit()].define(match_ivs, ebgraph, num_comps);

    // Build stencils for pure Dirichlet type irreg cells
#if 1 // This is what it is supposed to look like
    for (VoFIterator vofit(diri_ivs, ebgraph); vofit.ok(); ++vofit){
#else // Define for all irreg, this should now emulate DirichletConductivityEBBC::define
      for (VoFIterator vofit(irreg_ivs, ebgraph); vofit.ok(); ++vofit){
#endif
      const VolIndex& vof = vofit();

      Real& cur_weight        = m_irreg_weights[dit()](vof, comp);
      VoFStencil& cur_stencil = m_irreg_stencils[dit()](vof, comp);

      bool drop_order = false;

      //m_order = 2;
      if(m_order == 2){
	drop_order = this->get_second_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
      }
      else if(m_order == 1){
	this->get_first_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
      }

      if(m_order == 2 && drop_order){
	this->get_first_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
      }

      // Stencil should be scaled by bco*beta. We can do this with the weight as well, as long
      // as we remember to skip this during applyEBFlux. 
      const Real factor = m_beta*(*m_bcoe)[dit()](vof, comp);
      cur_stencil *= factor;
      cur_weight  *= factor;
    }

    // Stencils for matching cells
    for (VoFIterator vofit(match_ivs, ebgraph); vofit.ok(); ++vofit){
      const VolIndex& vof = vofit();

      Real& cur_weight        = m_matching_weights[dit()](vof, comp);
      VoFStencil& cur_stencil = m_matching_stencils[dit()](vof, comp);

      bool drop_order = false;
	
      if(m_order == 2){
	drop_order = this->get_second_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
      }
      else if(m_order == 1){
	this->get_first_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
      }

      if(m_order == 2 && drop_order){
	this->get_first_order_sten(cur_weight, cur_stencil, vof, ebisbox, cfivs);
      }



      // Should we also scale this stencil and weights? Not so sure....
      const Real factor = m_beta*(*m_bcoe)[dit()](vof, comp);
      cur_stencil *= factor;
      cur_weight  *= factor;
    }
  }
#endif
}

bool mfdirichletconductivityebbc::get_second_order_sten(Real&             a_weight,
							VoFStencil&       a_stencil,
							const VolIndex&   a_vof,
							const EBISBox&    a_ebisbox,
							const IntVectSet& a_cfivs){
  CH_TIME("mfdirichletconductivityebbc::get_second_order_sten");

  a_stencil.clear();
  bool drop_order = false;

  Vector<VoFStencil> point_stencils;
  Vector<Real> distance_along_lines;
  
  EBArith::johanStencil(drop_order, point_stencils, distance_along_lines, a_vof, a_ebisbox, m_dx, a_cfivs);
  if(drop_order){
    return true;
  }

  // If we got this far we have a stencil
  CH_assert(distance_along_lines.size() >= 2);
  CH_assert(point_stencils.size() >= 2);

  const Real& x1   = distance_along_lines[0];
  const Real& x2   = distance_along_lines[1];
  const Real denom = x2*x2*x1 - x1*x1*x2;

  VoFStencil& phi1Sten = point_stencils[0];
  VoFStencil& phi2Sten = point_stencils[1];
  
  phi1Sten *= -x2*x2/denom;
  phi2Sten *=  x1*x1/denom;

  a_weight   = -x1*x1/denom + x2*x2/denom;
  a_stencil +=  phi1Sten;
  a_stencil +=  phi2Sten;
}

void mfdirichletconductivityebbc::get_first_order_sten(Real&             a_weight,
						       VoFStencil&       a_stencil,
						       const VolIndex&   a_vof,
						       const EBISBox&    a_ebisbox,
						       const IntVectSet& a_cfivs){
  CH_TIME("mfdirichletconductivityebbc::get_first_order_sten");

  const RealVect& normal   = a_ebisbox.normal(a_vof);
  const RealVect& centroid = a_ebisbox.bndryCentroid(a_vof);

  if(s_quadrant_based){
    EBArith::getLeastSquaresGradSten(a_stencil, a_weight, a_vof, a_ebisbox, m_dx, m_domain, 0);
  }
  else{
    EBArith::getLeastSquaresGradStenAllVoFsRad(a_stencil,
					       a_weight,
					       normal,
					       centroid,
					       a_vof,
					       a_ebisbox,
					       m_dx*RealVect::Unit,
					       m_domain,
					       0,
					       s_lsq_radius);
  }
}

#if 0 // THIS IS OLD CODE
void mfdirichletconductivityebbc::applyEBFlux(EBCellFAB&                    a_lphi,
					      const EBCellFAB&              a_phi,
					      VoFIterator&                  a_vofit,
					      const LayoutData<IntVectSet>& a_cfivs,
					      const DataIndex&              a_dit,
					      const RealVect&               a_probLo,
					      const RealVect&               a_dx,
					      const Real&                   a_factor,
					      const bool&                   a_useHomogeneous,
					      const Real&                   a_time){

  const BaseIVFAB<Real>& poissWeight = (m_bc.getFluxWeight())[a_dit];

  Real value = 0.0;
  const EBISBox&   ebisBox = a_phi.getEBISBox();

#if 0
  pout() << "mfidirichletconductivityebbc::applyEBFlux - homogeneous = " << a_useHomogeneous << endl;
#endif

  // Homogeneous BCs are different for matching. We must use the supplied BC value through mfdirichletconductivityebbc
  if(a_useHomogeneous){ 
    for (a_vofit.reset(); a_vofit.ok(); ++a_vofit){
      const VolIndex& vof = a_vofit();

      if(m_dataBased){
	if(m_ivs[a_dit].contains(vof.gridIndex())){
	  value = (*m_data)[a_dit](vof, 0);  // Homogeneous version for variable-potential cells
	}
	else{
	  value = 0.; // Homogeneous version for fixed-potential cells
	}

	//	value = 0.;
#if 0
	if(Abs(value) > 1.E-4){
	  pout() << "value = " << value << endl;
	  MayDay::Abort("mfdirichletconductivityebbc::applyebflux - error");
	}
#endif
	const Real poissWeightPt = poissWeight(vof, 0);
	const Real& areaFrac     = ebisBox.bndryArea(vof);
	const Real& bcoef        = (*m_bcoe)[a_dit](vof,0);
	const Real flux          = poissWeightPt*value*areaFrac;
	const Real compFactor    = a_factor*bcoef*m_beta;
	a_lphi(vof,0)           += flux * compFactor;
      }
      else{
	MayDay::Abort("error");
      }
    }
  }
  else { // Non-homogeneous boundary conditions. We can use DirichletConductivityEBBC as long as the BCs are updated correctly. 
    DirichletConductivityEBBC::applyEBFlux(a_lphi,
					   a_phi,
					   a_vofit,
					   a_cfivs,
					   a_dit,
					   a_probLo,
					   a_dx,
					   a_factor,
					   a_useHomogeneous,
					   a_time);
  }
}
#else
void mfdirichletconductivityebbc::applyEBFlux(EBCellFAB&                    a_lphi,
					      const EBCellFAB&              a_phi,
					      VoFIterator&                  a_vofit,
					      const LayoutData<IntVectSet>& a_cfivs,
					      const DataIndex&              a_dit,
					      const RealVect&               a_probLo,
					      const RealVect&               a_dx,
					      const Real&                   a_factor,
					      const bool&                   a_useHomogeneous,
					      const Real&                   a_time){
  CH_assert(m_dataBased);
  
  const int comp = 0;
  
  if(a_useHomogeneous){
    // Do nothing for now. This changes for matched cells. 
  }
  else{

    const EBISBox& ebisbox = a_phi.getEBISBox();

    for (a_vofit.reset(); a_vofit.ok(); ++a_vofit){
      const VolIndex& vof = a_vofit();

      const Real& value     = (*m_data)[a_dit](vof, comp);
      const Real& weight    = m_irreg_weights[a_dit](vof, comp); // Already contains m_beta*m_bcoe
      const Real& area_frac = ebisbox.bndryArea(vof);
      const Real flux       = weight*value*area_frac;

      a_lphi(vof, comp) += flux*a_factor;
    }
  }
}
#endif
