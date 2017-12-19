/*!
  @file mfdirichletconductivityebbc.cpp
  @brief Implementation of mfdiricheltconductivityebbc.H
  @author Robert Marskar
  @date Dec. 2017
*/

#include "mfdirichletconductivityebbc.H"

mfdirichletconductivityebbc::~mfdirichletconductivityebbc(){

}

void mfdirichletconductivityebbc::define_ivs(const MFLevelGrid& a_mflg){
    
  m_ivs.define(a_mflg.get_grids());

  const DisjointBoxLayout& dbl = a_mflg.get_grids();

  int num = 0;
  m_ivs.define(dbl);
  for (DataIterator dit = dbl.dataIterator(); dit.ok(); ++dit){
    m_ivs[dit()] = a_mflg.interface_region(dbl.get(dit()), dit());
  }
}

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
  
  // Homogeneous boundary conditions
  if(!a_useHomogeneous){
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
  else{ // Homogeneous BCs are different for matching. We must use the supplied BC value through jump_bc with the homogeneous
        // version of the matching criteria; i.e. b1*dnPhi + b2*dnPhi = 0. We assume that m_data has already been filled
        // with these values. 

  }
}
