/*!
  @file   ConductivityDomainBcWrapper.cpp
  @brief  Implementation of ConductivityDomainBcWrapper.H
  @author Robert Marskar
  @date   June. 2018
*/

#include <CD_ConductivityDomainBcWrapper.H>
#include <CD_RobinConductivityDomainBc.H>

#include "CD_NamespaceHeader.H"

ConductivityDomainBcWrapper::ConductivityDomainBcWrapper(){
  m_bc.resize(2*SpaceDim);

  m_defined = false;
}

ConductivityDomainBcWrapper::~ConductivityDomainBcWrapper(){

}

void ConductivityDomainBcWrapper::setCoefficients(){
  for (int i = 0; i < m_bc.size(); i++){
    m_bc[i]->setCoef(this->m_eblg, this->m_beta, this->m_bcoef);
  }
  m_defined = true;
}

void ConductivityDomainBcWrapper::setPotentials(const Vector<RefCountedPtr<BaseBCFuncEval> >& a_potentials){
  m_potentials = a_potentials;
}

void ConductivityDomainBcWrapper::setRobinCoefficients(const Vector<RefCountedPtr<RobinCoefficients> >& a_robinco){
  m_robinco = a_robinco;
}

void ConductivityDomainBcWrapper::setWallBc(const Vector<RefCountedPtr<wall_bc> >& a_wallbc){
  for (int i = 0; i < a_wallbc.size(); i++){
    const int dir             = a_wallbc[i]->direction();
    const Side::LoHiSide side = a_wallbc[i]->side();
    const int idx             = wall_bc::map_bc(dir, side);
      
    if(a_wallbc[i]->which_bc() == wallbc::dirichlet){
      m_bc[idx] = RefCountedPtr<DirichletConductivityDomainBC> (new DirichletConductivityDomainBC());
      if(a_wallbc[i]->is_live()){
	m_bc[idx]->setFunction(m_potentials[i]);
      }
      else{
	m_bc[idx]->setValue(0.0);
      }
    }
    else if(a_wallbc[i]->which_bc() == wallbc::neumann){
      m_bc[idx] = RefCountedPtr<NeumannConductivityDomainBC> (new NeumannConductivityDomainBC());
      m_bc[idx]->setValue(a_wallbc[i]->get_value());
    }
    else if(a_wallbc[i]->which_bc() == wallbc::robin){
      // For the Poisson equation, the appropriate coefficients are (1,-1,0.0) for the robin BC class (see robincondu*.H).
      // Otherwise, use externally supplied values
      RobinConductivityDomainBc* robinbc = new RobinConductivityDomainBc();
      if(m_robinco[i] == NULL){
	robinbc->setCoefficients(1.0, -1.0, 0.0);
      }
      else{
	robinbc->setCoefficients(m_robinco[i]);
      }
	    
      m_bc[idx] = RefCountedPtr<RobinConductivityDomainBc> (robinbc);
      m_bc[idx]->setValue(0.0);
    }
    else{
      MayDay::Abort("ConductivityDomainBcWrapper::set_bc(wall) - unsupported bc type requested");
    }

  }
}

void ConductivityDomainBcWrapper::getFaceFlux(BaseFab<Real>&        a_faceFlux,
					       const BaseFab<Real>&  a_phi,
					       const RealVect&       a_probLo,
					       const RealVect&       a_dx,
					       const int&            a_idir,
					       const Side::LoHiSide& a_side,
					       const DataIndex&      a_dit,
					       const Real&           a_time,
					       const bool&           a_useHomogeneous){
  if(!m_defined){
    this->setCoefficients(); // I really hate that I have to do this, but there's no entry point in EbHelmholtzOp (yet) that
  }                   // allows me to do this at constructor level
  const int idx = wall_bc::map_bc(a_idir, a_side);
  m_bc[idx]->getFaceFlux(a_faceFlux,
			 a_phi,
			 a_probLo,
			 a_dx,
			 a_idir,
			 a_side,
			 a_dit,
			 a_time,
			 a_useHomogeneous);
}

void ConductivityDomainBcWrapper::getFaceFlux(Real&                 a_faceFlux,
					       const VolIndex&       a_vof, 
					       const int&            a_comp, 
					       const EBCellFAB&      a_phi, 
					       const RealVect&       a_probLo, 
					       const RealVect&       a_dx, 
					       const int&            a_idir, 
					       const Side::LoHiSide& a_side, 
					       const DataIndex&      a_dit, 
					       const Real&           a_time, 
					       const bool&           a_useHomogeneous){
  if(!m_defined){
    this->setCoefficients(); // I really hate that I have to do this, but there's no entry point in EbHelmholtzOp (yet) that
  }                   // allows me to do this at constructor level
  const int idx = wall_bc::map_bc(a_idir, a_side);
  m_bc[idx]->getFaceFlux(a_faceFlux,
			 a_vof,
			 a_comp,
			 a_phi,
			 a_probLo,
			 a_dx,
			 a_idir,
			 a_side,
			 a_dit,
			 a_time,
			 a_useHomogeneous);
}

void ConductivityDomainBcWrapper::getFaceGradPhi(Real&                 a_faceFlux,
						  const FaceIndex&      a_face,
						  const int&            a_comp,
						  const EBCellFAB&      a_phi,
						  const RealVect&       a_probLo,
						  const RealVect&       a_dx,
						  const int&            a_idir,
						  const Side::LoHiSide& a_side,
						  const DataIndex&      a_dit,
						  const Real&           a_time,
						  const bool&           a_useAreaFrac,
						  const RealVect&       a_centroid,
						  const bool&           a_useHomogeneous){
  if(!m_defined){
    this->setCoefficients(); // I really hate that I have to do this, but there's no entry point in EbHelmholtzOp (yet) that
  }                   // allows me to do this at constructor level
  const int idx = wall_bc::map_bc(a_idir, a_side);
  m_bc[idx]->getFaceGradPhi(a_faceFlux,
			    a_face,
			    a_comp,
			    a_phi,
			    a_probLo,
			    a_dx,
			    a_idir,
			    a_side,
			    a_dit,
			    a_time,
			    a_useAreaFrac,
			    a_centroid,
			    a_useHomogeneous);
}

void ConductivityDomainBcWrapper::fillPhiGhost(FArrayBox&     a_phi,
						const Box&     a_valid,
						const Box&     a_domain,
						Real           a_dx,
						bool           a_homogeneous){
  Box grownBox = a_valid;
  grownBox.grow(1);

  MayDay::Abort("ConductivityDomainBcWrapper::fillPhiGhost -- how did this get called..?");

  for (int idir=0; idir<CH_SPACEDIM; ++idir)
    {
      for(SideIterator sit; sit.ok(); ++sit)
	{
	  Box choppedBox = grownBox;
	  choppedBox.grow(idir,-1);
	  Box toRegion = adjCellBox(choppedBox, idir, sit(), 1);

	  if(!a_domain.contains(toRegion))
	    {
	      for (BoxIterator bit(toRegion); bit.ok(); ++bit)
		{
		  const IntVect& iv = bit();
		  //fake vof just to get the location
		  VolIndex vof(iv, 0);
		  RealVect loc = EBArith::getVoFLocation(vof, a_dx, RealVect::Zero);
		  int isign = sign(sit());
		  IntVect ivneigh = iv - isign*BASISV(idir);
		  Real value = bcvaluefunc(loc, idir, sit());
		  if(a_homogeneous) value = 0;
		  a_phi(iv, 0) = a_phi(ivneigh, 0)  + a_dx*value;
		}
	    }
	} 
    }//end loop over directions
}
#include "CD_NamespaceFooter.H"
