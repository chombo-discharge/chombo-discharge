/*!
  @file conductivitydomainbc_wrapper.H
  @brief Declaration of a class that wraps boundary Dirichlet and Neumann boundary conditions on individual domain edges
  @author Robert Marskar
  @date Jan. 2018
  @todo Segregate implementation file
*/

#ifndef _CONDUCTIVITYDOMAINBC_WRAPPER_
#define _CONDUCTIVITYDOMAINBC_WRAPPER_

#include "wall_bc.H"

#include <BaseDomainBC.H>
#include <DirichletConductivityDomainBC.H>
#include <NeumannConductivityDomainBC.H>

/*!
  @brief Wrapper class for imposing Dirichlet and Neumann on individual domain edges (faces in 3D)
  @details This class does not (currently) do Robin boundary conditions, although it can probably be extended to take care of 
  those cases as well. 
*/
class conductivitydomainbc_wrapper : public ConductivityBaseDomainBC {
public:
  /*!
    @brief Constructor
  */
  conductivitydomainbc_wrapper(){

    m_bc.resize(2*SpaceDim);

    m_defined = false;
  }

  /*!
    @brief Destructor
  */
  virtual ~conductivitydomainbc_wrapper(){}

  /*!
    @brief Dirichlet and Neumann need to know about bcoef
  */
  virtual void set_coef(){
    for (int i = 0; i < m_bc.size(); i++){
      m_bc[i]->setCoef(this->m_eblg, this->m_beta, this->m_bcoef);
    }
    m_defined = true;
  }

  virtual void set_potential(const RefCountedPtr<BaseBCFuncEval>& a_potential){
    m_potential = a_potential;
  }

  /*!
    @brief Set bc types. This should come through the factory. 
  */
  virtual void set_wallbc(const Vector<RefCountedPtr<wall_bc> >& a_wallbc){
    for (int i = 0; i < a_wallbc.size(); i++){
      const int dir             = a_wallbc[i]->direction();
      const Side::LoHiSide side = a_wallbc[i]->side();
      const int idx             = wall_bc::map_bc(dir, side);
      
      if(a_wallbc[i]->which_bc() == wallbc::dirichlet){
	m_bc[idx] = RefCountedPtr<DirichletConductivityDomainBC> (new DirichletConductivityDomainBC());

	if(a_wallbc[i]->is_live()){
	  m_bc[idx]->setFunction(m_potential);
	}
	else{
	  m_bc[idx]->setValue(0.0);
	}
      }
      else if(a_wallbc[i]->which_bc() == wallbc::neumann){
	m_bc[idx] = RefCountedPtr<NeumannConductivityDomainBC> (new NeumannConductivityDomainBC());
	m_bc[idx]->setValue(a_wallbc[i]->get_value());
      }
      else{
	MayDay::Abort("conductivitydomainbc_wrapper::set_bc(wall) - unsupported bc type requested");
      }

    }
  }

  /*! 
    @brief Get face flux
  */
  virtual void getFaceFlux(BaseFab<Real>&        a_faceFlux,
                           const BaseFab<Real>&  a_phi,
                           const RealVect&       a_probLo,
                           const RealVect&       a_dx,
                           const int&            a_idir,
                           const Side::LoHiSide& a_side,
                           const DataIndex&      a_dit,
                           const Real&           a_time,
                           const bool&           a_useHomogeneous){
    if(!m_defined){
      this->set_coef(); // I really hate that I have to do this, but there's no entry point in ebconductivityop (yet) that
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

  /*! 
    @brief Get face flux
  */
  virtual void getFaceFlux(Real&                 a_faceFlux,
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
      this->set_coef(); // I really hate that I have to do this, but there's no entry point in ebconductivityop (yet) that
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


  /*!
    @brief Get gradient
  */
  virtual void getFaceGradPhi(Real&                 a_faceFlux,
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
      this->set_coef(); // I really hate that I have to do this, but there's no entry point in ebconductivityop (yet) that
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

protected:

  /*!
    @brief Individual stuff
  */
  Vector<RefCountedPtr<ConductivityBaseDomainBC> > m_bc;

  /*!
    @brief Potential
  */
  RefCountedPtr<BaseBCFuncEval> m_potential;

  /*!
    @brief Is defined or not
  */
  bool m_defined;

};
#endif
