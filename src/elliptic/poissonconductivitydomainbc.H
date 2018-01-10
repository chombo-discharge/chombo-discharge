/*!
  @file poissonconductivitydomainbc.H
  @brief Declaration of a class that wraps boundary Dirichlet and Neumann boundary conditions on individual domain edges
  @author Robert Marskar
  @date Jan. 2018
  @todo Segregate implementation file
*/

#ifndef _POISSONCONDUCTIVITYDOMAINBC_
#define _POISSONCONDUCTIVITYDOMAINBC_

#include <BaseDomainBC.H>
#include <DirichletConductivityDomainBC.H>
#include <NeumannConductivityDomainBC.H>

/*!
  @brief Wrapper class for imposing Dirichlet and Neumann on individual domain edges (faces in 3D)
  @details This class does not (currently) do Robin boundary conditions, although it can probably be extended to take care of 
  those cases as well. 
*/
class poissonconductivitydomainbc : public ConductivityBaseDomainBC {
public:
  /*!
    @brief Constructor
  */
  poissonconductivitydomainbc(){

    m_diri = RefCountedPtr<DirichletConductivityDomainBC> (new DirichletConductivityDomainBC());
    m_neum = RefCountedPtr<NeumannConductivityDomainBC> (new NeumannConductivityDomainBC());

    this->set_value(0.0);

    m_defined = false;
  }

  /*!
    @brief Destructor
  */
  ~poissonconductivitydomainbc(){}

  /*!
    @brief Dirichlet and Neumann need to know about bcoef
  */
  void set_coef(){
    m_diri->setCoef(this->m_eblg, this->m_beta, this->m_bcoef);
    m_neum->setCoef(this->m_eblg, this->m_beta, this->m_bcoef);
    
    m_defined = true;
  }

  /*!
    @brief Set function
  */
  virtual void set_function(RefCountedPtr<BaseBCFuncEval> a_func){
    m_diri->setFunction(a_func);
    m_neum->setFunction(a_func);
  }

  /*!
    @brief Set value
  */
  virtual void set_value(Real a_value){
    m_diri->setValue(a_value);
    m_neum->setValue(a_value);
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
    m_diri->getFaceFlux(a_faceFlux,
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
    m_diri->getFaceFlux(a_faceFlux,
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
    m_diri->getFaceGradPhi(a_faceFlux,
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
    @brief Dirichlet boundary condition
  */
  RefCountedPtr<DirichletConductivityDomainBC> m_diri;

  /*!
    @brief Neumann boundary condition
  */
  RefCountedPtr<NeumannConductivityDomainBC> m_neum;

  /*!
    @brief Is defined or not
  */
  bool m_defined;

};
#endif
