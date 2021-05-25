/*!
  @file   ConductivityDomainBcWrapper.H
  @brief  Declaration of a class that wraps boundary Dirichlet and Neumann boundary conditions on individual domain edges
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _CONDUCTIVITYDOMAINBC_WRAPPER_
#define _CONDUCTIVITYDOMAINBC_WRAPPER_

#include <CD_WallBc.H>
#include <CD_RobinCoefficients.H>

#include <BaseDomainBC.H>
#include <DirichletConductivityDomainBC.H>
#include <NeumannConductivityDomainBC.H>

#include "CD_NamespaceHeader.H"

/*!
  @brief Wrapper class for imposing Dirichlet and Neumann on individual domain edges (faces in 3D)
  @details This class does not (currently) do Robin boundary conditions, although it can probably be extended to take care of 
  those cases as well. 
*/
class ConductivityDomainBcWrapper : public ConductivityBaseDomainBC {
public:
  /*!
    @brief Constructor
  */
  ConductivityDomainBcWrapper();

  /*!
    @brief Destructor
  */
  virtual ~ConductivityDomainBcWrapper();

  /*!
    @brief Dirichlet and Neumann need to know about bcoef
  */
  virtual void setCoefficients();
  
  /*!
    @brief Set the potential functions
  */
  virtual void setPotentials(const Vector<RefCountedPtr<BaseBCFuncEval> >& a_potentials);

  /*!
    @brief Set bc types. This should come through the factory. 
  */
  virtual void setWallBc(const Vector<RefCountedPtr<WallBc> >& a_wallbc);

  /*!
    @brief Set Robin coefficients
  */
  virtual void setRobinCoefficients(const Vector<RefCountedPtr<RobinCoefficients> >& a_robinco);  
  
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
			   const bool&           a_useHomogeneous);

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
			   const bool&           a_useHomogeneous);
  
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
			      const bool&           a_useHomogeneous);
  
  /*!
    @brief Ghost filling. This is for NWO.
  */
  virtual void fillPhiGhost(FArrayBox&     a_phi,
			    const Box&     a_valid,
			    const Box&     a_domain,
			    Real           a_dx,
			    bool           a_homogeneous);

protected:

  /*!
    @brief Individual stuff
  */
  Vector<RefCountedPtr<ConductivityBaseDomainBC> > m_bc;

  /*!
    @brief Potentials
  */
  Vector<RefCountedPtr<BaseBCFuncEval> > m_potentials;
  Vector<RefCountedPtr<RobinCoefficients> > m_robinCoefficients;         // Robin coefficients for Robin

  /*!
    @brief Is defined or not
  */
  bool m_defined;

};
#include "CD_NamespaceFooter.H"
#endif
