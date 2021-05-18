/*!
  @file   robinconductivitydomainbc.H
  @brief  Declaration of a Robin boundary condition class for domain faces
  @author Robert Marskar
  @date   Jan. 2018
*/

#ifndef _ROBINCONDUCTIVITYDOMAINBC_
#define _ROBINCONDUCTIVITYDOMAINBC_

#include "robin_coef.H"

#include <LevelData.H>
#include <BaseDomainBC.H>
#include <ConductivityBaseDomainBC.H>

#include "CD_NamespaceHeader.H"
  
/*!
  @brief Boundary condition class for Robin type boundary conditions
*/
class robinconductivitydomainbc : public ConductivityBaseDomainBC {
public:

  /*!
    @brief Constructor
  */
  robinconductivitydomainbc();

  /*!
    @brief Destructor
  */
  ~robinconductivitydomainbc();

  /*!
    @brief Function based coefficients and right hand sides
  */
  virtual void setCoefficientss(RefCountedPtr<robin_coef> a_robinco);

  /*!
    @brief Set constant coefficients
  */
  virtual void setCoefficientss(const Real a_aco, const Real b_bco, const Real a_rhs);

  /*!
    @brief Set data-based coefficients
  */
  virtual void setCoefficientss(const RefCountedPtr<LevelData<EBFluxFAB> >& a_aco,
			 const RefCountedPtr<LevelData<EBFluxFAB> >& a_bco,
			 const RefCountedPtr<LevelData<EBFluxFAB> >& a_rhs);

  /*!
    @brief Face flux
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
    @brief Compute face flux
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
    @brief Compute gradient
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
  
protected:
  
  RefCountedPtr<robin_coef> m_robinco;

  Real m_aCoefficient;

  Real m_bco;
  
  Real m_rhs;

  RefCountedPtr<LevelData<EBFluxFAB> > m_aCoefficientdata;

  RefCountedPtr<LevelData<EBFluxFAB> > m_bcodata;
  
  RefCountedPtr<LevelData<EBFluxFAB> > m_rhsdata;

  bool m_const_coeff;

  bool m_data_coeff;

  bool m_func_coeff;  
};
#include "CD_NamespaceFooter.H"
#endif
