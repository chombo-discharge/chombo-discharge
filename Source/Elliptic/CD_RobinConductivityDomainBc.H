/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_RobinConductivityDomainBc.H
  @brief  Declaration of a Robin boundary condition class for domain faces
  @author Robert Marskar
*/

#ifndef CD_RobinConductivityDomainBc_H
#define CD_RobinConductivityDomainBc_H

// Chombo includes
#include <LevelData.H>
#include <BaseDomainBC.H>
#include <ConductivityBaseDomainBC.H>

// Our includes
#include <CD_RobinCoefficients.H>
#include <CD_NamespaceHeader.H>
  
/*!
  @brief Boundary condition class for Robin type boundary conditions
*/
class RobinConductivityDomainBc : public ConductivityBaseDomainBC {
public:

  /*!
    @brief Constructor
  */
  RobinConductivityDomainBc();

  /*!
    @brief Destructor
  */
  ~RobinConductivityDomainBc();

  /*!
    @brief Function based coefficients and right hand sides
  */
  virtual void setCoefficients(RefCountedPtr<RobinCoefficients> a_robinco);

  /*!
    @brief Set constant coefficients
  */
  virtual void setCoefficients(const Real a_aco, const Real b_bco, const Real a_rhs);

  /*!
    @brief Set data-based coefficients
  */
  virtual void setCoefficients(const RefCountedPtr<LevelData<EBFluxFAB> >& a_aco,
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
  
  RefCountedPtr<RobinCoefficients> m_robinCoefficients;

  Real m_aco;

  Real m_bco;
  
  Real m_rhs;

  RefCountedPtr<LevelData<EBFluxFAB> > m_acodata;

  RefCountedPtr<LevelData<EBFluxFAB> > m_bcodata;
  
  RefCountedPtr<LevelData<EBFluxFAB> > m_rhsdata;

  bool m_const_coeff;

  bool m_data_coeff;

  bool m_func_coeff;  
};

#include <CD_NamespaceFooter.H>

#endif
