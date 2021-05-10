/*!
  @file   CD_ConductivityElectrostaticDomainBc.H
  @brief  Boundary condition class for using EBConductivityOp in field_solver_multigrid
  @author Robert Marskar
  @date   May 2021
*/

#ifndef _CD_ConductivityElectrostaticDomainBc_H_
#define _CD_ConductivityElectrostaticDomainBc_H_

#include <memory>

#include <BaseDomainBC.H>

#include "CD_ElectrostaticDomainBc.H"
#include "CD_ElectrostaticDomainBcFuncEval.H"
#include "CD_NamespaceHeader.H"

/*!
  @brief Class which interfaces a boundary condition wrapper to ConductivityBaseDomainBC.
  @details This class provides the functionality for setting Dirichlet and Neumann boundary conditions individually on each domain side, with different function. The
  boundary conditions are supplied by reference, so that the user can update the time for which the boundary condition is evaluated (through ElectrostaticDomainBcFuncEval).
  @note This class rolls DirichletConductivityDomainBC and NeumannConductivityDomainBC, and calls the flux functions from those. We have thought about bypassing those
  objects completely by just reimplementing those functinos here. That would make ElectrostaticDomainBcFuncEval pretty redundant and could roll ElectrostaticDomainBc::BcFunc
  directly. 
*/
class ConductivityElectrostaticDomainBc : public ConductivityBaseDomainBC {
public:

  /*!
    @brief Constructor. This does nothing. 
  */
  ConductivityElectrostaticDomainBc() = delete;

  /*!
    @brief Full constructor
  */
  ConductivityElectrostaticDomainBc(const ElectrostaticDomainBc& a_domainBc, const RealVect a_probLo);

  /*!
    @brief Destructor. This does nothing.
  */
  ~ConductivityElectrostaticDomainBc();
  
  /*! 
    @brief Get the face flux on a domain side and direction. 
    @details This is used for getting the conventional face flux. 
    @param[in] a_faceFlux       Regular face flux.
    @param[in] a_phi            Phi.
    @param[in] a_probLo         Lower-left corner of computational domain. 
    @param[in] a_dx             Resolution.
    @param[in] a_idir           Coordinate direction.
    @param[in] a_side           Lo/Hi side in coordinate direction.
    @param[in] a_dit            DataIndex for current box being evaluated. 
    @param[in] a_time           Time. 
    @param[in] a_useHomogeneous Use homogeneous BC or not. 
    @note Because of EBConductivityOp a_time and a_probLo are both zero, and these are instead passed in through the actual BC function objects. 
  */
  void getFaceFlux(BaseFab<Real>&        a_faceFlux,
		   const BaseFab<Real>&  a_phi,
		   const RealVect&       a_probLo,
		   const RealVect&       a_dx,
		   const int&            a_idir,
		   const Side::LoHiSide& a_side,
		   const DataIndex&      a_dit,
		   const Real&           a_time,
		   const bool&           a_useHomogeneous) override final;

  /*! 
    @brief Get the face flux on a domain side and direction. 
    @details This is used for getting the domain flux in irregular cells that are also cut-cells. 
    @param[in] a_faceFlux       Flux
    @param[in] a_vof            VolIndex being evaluated. 
    @param[in] a_comp           Component
    @param[in] a_probLo         Lower-left corner of computational domain. 
    @param[in] a_dx             Grid resolution.
    @param[in] a_idir           Coordinate direction.
    @param[in] a_side           Lo or Hi side in coordinate direction.
    @param[in] a_dit            Box currently being evaluated. 
    @param[in] a_time           Time
    @param[in] a_useHomogeneous Use homogeneous BC or not. 
    @note Because of EBConductivityOp a_time and a_probLo are both zero, and these are instead passed in through the actual BC function objects. 
  */
  void getFaceFlux(Real&                 a_faceFlux,
		   const VolIndex&       a_vof, 
		   const int&            a_comp, 
		   const EBCellFAB&      a_phi, 
		   const RealVect&       a_probLo, 
		   const RealVect&       a_dx, 
		   const int&            a_idir, 
		   const Side::LoHiSide& a_side, 
		   const DataIndex&      a_dit, 
		   const Real&           a_time, 
		   const bool&           a_useHomogeneous) override final;
  
  /*!
    @brief This function is not documented because it is not used in EBConductivityOp. It issues an error if it is called. 
  */
  void getFaceGradPhi(Real&                 a_faceFlux,
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
		      const bool&           a_useHomogeneous) override final;
  
  /*!
    @brief This function is not documented because it is not used in EBConductivityOp. It issues an error if it is called. 
  */
  void fillPhiGhost(FArrayBox&     a_phi,
		    const Box&     a_valid,
		    const Box&     a_domain,
		    Real           a_dx,
		    bool           a_homogeneous) override final;

  
protected:
  
  /*!
    @brief Flag for checking if coefficients have been set. 
  */
  bool m_hasCoeff;

  /*!
    @brief Boundary condition objects.
    @details Constructed during setBcFunctions.
  */
  std::map<ElectrostaticDomainBc::Wall, std::shared_ptr<ConductivityBaseDomainBC> > m_conductivityBaseDomainBcObjects;

  /*!
    @brief Set coefficients. This is simply done to ensure that the boundary condition object has access to the b-coefficient. 
    @details In EBConductivityOp, the boundary conditions need to know about the b-coefficient when BCs are supplied (for evaluating the flux). This function
    provides that functionality. 
  */
  void setCoefficients();

};

#include "CD_NamespaceFooter.H"

#endif
