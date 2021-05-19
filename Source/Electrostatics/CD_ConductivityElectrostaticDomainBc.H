/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_ConductivityElectrostaticDomainBc.H
  @brief  Boundary condition class for using EBConductivityOp in FieldSolverMultigrid
  @author Robert Marskar
*/

#ifndef CD_ConductivityElectrostaticDomainBc_H
#define CD_ConductivityElectrostaticDomainBc_H

// Std includes
#include <memory>

// Chombo includes
#include <BaseDomainBC.H>

// Our includes
#include <CD_ElectrostaticDomainBc.H>
#include <CD_NamespaceHeader.H>

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

private:

  /*!
    @brief Class which interfaces from more general boundary condition specification to something usable by the current Chombo multigrid solvers 
    that FieldSolverMultigrid uses. 
  */
  class ElectrostaticDomainBcFuncEval : public BaseBCFuncEval {
  public:
  
    /*!
      @brief Deleted initializing constructor. 
    */
    ElectrostaticDomainBcFuncEval() = delete;

    /*!
      @brief Full constructor. 
      @param[in] a_bcFunc The boundary condition function.
      @param[in] a_probLo Lower-left corner of computational domain. 
      @note The origin is sent in because EBConductivityOp calls the boundary condition function with no origin. So, this evaluation of BaseBCFuncEval rectifies that.
    */
    ElectrostaticDomainBcFuncEval(const ElectrostaticDomainBc::BcFunction a_bcFunc, const RealVect a_probLo);

    /*!
      @brief Destructor. This does nothing. 
    */
    ~ElectrostaticDomainBcFuncEval();

    /*!
      @brief Return functional evaluation in space and time of the boundary condition.
      @param[in] a_point Point in space
      @param[in] a_comp  Component. 
      @details The time-dependence is implicit - the user needs to set the time by calling setTime() first. 
      @note Calling this without having set the function first yields an error. 
    */
    Real value(const RealVect& a_point, const int& a_comp) const override final;

    /*!
      @brief Implementation of BaseBCFuncEval.
      @param[in] a_point    Position in space
      @param[in] a_comp     Component in elliptic equation.
      @param[in] a_derivDir Derivative direction.
      @return Throws an error if called since this function is not used in FieldSolverMultigrid code. 
    */
    Real derivative(const RealVect& a_point, const int& a_comp, const int& a_derivDir) const override final;

  protected:
  
    /*!
      @brief Lower left corner of computational boundary. 
    */
    RealVect m_probLo;

    /*!
      @brief Space- and time-dependent boundary condition function. 
      @note This is just an alias. 
    */
    ElectrostaticDomainBc::BcFunction m_bcFunc;
  };
};

#include <CD_NamespaceFooter.H>

#endif
