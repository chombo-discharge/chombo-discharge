/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzDirichletEBBC.H
  @brief  Declaration of a Dirichlet boundary condition class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzDirichletEBBC_H
#define CD_MFHelmholtzDirichletEBBC_H

// Our includes
#include <CD_MFHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making Dirichlet boundary conditions on the EB in a multifluid context. 
*/
class MFHelmholtzDirichletEBBC : public MFHelmholtzEBBC
{
public:
  /*!
    @brief Constructor. Must subsequently set value, order, and weight. 
    @param[in] a_phase  Phase
    @param[in] a_jumpBC Jump BC object
  */
  MFHelmholtzDirichletEBBC(const int a_phase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC);

  /*!
    @brief Default constructor
  */
  virtual ~MFHelmholtzDirichletEBBC();

  /*!
    @brief Set the value on the EB. 
    @param[in] a_value Constant value on EB
  */
  virtual void
  setValue(const Real a_value);

  /*!
    @brief Set the value on the EB. 
    @param[in] a_value Spatially dependent value on EB
  */
  virtual void
  setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Set stencil order
    @param[in] a_order Stencil order. Must be > 0
  */
  void
  setOrder(const int a_order);

  /*!
    @brief Set stencil weight
    @param[in] a_order Stencil weight. Must be >= 0
  */
  void
  setWeight(const int a_weight);

  /*!
    @brief Drop BC order if domain size is equal or below this.
  */
  virtual
  void setDomainDropOrder(const int a_domainSize);    

protected:
  /*!
    @brief Use constant for BC
  */
  bool m_useConstant;

  /*!
    @brief Use function for BC value
  */
  bool m_useFunction;

  /*!
    @brief Value on EB
  */
  Real m_constantValue;

  /*!
    @brief Stencil order
  */
  int m_order;

  /*!
    @brief Stencil weight
  */
  int m_weight;

  /*!
    @brief Special flag for dropping stencil order when domains become coarser than this.
  */
  int m_domainDropOrder;    

  /*!
    @brief Value on the EB. Does not apply to matching cells. 
  */
  std::function<Real(const RealVect& a_pos)> m_functionValue;

  /*!
    @brief User define function. 
  */
  virtual void
  defineSinglePhase() override;

  /*!
    @brief Apply flux single phase
    @param[inout] a_singlePhaseVofs   Iterator for single-phase cut-cells
    @param[inout] a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void
  applyEBFluxSinglePhase(VoFIterator&     a_singlePhaseVofs,
                         EBCellFAB&       a_Lphi,
                         const EBCellFAB& a_phi,
                         const DataIndex& a_dit,
                         const Real&      a_beta,
                         const bool&      a_homogeneousPhysBC) const override;
};

#include <CD_NamespaceFooter.H>

#endif
