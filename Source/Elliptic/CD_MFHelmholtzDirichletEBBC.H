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
#include <CD_EBHelmholtzDirichletEBBC.H>
#include <CD_JumpBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making Dirichlet boundary conditions on the EB in a multifluid context. 
*/
class MFHelmholtzDirichletEBBC : public EBHelmholtzEBBC {
public:

  /*!
    @brief Disallowed weak constructor. Factory must call the other version. 
  */
  MFHelmholtzDirichletEBBC() = delete;

  /*!
    @brief Full constructor
  */
  MFHelmholtzDirichletEBBC(const int a_phase, const RefCountedPtr<JumpBC>& a_jumpBC);

  /*!
    @brief Default constructor
  */
  ~MFHelmholtzDirichletEBBC();

  /*!
    @brief Set BC order. 
    @param[in] a_order Order
    @note The order is also the attempted radius of the stencil. 
  */
  virtual void setOrder(const int a_order);

  /*!
    @brief Set stencil weights
    @param[in] a_weight Weighting factor for least squares stencils
    @note a_weight=0 can lead to ill-conditioned systems.
  */
  virtual void setWeight(const int a_weight);

  /*!
    @brief Set the value on the EB. 
    @param[in] a_value Constant value on EB
  */
  virtual void setValue(const int a_value);

  /*!
    @brief Set the value on the EB. 
    @param[in] a_value Spatially dependent value on EB
  */
  virtual void setValue(const std::function<Real(const RealVect& a_pos)>& a_value);
  
  /*!
    @brief Apply the EB flux. This is the version that is called by EBHelmholtzOp.
    @param[inout] a_vofit             Iterator for cut-cells
    @param[inout] a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void applyEBFlux(VoFIterator&       a_vofit,
			   EBCellFAB&         a_Lphi,
			   const EBCellFAB&   a_phi,
			   const DataIndex&   a_dit,
			   const Real&        a_beta,
			   const bool&        a_homogeneousPhysBC) const override;

protected:

  /*!
    @brief Jump bc object
  */
  RefCountedPtr<JumpBC> m_jumpBC;

  /*!
    @brief Phase
  */
  int m_phase;

    /*!
    @brief Use constant for BC
  */
  bool m_useConstant;

  /*!
    @brief Use function for BC value
  */
  bool m_useFunction;

  /*!
    @brief Stencil approximation order,
  */
  int m_order;

  /*!
    @brief Distance weighting factor for stencils.
  */
  int m_weight;

  /*!
    @brief Value on EB
  */
  Real m_constantValue;

  /*!
    @brief Value on the EB. Does not apply to matching cells. 
  */
  std::function<Real(const RealVect& a_pos) > m_functionValue;

  /*!
    @brief Boundary weights for stencil. 
  */
  LayoutData<BaseIVFAB<Real> > m_boundaryWeights;

  /*!
    @brief Boundary weights for stencil. 
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_gradStencils;

    /*!
    @brief User define function. 
  */
  virtual void define() override;

  /*!
    @brief Get stencil with specified order. 
    @param[out] a_stencil      Weight and stencil
    @param[out] a_vof          Grid cell
    @param[out] a_neighborhood Neighboorhood specifier
    @param[in]  a_dit          Data index
    @param[in]  a_order        Attempted order for stencil. This is also the stencil radius. 
    @return If this returns true, a_stencil contains a valid stencil and weight. 
    @note Stencils are computed using weighted least squares with squared inverse distance weighting. 
  */
  virtual bool getLeastSquaresBoundaryGradStencil(std::pair<Real, VoFStencil>& a_stencil,
						  const VolIndex&              a_vof,
						  const VofUtils::Neighborhood a_neighborhood,
						  const DataIndex&             a_dit,
						  const int                    a_order) const;
};

#include <CD_NamespaceFooter.H>

#endif
