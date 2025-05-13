/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_EBHelmholtzDirichletEBBC.H
  @brief  Declaration of a boundary condition class for EBHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_EBHelmholtzDirichletEBBC_H
#define CD_EBHelmholtzDirichletEBBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_VofUtils.H>
#include <CD_EBHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing Dirichlet EB boundary conditions into EBHelmholtzOp. 
*/
class EBHelmholtzDirichletEBBC : public EBHelmholtzEBBC
{
public:
  /*!
    @brief Default constructor. 
  */
  EBHelmholtzDirichletEBBC();

  /*!
    @brief Destructor
  */
  virtual ~EBHelmholtzDirichletEBBC();

  /*!
    @brief Set BC order. 
    @param[in] a_order Order
    @note The order is also the attempted radius of the stencil. 
  */
  virtual void
  setOrder(const int a_order);

  /*!
    @brief Set stencil weights
    @param[in] a_weight Weighting factor for least squares stencils
    @note a_weight=0 can lead to ill-conditioned systems.
  */
  virtual void
  setWeight(const int a_weight);

  /*!
    @brief Set the value on the EB. 
    @param[in] a_value Constant value on EB
  */
  virtual void
  setValue(const int a_value);

  /*!
    @brief Set the value on the EB. 
    @param[in] a_value Spatially dependent value on EB
  */
  virtual void
  setValue(const std::function<Real(const RealVect& a_pos)>& a_value);

  /*!
    @brief Drop BC order if domain size is equal or below this.
  */
  virtual void
  setDomainDropOrder(const int a_domainSize);

  /*!
    @brief Signal that we should drop the stencil order if the cell is a coarsened cell
    @param[in] a_dropOrder
  */
  virtual void
  setCoarseGridDropOrder(const bool a_dropOrder);

  /*!
    @brief Apply the EB flux. This is the version that is called by EBHelmholtzOp.
    @param[inout] a_vofit             Iterator for cut-cells
    @param[in]    a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_Bcoef             EB-centered B-coefficient
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void
  applyEBFlux(VoFIterator&           a_vofit,
              EBCellFAB&             a_Lphi,
              const EBCellFAB&       a_phi,
              const BaseIVFAB<Real>& a_Bcoef,
              const DataIndex&       a_dit,
              const Real&            a_beta,
              const bool&            a_homogeneousPhysBC) const override;

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
    @brief Stencil approximation order,
  */
  int m_order;

  /*!
    @brief Distance weighting factor for stencils.
  */
  int m_weight;

  /*!
    @brief Special flag for dropping stencil order when domains become coarser than this.
  */
  int m_domainDropOrder;

  /*!
    @brief Special flag for dropping stencil order if the cell is a coarsened cell
  */
  bool m_dropOrder;

  /*!
    @brief Value on EB
  */
  Real m_constantValue;

  /*!
    @brief Value on the EB. 
  */
  std::function<Real(const RealVect& a_pos)> m_functionValue;

  /*!
    @brief Boundary weights
  */
  LayoutData<BaseIVFAB<Real>> m_boundaryWeights;

  /*!
    @brief User define function. 
  */
  virtual void
  define() override;

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
  virtual bool
  getLeastSquaresStencil(std::pair<Real, VoFStencil>& a_stencil,
                         const VolIndex&              a_vof,
                         const VofUtils::Neighborhood a_neighborhood,
                         const DataIndex&             a_dit,
                         const int                    a_order) const;

  /*!
    @brief Get johansen stencil with specified order. Not used but exists for backward compatibility testing and debugging. 
    @param[out] a_stencil Weight and stencil
    @param[out] a_vof     Grid cell
    @param[in]  a_dit     Data index
    @param[in]  a_order   Attempted order for stencil. This is also the stencil radius. 
    @return If this returns true, a_stencil contains a valid stencil and weight. 
  */
  virtual bool
  getJohansenStencil(std::pair<Real, VoFStencil>& a_stencil,
                     const VolIndex&              a_vof,
                     const DataIndex&             a_dit,
                     const int                    a_order) const;

  /*!
    @brief Get Chombo's least squares gradient stencil. Not used but kept for backward compatibility testing and debugging. 
    @param[out] a_stencil Weight and stencil
    @param[out] a_vof     Grid cell
    @param[in]  a_dit     Data index
    @param[in]  a_order   Attempted order for stencil. This is also the stencil radius. 
    @return If this returns true, a_stencil contains a valid stencil and weight. 
  */
  virtual bool
  getChomboLsqStencil(std::pair<Real, VoFStencil>& a_stencil, const VolIndex& a_vof, const DataIndex& a_dit) const;
};

#include <CD_NamespaceFooter.H>

#endif
