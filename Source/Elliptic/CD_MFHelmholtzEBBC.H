/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_MFHelmholtzEBBC.H
  @brief  Declaration of an EB boundary condition class for MFHelmholtzOp.
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzEBBC_H
#define CD_MFHelmholtzEBBC_H

// Our includes
#include <CD_EBHelmholtzEBBC.H>
#include <CD_MFHelmholtzJumpBC.H>
#include <CD_VofUtils.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for making boundary conditions on the EB in a multifluid context. 
  @details This class is almost like EBHelmholtzEBBC, except that it enforces matching on multiphase cells. The user can inherit from this
  layer if he needs additional BCs in a multiphase context. 
*/
class MFHelmholtzEBBC : public EBHelmholtzEBBC
{
public:
  /*!
    @brief Disallowed weak constructor. Factory must call the other version. 
  */
  MFHelmholtzEBBC() = delete;

  /*!
    @brief Disallowed copy constructor. Factory must call the other version. 
  */
  MFHelmholtzEBBC(const MFHelmholtzEBBC& a_other) = delete;

  /*!
    @brief Full constructor
  */
  MFHelmholtzEBBC(const int a_phase, const RefCountedPtr<MFHelmholtzJumpBC>& a_jumpBC);

  /*!
    @brief Default constructor
  */
  virtual ~MFHelmholtzEBBC();

  /*!
    @brief Apply the EB flux. This is the version that is called by EBHelmholtzOp.
    @details Users get a new interface to the subset of cells that are not multi-phase cells. 
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
    @brief Jump bc object
  */
  RefCountedPtr<MFHelmholtzJumpBC> m_jumpBC;

  /*!
    @brief Phase
  */
  int m_phase;

  /*!
    @brief Weights for matching cells.
  */
  LayoutData<BaseIVFAB<Real>> m_boundaryWeights;

  /*!
    @brief User define function. 
  */
  void
  define() override final;

  /*!
    @brief Define stuff for multi-phase
  */
  void
  defineMultiPhase();

  /*!
    @brief Define stuff for multi-phase
  */
  virtual void
  defineSinglePhase() = 0;

  /*!
    @brief Apply flux single phase
    @param[inout] a_singlePhaseVofs   Iterator for single-phase cut-cells
    @param[inout] a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_Bcoef             EB-centered B-coefficient
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void
  applyEBFluxSinglePhase(VoFIterator&           a_singlePhaseVofs,
                         EBCellFAB&             a_Lphi,
                         const EBCellFAB&       a_phi,
                         const BaseIVFAB<Real>& a_Bcoef,
                         const DataIndex&       a_dit,
                         const Real&            a_beta,
                         const bool&            a_homogeneousPhysBC) const = 0;

  /*!
    @brief Apply flux in multiphase cells.
    @param[inout] a_multiPhaseVofs    Iterator for multi-phase cut-cells
    @param[inout] a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_Bcoef             EB-centered B-coefficient
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  void
  applyEBFluxMultiPhase(VoFIterator&           a_multiPhaseVofs,
                        EBCellFAB&             a_Lphi,
                        const EBCellFAB&       a_phi,
                        const BaseIVFAB<Real>& a_Bcoef,
                        const DataIndex&       a_dit,
                        const Real&            a_beta,
                        const bool&            a_homogeneousPhysBC) const;

  /*!
    @brief Get stencil with specified order. 
    @param[out] a_stencil      Weight and stencil
    @param[out] a_vof          Grid cell
    @param[out] a_neighborhood Neighboorhood specifier
    @param[in]  a_dit          Data index
    @param[in]  a_order        Attempted order for stencil. This is also the stencil radius. 
    @param[in]  a_weight       Weighting factor for weighted least squares
    @return If this returns true, a_stencil contains a valid stencil and weight. 
    @note Stencils are computed using weighted least squares with squared inverse distance weighting. 
  */
  virtual bool
  getLeastSquaresBoundaryGradStencil(std::pair<Real, VoFStencil>& a_stencil,
                                     const VolIndex&              a_vof,
                                     const VofUtils::Neighborhood a_neighborhood,
                                     const DataIndex&             a_dit,
                                     const int                    a_order,
                                     const int                    a_weight) const;
};

#include <CD_NamespaceFooter.H>

#endif
