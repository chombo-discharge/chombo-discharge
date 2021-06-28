/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*
  @file   CD_MFHelmholtzRobinEBBC.H
  @brief  Declaration of a boundary condition class for MFHelmholtzOp
  @author Robert Marskar
*/

#ifndef CD_MFHelmholtzRobinEBBC_H
#define CD_MFHelmholtzRobinEBBC_H

// Std includes
#include <functional>

// Our includes
#include <CD_MFHelmholtzEBBC.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Base class for passing Robin MF boundary conditions into MFHelmholtzOp. 
  @note MFHelmholtzRobinEBBC ignores the a_order and a_weight argument passed into MFHelmholtzEBBC. It computes the interpolation of phi onto the boundary
  using a first order accurate stencil (with range 1). 
*/
class MFHelmholtzRobinEBBC : public MFHelmholtzEBBC {
public:

  /*!
    @brief Full constructor
  */
  MFHelmholtzRobinEBBC(const int a_phase, const RefCountedPtr<JumpBC>& a_jumpBC);

  /*!
    @brief Destructor
  */
  ~MFHelmholtzRobinEBBC();

  /*!
    @brief Set constant coefficients
    @param[in] a_A A-coefficient
    @param[in] a_B B-coefficient
    @param[in] a_C C-coefficient
  */
  void setCoefficients(const Real a_A, const Real a_B, const Real a_C);

  /*!
    @brief Set variable coefficients
    @param[in] a_A A-coefficient
    @param[in] a_B B-coefficient
    @param[in] a_C C-coefficient
  */
  void setCoefficients(const std::function<Real(const RealVect& a_pos) >& a_A,
		       const std::function<Real(const RealVect& a_pos) >& a_B,
		       const std::function<Real(const RealVect& a_pos) >& a_C);


protected:

  /*
    @brief Use constant for BC
  */
  bool m_useConstant;

  /*!
    @brief Use function for BC value
  */
  bool m_useFunction;

  /*!
    @brief Constant A-coefficient
  */
  Real m_constantA;

  /*!
    @brief Constant B-coefficient
  */
  Real m_constantB;

  /*!
    @brief Constant C-coefficient
  */
  Real m_constantC;

  /*!
    @brief Function-based A-coefficient
  */
  std::function<Real(const RealVect& a_pos)> m_functionA;

  /*!
    @brief Function-based B-coefficient
  */
  std::function<Real(const RealVect& a_pos)> m_functionB;

  /*!
    @brief Function-based C-coefficient
  */
  std::function<Real(const RealVect& a_pos)> m_functionC;

  /*!
    @brief Interpolation stencils
  */
  LayoutData<BaseIVFAB<VoFStencil> > m_interpolationStencils;

  /*!
    @brief User define function. 
  */
  virtual void defineSinglePhase() override;

  /*!
    @brief Apply flux single phase
    @param[inout] a_singlePhaseVofs   Iterator for single-phase cut-cells
    @param[inout] a_Lphi              Operator kappa*L(phi)
    @param[in]    a_phi               Phi, cell-centered. 
    @param[in]    a_dit               Data index    
    @param[in]    a_beta              Beta. Passed in from operator. 
    @param[in]    a_homogeneousPhysBC Homogeneous BC or not. 
  */
  virtual void applyEBFluxSinglePhase(VoFIterator&       a_singlePhaseVofs,
				      EBCellFAB&         a_Lphi,
				      const EBCellFAB&   a_phi,
				      const DataIndex&   a_dit,
				      const Real&        a_beta,
				      const bool&        a_homogeneousPhysBC) const override;

  /*!
    @brief Get an interpolation stencil 
    @param[in] a_vof Grid cell
    @param[in] a_dit Data index. 
    @note This extrapolates from cell-centers to the EB using a first order expression. 
  */
  VoFStencil getInterpolationStencil(const VolIndex& a_vof, const DataIndex& a_dit) const;

  /*!
    @brief Get an interpolation stencil using all vofs in a monotone path (avoids reaching over EBs)
    @param[in] a_vof Grid cell
    @param[in] a_dit Data index. 
  */
  VoFStencil getMonoPathStencil(const VolIndex& a_vof, const DataIndex& a_dit) const;

  /*!
    @brief Get an interpolation stencil using Taylor extrapolation. 
    @param[in] a_vof Grid cell
    @param[in] a_dit Data index. 
  */
  VoFStencil getTaylorStencil(const VolIndex& a_vof, const DataIndex& a_dit) const;

};

#include <CD_NamespaceFooter.H>

#endif
