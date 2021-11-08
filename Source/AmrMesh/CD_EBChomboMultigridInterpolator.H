/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBChomboMultigridInterpolator.H
  @brief  Declaration of an interpolator class which uses the Chombo EB interpolator for multigrid. 
  @author Robert Marskar
*/

#ifndef CD_EBChomboMultigridInterpolator_H
#define CD_EBChomboMultigridInterpolator_H

// Chombo includes
#include <EBQuadCFInterp.H>

// Our includes
#include <CD_EBMultigridInterpolator.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Interpolator class which uses the Chombo EB interpolator for multigrid. 
  @details This class fills only a single layer of ghost cells. The user must make sure that the cut-cell stencils only has a radius of one if there's an EBCF situation. If
  this class is used together with EBHelmholtzOp, EBHelmholtzOp will call getGhostCF which signals how many ghost cells were filled and automatically adjust the stencils
  from that. 
*/
class EBChomboMultigridInterpolator : public EBMultigridInterpolator {
public:

  /*!
    @brief Disallowed constructor, use the full constructor. 
  */
  EBChomboMultigridInterpolator() = delete;

  /*!
    @brief Disallowed copy constructor.
  */
  EBChomboMultigridInterpolator(const EBChomboMultigridInterpolator& a_other) = delete;

  /*!
    @brief Full constructor
    @param[in] a_eblgFine The fine level
    @param[in] a_eblgCoar The coarse level
    @param[in] a_refRat   Refinement ratio between levels
    @param[in] a_nVar     Number of variables
  */
  EBChomboMultigridInterpolator(const EBLevelGrid& a_eblgFine,
				const EBLevelGrid& a_eblgCoar,
				const int          a_refRat,
				const int          a_nVar);

  /*!
    @brief Destructor (does nothing)
  */
  virtual ~EBChomboMultigridInterpolator();

  /*!
    @brief Return number of ghost cells across CF interface.
  */
  virtual int getGhostCF() const override;

  /*!
    @brief Do inhomogeneous interpolation
    @param[inout] a_phiFine   Fine phi
    @param[in]    a_phiCoar   Coarse phi
    @param[in]    a_variables Variables to interpolate
  */
  virtual void coarseFineInterp(LevelData<EBCellFAB>& a_phiFine, const LevelData<EBCellFAB>& a_phiCoar, const Interval a_variables) override;

  /*!
    @brief Do homogeneous interpolation
    @param[in] a_phiFine   Fine phi
    @param[in] a_variables Variables to interpolate
  */
  virtual void coarseFineInterpH(LevelData<EBCellFAB>& a_phiFine, const Interval a_variables) const override;

  /*!
    @brief Do homogeneous coarse-fine interpolation on a grid patch
    @param[in] a_phiFine   Fine-level data
    @param[in] a_variables Variables to be interpolated
    @param[in] a_dit       Data index
    @note This is the kernel version. The constructor guarantees that the fine-grid stencil will not out of the ghost region defined by a_ghostVector. However, this is
    the direct kernel version which is public because multigrid operators might need it for performance reasons. This routine does NOT warn against input data that do not
    contain enough ghost cells. 
  */
  virtual void coarseFineInterpH(EBCellFAB& a_phiFine, const Interval a_variables, const DataIndex& a_dit) const override;

protected:

  /*!
    @brief For EB-aware multigrid interpolation of scalar quantities
  */
  RefCountedPtr<EBQuadCFInterp> m_scalarInterpolator;
};

#include <CD_NamespaceFooter.H>

#endif
