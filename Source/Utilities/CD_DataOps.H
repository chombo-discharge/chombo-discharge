/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_DataOps.H
  @brief  Agglomeration of useful data operations
  @author Robert Marskar
*/

#ifndef CD_DATAOPS_H
#define CD_DATAOPS_H

// Std includes
#include <array>

// Chombo includes
#include <FaceIterator.H>
#include <VoFIterator.H>

// Our includes
#include <CD_Average.H>
#include <CD_EBAMRData.H>
#include <CD_Decorations.H>
#include <CD_MFInterfaceFAB.H>
#include <CD_NamespaceHeader.H>

/**
  @brief Agglomeration of useful data operations.
  @note All methods in this class assumes that data are allocated over the same realm. If you are trying to do something across realms, you need to copy
  the data onto realms that lie on the same data holders first!.
*/
class DataOps
{
public:
  /**
    @brief Set value in an MFInterfaceFAB data holder
    @param[in] a_lhs   Data holder
    @param[in] a_value Value to set
  */
  template <typename T>
  static void
  setValue(LevelData<MFInterfaceFAB<T>>& a_lhs, const T& a_value);

  /**
    @brief Sign function. Returns +/- if the value is > 0 or < 0
    @param[in] a_value Value to evaluate. 
    @return Return value
  */
  template <typename T>
  static int
  sgn(const T a_value);

  /**
    @brief Routine which computes the average of a cell-centered quantity on faces for the normal component only.
    @param[out] a_faceData Face data. Must have one component. 
    @param[in]  a_cellData Cell data. Must have SpaceDim components. 
    @param[in]  a_domains  AMR domains.
    @param[in] a_tanGhosts Tan ghosts
  */
  static void
  averageCellVelocityToFaceVelocity(EBAMRFluxData&               a_faceData,
                                    const EBAMRCellData&         a_cellData,
                                    const Vector<ProblemDomain>& a_domains,
                                    const int                    a_tanGhosts);

  /**
    @brief Routine which computes the average of a cell-centered quantity on faces using pre-built face iterators.
    @param[out] a_faceData  Face data. Must have one component.
    @param[in]  a_cellData  Cell data. Must have SpaceDim components.
    @param[in]  a_domains   AMR domains.
    @param[in]  a_tanGhosts Tangential ghost faces.
    @param[in]  a_faceIter  Pre-built tangential-ghost face iterators for each level.
  */
  static void
  averageCellVelocityToFaceVelocity(EBAMRFluxData&                                                                  a_faceData,
                                    const EBAMRCellData&                                                            a_cellData,
                                    const Vector<ProblemDomain>&                                                    a_domains,
                                    const int                                                                       a_tanGhosts,
                                    Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>& a_faceIter);

  /**
    @brief Routine which computes the average of a cell-centered quantity on faces for the normal component only.
    @param[out] a_faceData Face data. Must have one component.
    @param[in]  a_cellData Cell data. Must have SpaceDim components.
    @param[in]  a_domain   Problem domain
    @param[in] a_tanGhosts Tan ghosts
  */
  static void
  averageCellVelocityToFaceVelocity(LevelData<EBFluxFAB>&       a_faceData,
                                    const LevelData<EBCellFAB>& a_cellData,
                                    const ProblemDomain&        a_domain,
                                    const int                   a_tanGhosts);

  /**
    @brief Routine which computes the average of a cell-centered quantity on faces using pre-built face iterators.
    @param[out] a_faceData  Face data. Must have one component.
    @param[in]  a_cellData  Cell data. Must have SpaceDim components.
    @param[in]  a_domain    Problem domain.
    @param[in]  a_tanGhosts Tangential ghost faces.
    @param[in]  a_faceIter  Pre-built tangential-ghost face iterators.
  */
  static void
  averageCellVelocityToFaceVelocity(LevelData<EBFluxFAB>&                                   a_faceData,
                                    const LevelData<EBCellFAB>&                             a_cellData,
                                    const ProblemDomain&                                    a_domain,
                                    const int                                               a_tanGhosts,
                                    LayoutData<std::array<FaceIterator, SpaceDim>>&         a_faceIter);

  /**
    @brief Average all components of the cell-centered data to faces.
    @details Calls the other version; uses arihtmetic averaging and does not fill ghost faces. 
    @param[out] a_faceData Face data. 
    @param[in]  a_cellData Cell data. 
    @param[in]  a_domains  AMR domains.
    @note This does not fill "tangential" ghost faces.
  */
  static void
  averageCellToFace(EBAMRFluxData& a_faceData, const EBAMRCellData& a_cellData, const Vector<ProblemDomain>& a_domains);

  /**
    @brief Routine which averages a cell-centered component to faces. 
    @details User can choose what type of averaging to use, as well as choosing the number of ghost face to fill. These
    faces are "tangential ghost faces", i.e. if the face direction is along +x, we extend the box with faces in the y- and z-
    directions.
    @param[out] a_faceData     Face data.
    @param[in]  a_cellData     Cell data.
    @param[in]  a_domains      AMR domains.
    @param[in]  a_tanGhosts    Tangential ghost faces to fill.
    @param[in]  a_faceInterval Face variables
    @param[in]  a_cellInterval Cell variables
    @param[in]  a_average      Averaging method
  */
  static void
  averageCellToFace(EBAMRFluxData&               a_faceData,
                    const EBAMRCellData&         a_cellData,
                    const Vector<ProblemDomain>& a_domains,
                    const int                    a_tanGhosts,
                    const Interval&              a_faceInterval,
                    const Interval&              a_cellInterval,
                    const Average&               a_average);

  /**
    @brief Routine which averages a cell-centered component to faces using pre-built face iterators.
    @param[out] a_faceData     Face data.
    @param[in]  a_cellData     Cell data.
    @param[in]  a_domains      AMR domains.
    @param[in]  a_tanGhosts    Tangential ghost faces to fill.
    @param[in]  a_faceInterval Face variables.
    @param[in]  a_cellInterval Cell variables.
    @param[in]  a_average      Averaging method.
    @param[in]  a_faceIter     Pre-built tangential-ghost face iterators for each level.
  */
  static void
  averageCellToFace(EBAMRFluxData&                                                                  a_faceData,
                    const EBAMRCellData&                                                            a_cellData,
                    const Vector<ProblemDomain>&                                                    a_domains,
                    const int                                                                       a_tanGhosts,
                    const Interval&                                                                 a_faceInterval,
                    const Interval&                                                                 a_cellInterval,
                    const Average&                                                                  a_average,
                    Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>& a_faceIter);

  /**
    @brief Routine which averages a cell-centered component to faces.
    @details User can choose what type of averaging to use, as well as choosing the number of ghost face to fill. These
    faces are "tangential ghost faces", i.e. if the face direction is along +x, we extend the box with faces in the y- and z-
    directions.
    @param[out] a_faceData     Face data.
    @param[in]  a_cellData     Cell data.
    @param[in]  a_domain       Problem domain.
    @param[in]  a_tanGhosts    Tangential ghost faces to fill.
    @param[in]  a_faceInterval Face variables
    @param[in]  a_cellInterval Cell variables
    @param[in]  a_average      Averaging method
  */
  static void
  averageCellToFace(LevelData<EBFluxFAB>&       a_faceData,
                    const LevelData<EBCellFAB>& a_cellData,
                    const ProblemDomain&        a_domain,
                    const int                   a_tanGhosts,
                    const Interval&             a_faceInterval,
                    const Interval&             a_cellInterval,
                    const Average&              a_average);

  /**
    @brief Routine which averages a cell-centered component to faces using pre-built face iterators.
    @param[out] a_faceData     Face data.
    @param[in]  a_cellData     Cell data.
    @param[in]  a_domain       Problem domain.
    @param[in]  a_tanGhosts    Tangential ghost faces to fill.
    @param[in]  a_faceInterval Face variables.
    @param[in]  a_cellInterval Cell variables.
    @param[in]  a_average      Averaging method.
    @param[in]  a_faceIter     Pre-built tangential-ghost face iterators.
  */
  static void
  averageCellToFace(LevelData<EBFluxFAB>&                                   a_faceData,
                    const LevelData<EBCellFAB>&                             a_cellData,
                    const ProblemDomain&                                    a_domain,
                    const int                                               a_tanGhosts,
                    const Interval&                                         a_faceInterval,
                    const Interval&                                         a_cellInterval,
                    const Average&                                          a_average,
                    LayoutData<std::array<FaceIterator, SpaceDim>>&         a_faceIter);

  /**
    @brief Average all components of face centered data to cell centers.
    @param[out] a_cellData Cell data.
    @param[in]  a_faceData Face data.
    @param[in]  a_domains  AMR domains.
    @param[in]  a_vofIter  Pre-built VoF iterators (all cut-cells) for each level.
    @note Only valid cells are filled (not ghost cells)
  */
  static void
  averageFaceToCell(EBAMRCellData&                                        a_cellData,
                    const EBAMRFluxData&                                  a_faceData,
                    const Vector<ProblemDomain>&                          a_domains,
                    const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Average all components of face centered data to cell centers.
    @param[out] a_cellData Cell data.
    @param[in]  a_fluxData Face data.
    @param[in]  a_domain   Problem domain.
    @param[in]  a_vofIter  Pre-built VoF iterator (all cut-cells).
    @note Only valid cells are filled (not ghost cells)
  */
  static void
  averageFaceToCell(LevelData<EBCellFAB>&       a_cellData,
                    const LevelData<EBFluxFAB>& a_fluxData,
                    const ProblemDomain&        a_domain,
                    LayoutData<VoFIterator>&    a_vofIter);

  /**
    @brief Compute lhs = a * x + b * y
    @param[out] a_lhs Result
    @param[in] a_x First FAB
    @param[in] a_y Second FAB
    @param[in] a_a Scaling factor
    @param[in] a_b Scaling factor
  */
  static void
  axby(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_x,
       const LevelData<EBCellFAB>& a_y,
       const Real                  a_a,
       const Real                  a_b) noexcept;

  /**
    @brief Compute a new value given the old cell value.
    @details Useful for quickly calculating a cell value as e.g. exp(x) where the original data was only x
    @param[in,out] a_data    Input data
    @param[in]    a_func    Function for transforming the data
    @param[in]    a_vofIter Pre-built VoF iterators (multiply-cut cells) for each level.
  */
  static void
  compute(EBAMRCellData&                                        a_data,
          const std::function<Real(const Real a_cellValue)>&    a_func,
          const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter) noexcept;

  /**
    @brief Compute a new value given the old cell value.
    @details Useful for quickly calculating a cell value as e.g. exp(x) where the original data was only x
    @param[in,out] a_data    Input data
    @param[in]    a_func    Function for transforming the data
    @param[in]    a_vofIter Pre-built VoF iterator (multiply-cut cells).
  */
  static void
  compute(LevelData<EBCellFAB>&                              a_data,
          const std::function<Real(const Real a_cellValue)>& a_func,
          LayoutData<VoFIterator>&                           a_vofIter) noexcept;

  /**
    @brief Compote the cell-wise dot product between two data holders. 
    @param[out] a_result Result. Holds the dot product in each cell.
    @param[in]  a_data1 First data holder.
    @param[in]  a_data2 Other data holder.
    @note The input data data1 and data2 must have the same number of input components. Summation is done over all components (does not have to be SpaceDim).
  */
  static void
  dotProduct(MFAMRCellData& a_result, const MFAMRCellData& a_data1, const MFAMRCellData& a_data2);

  /**
    @brief Compote the cell-wise dot product between two data holders. 
    @details This calls the single-phase dotProduct patch versions with EBCellFABs on each phase.
    @param[out] a_result Result. Holds the dot product in each cell.
    @param[in]  a_data1 First data holder.
    @param[in]  a_data2 Other data holder.
    @note The input data data1 and data2 must have the same number of input components. Summation is done over all components (does not have to be SpaceDim).
  */
  static void
  dotProduct(LevelData<MFCellFAB>& a_result, const LevelData<MFCellFAB>& a_data1, const LevelData<MFCellFAB>& a_data2);

  /**
    @brief Compute the cell-wise dot product between two data holders.
    @param[out] a_result  Result. Holds the dot product in each cell.
    @param[in]  a_data1   First data holder.
    @param[in]  a_data2   Other data holder.
    @param[in]  a_vofIter Pre-built VoF iterators (multiply-cut cells) for each level.
    @note The input data data1 and data2 must have the same number of input components. Summation is done over all components (does not have to be SpaceDim).
  */
  static void
  dotProduct(EBAMRCellData&                                        a_result,
             const EBAMRCellData&                                  a_data1,
             const EBAMRCellData&                                  a_data2,
             const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Compute the cell-wise dot product between two data holders.
    @param[out] a_result  Result. Holds the dot product in each cell.
    @param[in]  a_data1   First data holder.
    @param[in]  a_data2   Other data holder.
    @param[in]  a_vofIter Pre-built VoF iterator (multiply-cut cells).
    @note The input data data1 and data2 must have the same number of input components. Summation is done over all components (does not have to be SpaceDim).
  */
  static void
  dotProduct(LevelData<EBCellFAB>&       a_result,
             const LevelData<EBCellFAB>& a_data1,
             const LevelData<EBCellFAB>& a_data2,
             LayoutData<VoFIterator>&    a_vofIter);

  /**
    @brief Compute the cell-wise dot product between two data holders.
    @param[out] a_result  Result. Holds the dot product in each cell.
    @param[in]  a_data1   First data holder.
    @param[in]  a_data2   Other data holder.
    @param[in]  a_box     Cell box.
    @param[in]  a_vofIter Pre-built VoF iterator (multiply-cut cells) for this patch.
    @note The input data data1 and data2 must have the same number of input components. Summation is done over all components (does not have to be SpaceDim).
  */
  static void
  dotProduct(EBCellFAB&       a_result,
             const EBCellFAB& a_data1,
             const EBCellFAB& a_data2,
             const Box&       a_box,
             VoFIterator&     a_vofIter);

  /**
    @brief Apply a convolved filter phi = alpha * phi_i + 0.5*(1-alpha) * [phi_(i+s) + phi_(i-s)] in each direction. 
    @details Ghost cells and coarsening must be done before entering this routine
    @param[in,out] a_data   Data to be filtered.
    @param[in]    a_alpha  Filter factor
    @param[in]    a_stride Cell stride s
    @param[in]    a_zeroEB Set value in covered cells to zero or not
    @note This routine will fill ghost cells outside the domain boundary, but ignore the corner ghost cells at the domain corners. This routine may
    be dangerous to use with EBs as the data inside the EB is either kept as-is or replaced by zero.
  */
  static void
  filterSmooth(EBAMRCellData& a_data, const Real a_alpha, const int a_stride, const bool a_zeroEB) noexcept;

  /**
    @brief Apply a convolved filter phi = alpha * phi_i + 0.5*(1-alpha) * [phi_(i+s) + phi_(i-s)] in each direction. 
    @details Ghost cells must be updated before entering this routine. 
    @param[in,out] a_data   Data to be filtered.
    @param[in]    a_alpha  Filter factor
    @param[in]    a_stride Cell stride s
    @param[in]    a_zeroEB Set value in covered cells to zero or not
    @note This routine will fill ghost cells outside the domain boundary, but ignore the corner ghost cells at the domain corners. This routine may
    be dangerous to use with EBs as the data inside the EB is either kept as-is or replaced by zero.
  */
  static void
  filterSmooth(LevelData<EBCellFAB>& a_data, const Real a_alpha, const int a_stride, const bool a_zeroEB) noexcept;

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(MFAMRCellData& a_lhs, const MFAMRCellData& a_rhs, const Real a_scale) noexcept;

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const Real a_scale) noexcept;

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale) noexcept;

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real& a_scale) noexcept;

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRFluxData& a_lhs, const EBAMRFluxData& a_rhs, const Real& a_scale);

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<EBFluxFAB>& a_lhs, const LevelData<EBFluxFAB>& a_rhs, const Real& a_scale);

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBFluxFAB& a_lhs, const EBFluxFAB& a_rhs, const Real& a_scale);

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs, const Real& a_scale);

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs, const Real& a_scale);

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRIFData& a_lhs, const EBAMRIFData& a_rhs, const Real& a_scale);

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<DomainFluxIFFAB>& a_lhs, const LevelData<DomainFluxIFFAB>& a_rhs, const Real& a_scale);

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @details This versions only adds irregular data from an EBAMRIVFAB. It ignores regular cells. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRCellData& a_lhs, const EBAMRIVData& a_rhs, const Real a_scale);

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @details This versions only adds irregular data from an EBAMRIVFAB. It ignores regular cells. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<EBCellFAB>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs, const Real a_scale);

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @details This versions only adds irregular data form an EBCellFAB to an BaseIVFAB. It ignores regular cells. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRIVData& a_lhs, const EBAMRCellData& a_rhs, const Real a_scale);

  /**
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @details This versions only adds irregular data form an EBCellFAB to an BaseIVFAB. It ignores regular cells. 
    @param[in,out] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale);

  /**
    @brief Copy data from one data holder to another.
    @param[out] a_dst Destination data holder
    @param[in ] a_src Source data holder
  */
  static void
  copy(MFAMRCellData& a_dst, const MFAMRCellData& a_src);

  /**
    @brief Copy data from one data holder to another.
    @param[out] a_dst Destination data holder
    @param[in ] a_src Source data holder
  */
  static void
  copy(EBAMRCellData& a_dst, const EBAMRCellData& a_src);

  /**
    @brief Copy data from one data holder to another.
    @param[out] a_dst Destination data holder
    @param[in ] a_src Source data holder
  */
  static void
  copy(EBAMRIVData& a_dst, const EBAMRIVData& a_src);

  /**
    @brief Copy data from one data holder to another. Allows specification of components. 
    @param[out] a_dst       Destination data holder
    @param[in ] a_src       Source data holder
    @param[out] a_dstInterv Destination data component
    @param[in ] a_srcInterv Source data components
  */
  template <typename T>
  static void
  copy(EBAMRData<T>& a_dst, const EBAMRData<T>& a_src, const Interval a_dstInterv, const Interval a_srcInterv);

  /**
    @brief Divide one data holder by another. Does specified components.
    @param[in,out] a_lhs Data to be divided
    @param[in]    a_rhs Divisor.
    @param[in]    a_lhsComp Destination data component
    @param[in]    a_rhsComp Source data component
  */
  static void
  divide(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_lhsComp, const int a_rhsComp);

  /**
    @brief Divide one data holder by another. Does specified components.
    @param[in,out] a_lhs Data to be divided
    @param[in]    a_rhs Divisor.
    @param[in]    a_lhsComp Destination data component
    @param[in]    a_rhsComp Source data component
  */
  static void
  divide(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int a_lhsComp, const int a_rhsComp);

  /**
    @brief Divide all components in a data holder by a scalar value
    @param[in,out] a_lhs Data to be divided
    @param[in]    a_rhs Divisor. Must have exactly one component. 
  */
  static void
  divideByScalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  /**
    @brief Divide all components in a data holder by a scalar value
    @param[in,out] a_lhs Data to be divided
    @param[in]    a_rhs Divisor. Must have exactly one component. 
  */
  static void
  divideByScalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  /**
    @brief Divide data. If the denominator is zero, set the value to a fallback option.
    @param[in,out] a_numerator   Dividend
    @param[in]    a_denominator Divisor
    @param[in]    a_fallback    Fallback value in case denominator is zero.
    @param[in]    a_vofIter     Pre-built VoF iterators (multiply-cut cells) for each level.
  */
  static void
  divideFallback(EBAMRCellData&                                        a_numerator,
                 const EBAMRCellData&                                  a_denominator,
                 const Real                                            a_fallback,
                 const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Divide data. If the denominator is zero, set the value to a fallback option.
    @param[in,out] a_numerator   Dividend
    @param[in]    a_denominator Divisor
    @param[in]    a_fallback    Fallback value in case denominator is zero.
    @param[in]    a_vofIter     Pre-built VoF iterator (multiply-cut cells).
    @note Does all components.
  */
  static void
  divideFallback(LevelData<EBCellFAB>&       a_numerator,
                 const LevelData<EBCellFAB>& a_denominator,
                 const Real                  a_fallback,
                 LayoutData<VoFIterator>&    a_vofIter);

  /**
    @brief Divide data. If the denominator is zero, set the value to a fallback option.
    @param[in,out] a_numerator   Dividend
    @param[in]    a_denominator Divisor
    @param[in]    a_fallback    Fallback data
    @param[in]    a_vofIter     Pre-built VoF iterators (multiply-cut cells) for each level.
  */
  static void
  divideFallback(EBAMRCellData&                                        a_numerator,
                 const EBAMRCellData&                                  a_denominator,
                 const EBAMRCellData&                                  a_fallback,
                 const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Divide data. If the denominator is zero, set the value to a fallback option.
    @param[in,out] a_numerator   Dividend
    @param[in]    a_denominator Divisor
    @param[in]    a_fallback    Fallback data
    @param[in]    a_vofIter     Pre-built VoF iterator (multiply-cut cells).
    @note Does all components.
  */
  static void
  divideFallback(LevelData<EBCellFAB>&       a_numerator,
                 const LevelData<EBCellFAB>& a_denominator,
                 const LevelData<EBCellFAB>& a_fallback,
                 LayoutData<VoFIterator>&    a_vofIter);

  /**
    @brief Floor values in data holder. This sets all values below a_value to a_value.
    @param[in,out] a_lhs     Data to be floored.
    @param[in]    a_value   Floor value.
    @param[in]    a_vofIter Pre-built VoF iterators (all cut-cells, including ghost region) for each level.
    @note Also floors ghost cells.
  */
  static void
  floor(EBAMRCellData& a_lhs, const Real a_value, const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Floor values in data holder. This sets all values below a_value to a_value.
    @param[in,out] a_lhs     Data to be floored.
    @param[in]    a_value   Floor value.
    @param[in]    a_vofIter Pre-built VoF iterator (all cut-cells, including ghost region).
    @note Also floors ghost cells.
  */
  static void
  floor(LevelData<EBCellFAB>& a_lhs, const Real a_value, LayoutData<VoFIterator>& a_vofIter);

  /**
    @brief Floor values in data holder. This sets all values below a_value to a_value.
    @param[in,out] a_lhs     Data to be floored.
    @param[in]    a_value   Floor value.
    @param[in]    a_vofIter Pre-built VoF iterators (all cut-cells) for each level.
  */
  static void
  floor(EBAMRIVData& a_lhs, const Real a_value, const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Floor values in data holder. This sets all values below a_value to a_value.
    @param[in,out] a_lhs     Data to be floored.
    @param[in]    a_value   Floor value.
    @param[in]    a_vofIter Pre-built VoF iterator (all cut-cells).
  */
  static void
  floor(LevelData<BaseIVFAB<Real>>& a_lhs, const Real a_value, LayoutData<VoFIterator>& a_vofIter);

  /**
    @brief Roof values in data holder. This sets all values above a_value to a_value.
    @param[in,out] a_lhs     Data to be roofed.
    @param[in]    a_value   Roof value.
    @param[in]    a_vofIter Pre-built VoF iterators (all cut-cells, including ghost region) for each level.
    @note Also roofs ghost cells.
  */
  static void
  roof(EBAMRCellData& a_lhs, const Real a_value, const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Roof values in data holder. This sets all values above a_value to a_value.
    @param[in,out] a_lhs     Data to be roofed.
    @param[in]    a_value   Roof value.
    @param[in]    a_vofIter Pre-built VoF iterator (all cut-cells, including ghost region).
    @note Also roofs ghost cells.
  */
  static void
  roof(LevelData<EBCellFAB>& a_lhs, const Real a_value, LayoutData<VoFIterator>& a_vofIter);

  /**
    @brief Roof values in data holder. This sets all values above a_value to a_value
    @param[in,out] a_lhs Data to be roofed.
    @param[in]    a_value Roof value
    @note Also roofs ghost cells. 
  */
  static void
  roof(EBAMRIVData& a_lhs, const Real a_value, const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Roof values in data holder. This sets all values above a_value to a_value.
    @param[in,out] a_lhs     Data to be roofed.
    @param[in]    a_value   Roof value.
    @param[in]    a_vofIter Pre-built VoF iterator (all cut-cells).
  */
  static void
  roof(LevelData<BaseIVFAB<Real>>& a_lhs, const Real a_value, LayoutData<VoFIterator>& a_vofIter);

  /**
    @brief Make a_data = max(a_data1, a_data2).
    @param[out] a_data    Output data. Will be max(a_data1, a_data2) in every cell.
    @param[in]  a_data1   Input data holder.
    @param[in]  a_data2   Input data holder.
    @param[in]  a_vofIter Pre-built VoF iterators (multiply-cut cells) for each level.
  */
  static void
  max(EBAMRCellData&                                        a_data,
      const EBAMRCellData&                                  a_data1,
      const EBAMRCellData&                                  a_data2,
      const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Make a_data = max(a_data1, a_data2).
    @param[out] a_data    Output data. Will be max(a_data1, a_data2) in every cell.
    @param[in]  a_data1   Input data holder.
    @param[in]  a_data2   Input data holder.
    @param[in]  a_vofIter Pre-built VoF iterator (multiply-cut cells).
  */
  static void
  max(LevelData<EBCellFAB>&       a_data,
      const LevelData<EBCellFAB>& a_data1,
      const LevelData<EBCellFAB>& a_data2,
      LayoutData<VoFIterator>&    a_vofIter);

  /**
    @brief Get maximum and minimum value of specified component.
    @param[out] max       Maximum value.
    @param[out] min       Minimum value.
    @param[in]  a_data    Cell-centered data.
    @param[in]  a_comp    Component.
    @param[in]  a_vofIter Pre-built VoF iterators (multiply-cut cells) for each level.
    @note This does the calculation over all levels, including grids that is covered by other grids.
  */
  static void
  getMaxMin(Real&                                                 max,
            Real&                                                 min,
            EBAMRCellData&                                        a_data,
            const int                                             a_comp,
            const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Get maximum and minimum value of specified component.
    @param[out] a_max     Maximum value.
    @param[out] a_min     Minimum value.
    @param[in]  a_data    Cell-centered data.
    @param[in]  a_comp    Component.
    @param[in]  a_vofIter Pre-built VoF iterator (multiply-cut cells).
  */
  static void
  getMaxMin(Real&                    a_max,
            Real&                    a_min,
            LevelData<EBCellFAB>&    a_data,
            const int                a_comp,
            LayoutData<VoFIterator>& a_vofIter);

  /**
    @brief Get maximum and minimum value of specified component
    @param[out] a_max  Maximum value
    @param[out] a_min  Minimum value
    @param[in]  a_data Cell-centered data
    @param[in]  a_comp Component
    @note This does the calculation over all levels, including grids that is covered by other grids. 
  */
  static void
  getMaxMin(Real& a_max, Real& a_min, EBAMRFluxData& a_data, const int a_comp) noexcept;

  /**
    @brief Get maximum and minimum value of specified component using pre-built face iterators.
    @param[out] a_max      Maximum value
    @param[out] a_min      Minimum value
    @param[in]  a_data     Face-centered data
    @param[in]  a_comp     Component
    @param[in]  a_faceIter Pre-built face iterators (cut-cell faces, all directions) for each level.
    @note This does the calculation over all levels, including grids that are covered by other grids.
  */
  static void
  getMaxMin(Real&                                                                   a_max,
            Real&                                                                   a_min,
            EBAMRFluxData&                                                          a_data,
            const int                                                               a_comp,
            Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>& a_faceIter) noexcept;

  /**
    @brief Get maximum and minimum value of specified component
    @param[out] a_max  Maximum value
    @param[out] a_min  Minimum value
    @param[in]  a_data Cell-centered data
    @param[in]  a_comp Component
  */
  static void
  getMaxMin(Real& a_max, Real& a_min, LevelData<EBFluxFAB>& a_data, const int a_comp) noexcept;

  /**
    @brief Get maximum and minimum value of specified component using pre-built face iterators.
    @param[out] a_max      Maximum value
    @param[out] a_min      Minimum value
    @param[in]  a_data     Face-centered data
    @param[in]  a_comp     Component
    @param[in]  a_faceIter Pre-built face iterators (cut-cell faces, all directions).
  */
  static void
  getMaxMin(Real&                                           a_max,
            Real&                                           a_min,
            LevelData<EBFluxFAB>&                           a_data,
            const int                                       a_comp,
            LayoutData<std::array<FaceIterator, SpaceDim>>& a_faceIter) noexcept;

  /**
    @brief Get maximum and minimum value. Assumes single-component data.
    @param[out] a_max  Maximum value
    @param[out] a_min  Minimum value
    @param[in]  a_data Cell-centered data
    @note This does the calculation over all levels, including grids that is covered by other grids. 
  */
  static void
  getMaxMin(Vector<Real>& a_max, Vector<Real>& a_min, Vector<EBAMRCellData>& a_data);

  /**
    @brief Get maximum and minimum value of normed data.
    @param[out] a_max     Maximum value.
    @param[out] a_min     Minimum value.
    @param[in]  data      Data.
    @param[in]  a_vofIter Pre-built VoF iterators (multiply-cut cells) for each level.
    @note This does the calculation over all levels, including grids that is covered by other grids.
  */
  static void
  getMaxMinNorm(Real&                                                 a_max,
                Real&                                                 a_min,
                EBAMRCellData&                                        data,
                const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Get maximum and minimum value of normed data.
    @param[out] a_max     Maximum value.
    @param[out] a_min     Minimum value.
    @param[in]  data      Data.
    @param[in]  a_vofIter Pre-built VoF iterator (multiply-cut cells).
  */
  static void
  getMaxMinNorm(Real& a_max, Real& a_min, LevelData<EBCellFAB>& data, LayoutData<VoFIterator>& a_vofIter);

  /**
    @brief Get maximum and minimum value of normed data.
    @param[out] a_max     Maximum value.
    @param[out] a_min     Minimum value.
    @param[in]  a_data    Data.
    @param[in]  a_vofIter Pre-built VoF iterators (all cut-cells) for each level.
    @note This does the calculation over all levels, including grids that are covered by other grids.
  */
  static void
  getMaxMinNorm(Real&                                                 a_max,
                Real&                                                 a_min,
                EBAMRIVData&                                          a_data,
                const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Get maximum and minimum value of normed data.
    @param[out] a_max     Maximum value.
    @param[out] a_min     Minimum value.
    @param[in]  a_data    Data.
    @param[in]  a_vofIter Pre-built VoF iterator (all cut-cells).
  */
  static void
  getMaxMinNorm(Real& a_max, Real& a_min, LevelData<BaseIVFAB<Real>>& a_data, LayoutData<VoFIterator>& a_vofIter);

  /**
    @brief Invert data
    @note Does all components and does not guard against division by zero.
    @param[in] a_data Data
  */
  static void
  invert(EBAMRFluxData& a_data);

  /**
    @brief Invert data using pre-built face iterators.
    @note Does all components and does not guard against division by zero.
    @param[in,out] a_data     Data
    @param[in]     a_faceIter Pre-built face iterators (cut-cell faces, all directions) for each level.
  */
  static void
  invert(EBAMRFluxData&                                                                  a_data,
         Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>& a_faceIter);

  /**
    @brief Invert data
    @note Does all components and does not guard against division by zero.
    @param[in] a_data Data
  */
  static void
  invert(LevelData<EBFluxFAB>& a_data);

  /**
    @brief Invert data using pre-built face iterators.
    @note Does all components and does not guard against division by zero.
    @param[in,out] a_data     Data
    @param[in]     a_faceIter Pre-built face iterators (cut-cell faces, all directions).
  */
  static void
  invert(LevelData<EBFluxFAB>& a_data, LayoutData<std::array<FaceIterator, SpaceDim>>& a_faceIter);

  /**
    @brief Invert data
    @note Does all components and does not guard against division by zero.
    @param[in] a_data Data
  */
  static void
  invert(EBFluxFAB& a_data);

  /**
    @brief Compute the kappa-weighted sum of the input data. Assumes one component.
    @param[out] a_mass    Total mass, computed as kappa * lhs.
    @param[in]  a_lhs     Cell data.
    @param[in]  a_comp    Component.
    @param[in]  a_vofIter Pre-built VoF iterator (all cut-cells).
  */
  static void
  kappaSum(Real& a_mass, const LevelData<EBCellFAB>& a_lhs, const int a_comp, LayoutData<VoFIterator>& a_vofIter);

  /**
    @brief Scale data by volume fraction.
    @param[in,out] a_data    Data to be scaled. On output we have a_data = kappa * a_data.
    @param[in]     a_vofIter Pre-built VoF iterators (all cut-cells) for each level.
  */
  static void
  kappaScale(EBAMRCellData& a_data, const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter) noexcept;

  /**
    @brief Scale data by volume fraction.
    @param[in,out] a_data    Data to be scaled. On output we have a_data = kappa * a_data.
    @param[in]     a_vofIter Pre-built VoF iterator (all cut-cells).
  */
  static void
  kappaScale(LevelData<EBCellFAB>& a_data, LayoutData<VoFIterator>& a_vofIter) noexcept;

  /**
    @brief Scale data by volume fraction
    @param[in,out] a_data Data to be scaled. On output we have a_data = kappa * a_data
  */
  static void
  kappaScale(MFAMRCellData& a_data) noexcept;

  /**
    @brief Scale data by volume fraction
    @param[in,out] a_data Data to be scaled. On output we have a_data = kappa * a_data
  */
  static void
  kappaScale(LevelData<MFCellFAB>& a_data) noexcept;

  /**
    @brief Scale data by dx^SpaceDim
    @param[in,out] a_data Data to be scaled. On output we have a_data = dx^SpaceDim * a_data
    @param[in] a_dx Description
  */
  static void
  volumeScale(EBAMRCellData& a_data, const Vector<Real>& a_dx);

  /**
    @brief Multiply data holder by another data holder
    @param[in,out] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  /**
    @brief Multiply data holder by another data holder
    @param[in,out] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  /**
    @brief Multiply data holder by another data holder
    @param[in,out] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(EBAMRFluxData& a_lhs, const EBAMRFluxData& a_rhs);

  /**
    @brief Multiply data holder by another data holder
    @param[in,out] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(LevelData<EBFluxFAB>& a_lhs, const LevelData<EBFluxFAB>& a_rhs);

  /**
    @brief Multiply data holder by another data holder
    @param[in,out] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs);

  /**
    @brief Multiply data holder by another data holder
    @param[in,out] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs);

  /**
    @brief Multiply data holder by another data holder
    @param[in,out] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note This multiplies all components in a_lhs with the single component in a_rhs
  */
  static void
  multiplyScalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  /**
    @brief Multiply data holder by another data holder
    @param[in,out] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note This multiplies all components in a_lhs with the single component in a_rhs
  */
  static void
  multiplyScalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  /**
    @brief Multiply data holder by another data holder
    @param[in,out] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note This multiplies all components in a_lhs with the single component in a_rhs
  */
  static void
  multiplyScalar(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs);

  /**
    @brief Multiply data holder by another data holder
    @param[in,out] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note This multiplies all components in a_lhs with the single component in a_rhs
  */
  static void
  multiplyScalar(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs);

  /**
    @brief Compute specified norm of data.
    @param[in]  a_data    Data to compute the norm of.
    @param[in]  a_p       Norm factor, e.g. a_p = 2 => squared norm. 0 => inf norm.
    @param[in]  a_comp    Component.
    @param[in]  a_vofIter Pre-built VoF iterator (all cut-cells).
    @note Does not do kappa-weighting.
    @return Return value.
  */
  static Real
  norm(const LevelData<EBCellFAB>& a_data, const int a_p, const int a_comp, LayoutData<VoFIterator>& a_vofIter);

  /**
    @brief General addition operator for adding together data. The user can choose which components to add. 
    @param[in,out] a_lhs     Output data component
    @param[in,out] a_rhs     Data to be added to a_lhs
    @param[in]    a_srcComp Component in a_rhs to add from.
    @param[in]    a_dstComp Component in a_lhs to add into.
    @param[in]    a_numComp Number of components to add. 
    @note This will add a_lhs(a_dstComp + icomp) += a_rhs(a_srcComp + icomp) for icomp = 0,1,2,..numComp-1.
  */
  static void
  plus(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_srcComp, const int a_dstComp, const int a_numComp);

  /**
    @brief General addition operator for adding together data. The user can choose which components to add. 
    @param[in,out] a_lhs     Output data component
    @param[in,out] a_rhs     Data to be added to a_lhs
    @param[in]    a_srcComp Component in a_rhs to add from.
    @param[in]    a_dstComp Component in a_lhs to add into.
    @param[in]    a_numComp Number of components to add. 
    @note This will add a_lhs(a_dstComp + icomp) += a_rhs(a_srcComp + icomp) for icomp = 0,1,2,..numComp-1.
  */
  static void
  plus(LevelData<EBCellFAB>&       a_lhs,
       const LevelData<EBCellFAB>& a_rhs,
       const int                   a_srcComp,
       const int                   a_dstComp,
       const int                   a_numComp);

  /**
    @brief Scale data by factor
    @param[in,out] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(MFAMRCellData& a_lhs, const Real& a_scale) noexcept;

  /**
    @brief Scale data by factor
    @param[in,out] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale) noexcept;

  /**
    @brief Scale data by factor
    @param[in,out] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(MFAMRFluxData& a_lhs, const Real& a_scale);

  /**
    @brief Scale data by factor
    @param[in,out] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(LevelData<MFFluxFAB>& a_lhs, const Real& a_scale);

  /**
    @brief Scale data by factor
    @param[in,out] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(EBAMRIVData& a_lhs, const Real& a_scale);

  /**
    @brief Scale data by factor
    @param[in,out] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(EBAMRCellData& a_lhs, const Real a_scale) noexcept;

  /**
    @brief Scale data by factor
    @param[in,out] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(LevelData<EBCellFAB>& a_lhs, const Real a_scale) noexcept;

  /**
    @brief Scale data by factor
    @param[in,out] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(EBAMRFluxData& a_lhs, const Real a_scale);

  /**
    @brief Scale data by factor
    @param[in,out] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(LevelData<EBFluxFAB>& a_lhs, const Real a_scale);

  /**
    @brief Scale data by factor
    @param[in,out] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(LevelData<BaseIVFAB<Real>>& a_lhs, const Real& a_scale);

  /**
    @brief Set value in covered cells. Does specified component
    @param[in,out] a_lhs   Input/output data. 
    @param[in]   a_comp  Component to set
    @param[in]   a_value Value to set
  */
  static void
  setCoveredValue(EBAMRCellData& a_lhs, const int a_comp, const Real a_value);

  /**
    @brief Set value in covered cells. Does specified component
    @param[in,out] a_lhs   Input/output data. 
    @param[in]   a_comp  Component to set
    @param[in]   a_value Value to set
  */
  static void
  setCoveredValue(LevelData<EBCellFAB>& a_lhs, const int a_comp, const Real a_value);

  /**
    @brief Version of setCoveredValue which does all components
    @param[in,out] a_lhs   Data 
    @param[in]    a_value Value in covered grid cells. 
  */
  static void
  setCoveredValue(EBAMRCellData& a_lhs, const Real a_value);

  /**
    @brief Version of setCoveredValue which does all component
    @param[in,out] a_lhs   Data 
    @param[in]    a_value Value in covered grid cells. 
  */
  static void
  setCoveredValue(LevelData<EBCellFAB>& a_lhs, const Real a_value);

  /**
    @brief Set value on covered faces.
    @param[in,out] a_lhs   Input/output data. 
    @param[in]   a_comp  Component to set
    @param[in]   a_value Value to set
  */
  static void
  setCoveredValue(EBAMRFluxData& a_lhs, const int a_comp, const Real a_value) noexcept;

  /**
    @brief Set value on covered faces.
    @param[in,out] a_lhs   Input/output data. 
    @param[in]   a_comp  Component to set
    @param[in]   a_value Value to set
  */
  static void
  setCoveredValue(LevelData<EBFluxFAB>& a_lhs, const int a_comp, const Real a_value) noexcept;

  /**
    @brief Set all covered face values.
    @param[in,out] a_lhs   Data 
    @param[in]    a_value Value in covered grid cells. 
  */
  static void
  setCoveredValue(EBAMRFluxData& a_lhs, const Real a_value) noexcept;

  /**
    @brief Set all covered face values.
    @param[in,out] a_lhs   Data 
    @param[in]    a_value Value in covered grid cells. 
  */
  static void
  setCoveredValue(LevelData<EBFluxFAB>& a_lhs, const Real a_value) noexcept;

  /**
    @brief Set value in cells that are covered by a finer grid.
    @param[in,out] a_lhs    AMR data
    @param[in]    a_refRat Refinement ratios between levels.
    @param[in]    a_value  Value to set in cells covered by a finer grid. 
  */
  static void
  setInvalidValue(EBAMRCellData& a_lhs, const Vector<int>& a_refRat, const Real a_value);

  /**
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set
    @param[in]  a_function Function to use for setting the value. 
    @param[in]  a_probLo   Lower-left corner of physical domain. 
    @param[in]  a_dx       Grid resolutions
    @param[in]  a_comp     Component to set
    @note Uses a VofIterator everywhere to this might be slow!
  */
  static void
  setValue(MFAMRCellData&                             a_lhs,
           const std::function<Real(const RealVect)>& a_function,
           const RealVect&                            a_probLo,
           const Vector<Real>&                        a_dx,
           const int                                  a_comp);

  /**
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set
    @param[in]  a_function Function to use for setting the value. 
    @param[in]  a_probLo   Lower-left corner of physical domain. 
    @param[in]  a_dx       Grid resolutions
    @param[in]  a_comp     Component to set
    @note Uses a VofIterator everywhere to this might be slow!
  */
  static void
  setValue(LevelData<MFCellFAB>&                      a_lhs,
           const std::function<Real(const RealVect)>& a_function,
           const RealVect                             a_probLo,
           const Real                                 a_dx,
           const int                                  a_comp);

  /**
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set.
    @param[in]  a_function Function to use for setting the value.
    @param[in]  a_probLo   Lower-left corner of physical domain.
    @param[in]  a_dx       Grid resolutions.
    @param[in]  a_comp     Component to set.
    @param[in]  a_vofIter  Pre-built VoF iterators (multiply-cut cells) for each level.
  */
  static void
  setValue(EBAMRCellData&                                        a_lhs,
           const std::function<Real(const RealVect)>&            a_function,
           const RealVect&                                       a_probLo,
           const Vector<Real>&                                   a_dx,
           const int                                             a_comp,
           const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set.
    @param[in]  a_function Function to use for setting the value.
    @param[in]  a_probLo   Lower-left corner of physical domain.
    @param[in]  a_dx       Grid resolutions.
    @param[in]  a_comp     Component to set.
    @param[in]  a_vofIter  Pre-built VoF iterator (multiply-cut cells).
  */
  static void
  setValue(LevelData<EBCellFAB>&                      a_lhs,
           const std::function<Real(const RealVect)>& a_function,
           const RealVect                             a_probLo,
           const Real                                 a_dx,
           const int                                  a_comp,
           LayoutData<VoFIterator>&                   a_vofIter);

  /**
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set
    @param[in]  a_function Function to use for setting the value. 
    @param[in]  a_probLo   Lower-left corner of physical domain. 
    @param[in]  a_dx       Grid resolutions
    @param[in]  a_comp     Component to set
    @note Uses a VofIterator everywhere to this might be slow!
  */
  static void
  setValue(EBAMRFluxData&                             a_lhs,
           const std::function<Real(const RealVect)>& a_function,
           const RealVect&                            a_probLo,
           const Vector<Real>&                        a_dx,
           const int                                  a_comp);

  /**
    @brief Polymorphic set value function using pre-built face iterators.
    @param[out] a_lhs      Data to set.
    @param[in]  a_function Function to use for setting the value.
    @param[in]  a_probLo   Lower-left corner of physical domain.
    @param[in]  a_dx       Grid resolutions.
    @param[in]  a_comp     Component to set.
    @param[in]  a_faceIter Pre-built face iterators (cut-cell faces, all directions) for each level.
  */
  static void
  setValue(EBAMRFluxData&                                                                  a_lhs,
           const std::function<Real(const RealVect)>&                                     a_function,
           const RealVect&                                                                 a_probLo,
           const Vector<Real>&                                                             a_dx,
           const int                                                                       a_comp,
           Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>& a_faceIter);

  /**
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set
    @param[in]  a_function Function to use for setting the value.
    @param[in]  a_probLo   Lower-left corner of physical domain.
    @param[in]  a_dx       Grid resolutions
    @param[in]  a_comp     Component to set
    @note Uses a VofIterator everywhere to this might be slow!
  */
  static void
  setValue(LevelData<EBFluxFAB>&                      a_lhs,
           const std::function<Real(const RealVect)>& a_function,
           const RealVect                             a_probLo,
           const Real                                 a_dx,
           const int                                  a_comp);

  /**
    @brief Polymorphic set value function using pre-built face iterators.
    @param[out] a_lhs      Data to set.
    @param[in]  a_function Function to use for setting the value.
    @param[in]  a_probLo   Lower-left corner of physical domain.
    @param[in]  a_dx       Grid resolution.
    @param[in]  a_comp     Component to set.
    @param[in]  a_faceIter Pre-built face iterators (cut-cell faces, all directions).
  */
  static void
  setValue(LevelData<EBFluxFAB>&                                   a_lhs,
           const std::function<Real(const RealVect)>&              a_function,
           const RealVect                                          a_probLo,
           const Real                                              a_dx,
           const int                                               a_comp,
           LayoutData<std::array<FaceIterator, SpaceDim>>&         a_faceIter);

  /**
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set.
    @param[in]  a_function Function to use for setting the value.
    @param[in]  a_probLo   Lower-left corner of physical domain.
    @param[in]  a_dx       Grid resolutions.
    @param[in]  a_comp     Component to set.
    @param[in]  a_vofIter  Pre-built VoF iterators (all cut-cells) for each level.
  */
  static void
  setValue(EBAMRIVData&                                          a_lhs,
           const std::function<Real(const RealVect)>&            a_function,
           const RealVect&                                       a_probLo,
           const Vector<Real>&                                   a_dx,
           const int                                             a_comp,
           const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set.
    @param[in]  a_function Function to use for setting the value.
    @param[in]  a_probLo   Lower-left corner of physical domain.
    @param[in]  a_dx       Grid resolution.
    @param[in]  a_comp     Component to set.
    @param[in]  a_vofIter  Pre-built VoF iterator (all cut-cells).
  */
  static void
  setValue(LevelData<BaseIVFAB<Real>>&                a_lhs,
           const std::function<Real(const RealVect)>& a_function,
           const RealVect                             a_probLo,
           const Real                                 a_dx,
           const int                                  a_comp,
           LayoutData<VoFIterator>&                   a_vofIter);

  /**
    @brief Polymorphic set value function. Assumes that a_lhs has SpaceDim components and sets all those components from the input function.
    @param[out] a_lhs      Data to set.
    @param[in]  a_function Function to use for setting the value.
    @param[in]  a_probLo   Lower-left corner of physical domain.
    @param[in]  a_dx       Grid resolutions.
    @param[in]  a_vofIter  Pre-built VoF iterators (multiply-cut cells) for each level.
  */
  static void
  setValue(EBAMRCellData&                                        a_lhs,
           const std::function<RealVect(const RealVect)>&        a_function,
           const RealVect&                                       a_probLo,
           const Vector<Real>&                                   a_dx,
           const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Polymorphic set value function. Assumes that a_lhs has SpaceDim components and sets all those components from the input function.
    @param[out] a_lhs      Data to set.
    @param[in]  a_function Function to use for setting the value.
    @param[in]  a_probLo   Lower-left corner of physical domain.
    @param[in]  a_dx       Grid resolutions.
    @param[in]  a_vofIter  Pre-built VoF iterator (multiply-cut cells).
  */
  static void
  setValue(LevelData<EBCellFAB>&                          a_lhs,
           const std::function<RealVect(const RealVect)>& a_function,
           const RealVect                                 a_probLo,
           const Real                                     a_dx,
           LayoutData<VoFIterator>&                       a_vofIter);

  /**
    @brief Set value in data
    @param[out] a_data  Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(EBAMRCellData& a_data, const Real& a_value);

  /**
    @brief Set value in data with specified component. 
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
    @param[in] a_comp Description
  */
  static void
  setValue(EBAMRCellData& a_lhs, const Real a_value, const int a_comp);

  /**
    @brief Set value in data with specified component. 
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
    @param[in] a_comp Description
  */
  static void
  setValue(LevelData<EBCellFAB>& a_lhs, const Real a_value, const int a_comp);

  /**
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<EBCellFAB>& a_lhs, const Real a_value);

  /**
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<EBFluxFAB>& a_lhs, const Real a_value);

  /**
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<BaseIVFAB<Real>>& a_lhs, const Real a_value);

  /**
    @brief Set value in data
    @param[out] a_data  Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(EBAMRFluxData& a_data, const Real& a_value);

  /**
    @brief Set value in data
    @param[out] a_data  Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(EBAMRIVData& a_data, const Real& a_value);

  /**
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(MFAMRCellData& a_lhs, const Real& a_value) noexcept;

  /**
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<MFCellFAB>& a_lhs, const Real& a_value) noexcept;

  /**
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(MFAMRFluxData& a_lhs, const Real& a_value);

  /**
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<MFFluxFAB>& a_lhs, const Real& a_value);

  /**
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(MFAMRIVData& a_lhs, const Real& a_value);

  /**
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<MFBaseIVFAB>& a_lhs, const Real& a_value);

  /**
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(EBAMRIFData& a_lhs, const Real a_value);

  /**
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<DomainFluxIFFAB>& a_lhs, const Real a_value);

  /**
    @brief Compute the sum of the input value.
    @param[in,out] a_value On output, contains the sum of a_value over all MPI processes.
  */
  static void
  sum(Real& a_value);

  /**
    @brief Compute the square root of the input data
    @param[in,out] a_lhs On output, contains sqrt(a_lhs)
    @note Does not do ghosted faces. 
  */
  static void
  squareRoot(EBAMRFluxData& a_lhs);

  /**
    @brief Compute the square root of the input data using pre-built face iterators.
    @param[in,out] a_lhs      On output, contains sqrt(a_lhs).
    @param[in]     a_faceIter Pre-built face iterators (cut-cell faces, all directions) for each level.
    @note Does not do ghosted faces.
  */
  static void
  squareRoot(EBAMRFluxData&                                                                  a_lhs,
             Vector<RefCountedPtr<LayoutData<std::array<FaceIterator, SpaceDim>>>>& a_faceIter);

  /**
    @brief Compute the square root of the input data
    @param[in,out] a_lhs On output, contains sqrt(a_lhs)
    @note Does not do ghosted faces.
  */
  static void
  squareRoot(LevelData<EBFluxFAB>& a_lhs);

  /**
    @brief Compute the square root of the input data using pre-built face iterators.
    @param[in,out] a_lhs      On output, contains sqrt(a_lhs).
    @param[in]     a_faceIter Pre-built face iterators (cut-cell faces, all directions).
    @note Does not do ghosted faces.
  */
  static void
  squareRoot(LevelData<EBFluxFAB>& a_lhs, LayoutData<std::array<FaceIterator, SpaceDim>>& a_faceIter);

  /**
    @brief Take the square root of the input data
    @param[in,out] a_lhs Input data
  */
  static void
  squareRoot(MFAMRCellData& a_lhs);

  /**
    @brief Take the square root of the input data
    @param[in,out] a_lhs Input data
  */
  static void
  squareRoot(LevelData<MFCellFAB>& a_lhs);

  /**
    @brief Compute the vector length of a data holder. Sets a_lhs = |a_rhs| where a_rhs contains SpaceDim components.
    @param[out] a_lhs     Vector length.
    @param[in]  a_rhs     Input vector.
    @param[in]  a_vofIter Pre-built VoF iterators (multiply-cut cells) for each level.
    @note a_lhs must have one component and a_rhs must have SpaceDim components.
  */
  static void
  vectorLength(EBAMRCellData&                                        a_lhs,
               const EBAMRCellData&                                  a_rhs,
               const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Compute the vector length of a data holder. Sets a_lhs = |a_rhs| where a_rhs contains SpaceDim components.
    @param[out] a_lhs     Vector length.
    @param[in]  a_rhs     Input vector.
    @param[in]  a_vofIter Pre-built VoF iterator (multiply-cut cells).
    @note a_lhs must have one component and a_rhs must have SpaceDim components.
  */
  static void
  vectorLength(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, LayoutData<VoFIterator>& a_vofIter);

  /**
    @brief Compute the vector length of a data holder. Sets a_lhs = |a_rhs| where a_rhs contains SpaceDim components.
    @param[out] a_lhs     Vector length.
    @param[in]  a_rhs     Input vector.
    @param[in]  box       Cell box.
    @param[in]  a_vofIter Pre-built VoF iterator (multiply-cut cells) for this patch.
    @note a_lhs must have one component and a_rhs must have SpaceDim components.
  */
  static void
  vectorLength(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& box, VoFIterator& a_vofIter);

  /**
    @brief Compute the squared vector length of a data holder. Sets a_lhs = |a_rhs|^2 where a_rhs contains SpaceDim components.
    @param[out] a_lhs     Vector length.
    @param[in]  a_rhs     Input vector.
    @param[in]  a_vofIter Pre-built VoF iterators (multiply-cut cells) for each level.
    @note a_lhs must have one component and a_rhs must have SpaceDim components.
  */
  static void
  vectorLength2(EBAMRCellData&                                        a_lhs,
                const EBAMRCellData&                                  a_rhs,
                const Vector<RefCountedPtr<LayoutData<VoFIterator>>>& a_vofIter);

  /**
    @brief Compute the squared vector length of a data holder. Sets a_lhs = |a_rhs|^2 where a_rhs contains SpaceDim components.
    @param[out] a_lhs     Vector length.
    @param[in]  a_rhs     Input vector.
    @param[in]  a_vofIter Pre-built VoF iterator (multiply-cut cells).
    @note a_lhs must have one component and a_rhs must have SpaceDim components.
  */
  static void
  vectorLength2(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, LayoutData<VoFIterator>& a_vofIter);

  /**
    @brief Compute the squared vector length of a data holder. Sets a_lhs = |a_rhs|^2 where a_rhs contains SpaceDim components.
    @param[out] a_lhs     Vector length.
    @param[in]  a_rhs     Input vector.
    @param[in]  box       Cell box.
    @param[in]  a_vofIter Pre-built VoF iterator (multiply-cut cells) for this patch.
    @note a_lhs must have one component and a_rhs must have SpaceDim components.
  */
  static void
  vectorLength2(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& box, VoFIterator& a_vofIter);

  /**
    @brief Compute the tightest possible valid box around a cut-cell volume
    @param[out] a_lo       Low corner
    @param[out] a_hi       High corner
    @param[in]  a_normal   EB face normal
    @param[in]  a_centroid EB centroid
  */
  static void
  computeMinValidBox(RealVect& a_lo, RealVect& a_hi, const RealVect& a_normal, const RealVect& a_centroid);

  /**
    @brief Check if all corners/points are inside an EB
    @param[in]  a_normal   EB face normal
    @param[in]  a_centroid EB centroid
    @return Returns true if all points are inside the EB and false otherwise. 
    @param[in] a_corners Corners
  */
  static bool
  allCornersInsideEb(const Vector<RealVect>& a_corners, const RealVect& a_normal, const RealVect& a_centroid);

  /**
    @brief Shift points by a specified distance
    @details On output, every point in a_corners is shifted by a_distance. 
    @param[in,out] a_corners  Points to be shifted
    @param[in]    a_distance Shift vector
  */
  static void
  shiftCorners(Vector<RealVect>& a_corners, const RealVect& a_distance);
};

#include <CD_NamespaceFooter.H>

#include <CD_DataOpsImplem.H>

#endif
