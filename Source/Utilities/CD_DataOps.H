/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DataOps.H
  @brief  Agglomeration of useful data operations
  @author Robert Marskar
*/

#ifndef CD_DataOps_H
#define CD_DataOps_H

// Our includes
#include <CD_Average.H>
#include <CD_EBAMRData.H>
#include <CD_Decorations.H>
#include <CD_MFInterfaceFAB.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Agglomeration of useful data operations.
  @note All methods in this class assumes that data are allocated over the same realm. If you are trying to do something across realms, you need to copy
  the data onto realms that lie on the same data holders first!.
*/
class DataOps
{
public:
  /*!
    @brief Set value in an MFInterfaceFAB data holder
    @param[in] a_lhs   Data holder
    @param[in] a_value Value to set
  */
  template <typename T>
  static void
  setValue(LevelData<MFInterfaceFAB<T>>& a_lhs, const T& a_value);

  /*!
    @brief Sign function. Returns +/- if the value is > 0 or < 0
    @param[in] a_value Value to evaluate. 
  */
  template <typename T>
  static int
  sgn(const T a_value);

  /*!
    @brief Routine which computes the average of a cell-centered quantity on faces for the normal component only.
    @param[out] a_faceData Face data. Must have one component. 
    @param[in]  a_cellData Cell data. Must have SpaceDim components. 
    @param[in]  a_domains  AMR domains.
    @note This routine also fills one tangential ghost face.
  */
  static void
  averageCellVelocityToFaceVelocity(EBAMRFluxData&               a_faceData,
                                    const EBAMRCellData&         a_cellData,
                                    const Vector<ProblemDomain>& a_domains);

  /*!
    @brief Routine which computes the average of a cell-centered quantity on faces for the normal component only.
    @param[out] a_faceData Face data. Must have one component. 
    @param[in]  a_cellData Cell data. Must have SpaceDim components. 
    @param[in]  a_domain   Problem domain
    @param[in]  a_tanGhost Number of ghost faces to fill
    @param[in]  a_average  Averaging method
    @note This routine also fills one tangential ghost face.
  */
  static void
  averageCellVelocityToFaceVelocity(LevelData<EBFluxFAB>&       a_faceData,
                                    const LevelData<EBCellFAB>& a_cellData,
                                    const ProblemDomain&        a_domain);

  /*!
    @brief Average all components of the cell-centered data to faces.
    @details Calls the other version; uses arihtmetic averaging and does not fill ghost faces. 
    @param[out] a_faceData Face data. 
    @param[in]  a_cellData Cell data. 
    @param[in]  a_domains  AMR domains.
    @note This does not fill "tangential" ghost faces.
  */
  static void
  averageCellToFace(EBAMRFluxData& a_faceData, const EBAMRCellData& a_cellData, const Vector<ProblemDomain>& a_domains);

  /*!
    @brief Routine which averages a cell-centered component to faces. 
    @details User can choose what type of averaging to use, as well as choosing the number of ghost face to fill. These
    faces are "tangential ghost faces", i.e. if the face direction is along +x, we extend the box with faces in the y- and z-
    directions.
    @param[out] a_faceData Face data. 
    @param[in]  a_cellData Cell data. 
    @param[in]  a_domain   Problem domain
    @param[in]  a_tanGhost Number of ghost faces to fill
    @param[in]  a_faceVars Face variables
    @param[in]  a_cellVars Cell variables
    @param[in]  a_average  Averaging method
  */
  static void
  averageCellToFace(EBAMRFluxData&               a_faceData,
                    const EBAMRCellData&         a_cellData,
                    const Vector<ProblemDomain>& a_domains,
                    const int                    a_tanGhosts,
                    const Interval&              a_faceVars,
                    const Interval&              a_cellVars,
                    const Average&               a_average);

  /*!
    @brief Routine which averages a cell-centered component to faces. 
    @details User can choose what type of averaging to use, as well as choosing the number of ghost face to fill. These
    faces are "tangential ghost faces", i.e. if the face direction is along +x, we extend the box with faces in the y- and z-
    directions.
    @param[out] a_faceData Face data. 
    @param[in]  a_cellData Cell data. 
    @param[in]  a_domain   Problem domain
    @param[in]  a_tanGhost Number of ghost faces to fill
    @param[in]  a_faceVars Face variables
    @param[in]  a_cellVars Cell variables
    @param[in]  a_average  Averaging method
  */
  static void
  averageCellToFace(LevelData<EBFluxFAB>&       a_faceData,
                    const LevelData<EBCellFAB>& a_cellData,
                    const ProblemDomain&        a_domain,
                    const int                   a_tanGhosts,
                    const Interval&             a_faceVars,
                    const Interval&             a_cellVars,
                    const Average&              a_average);

  /*!
    @brief Average all components of face centered data to cell centers.
    @param[out] a_cellData Cell data. 
    @param[int] a_faceData Face data. 
    @param[in]  a_domains  AMR domains.
    @note Only valid cells are filled (not ghost cells)
  */
  static void
  averageFaceToCell(EBAMRCellData& a_cellData, const EBAMRFluxData& a_fluxData, const Vector<ProblemDomain>& a_domains);

  /*!
    @brief Average all components of face centered data to cell centers.
    @param[out] a_cellData Cell data. 
    @param[int] a_faceData Face data. 
    @param[in]  a_domains  Domain
    @note Only valid cells are filled (not ghost cells)
  */
  static void
  averageFaceToCell(LevelData<EBCellFAB>&       a_cellData,
                    const LevelData<EBFluxFAB>& a_fluxData,
                    const ProblemDomain&        a_domain);

  /*!
    @brief Compote the cell-wise dot product between two data holders. 
    @param[out] a_result Result. Holds the dot product in each cell.
    @param[in]  a_data1 First data holder.
    @param[in]  a_data2 Other data holder.
    @note The input data data1 and data2 must have the same number of input components. Summation is done over all components (does not have to be SpaceDim).
  */
  static void
  dotProduct(MFAMRCellData& a_result, const MFAMRCellData& a_data1, const MFAMRCellData& a_data2);

  /*!
    @brief Compote the cell-wise dot product between two data holders. 
    @details This calls the single-phase dotProduct patch versions with EBCellFABs on each phase.
    @param[out] a_result Result. Holds the dot product in each cell.
    @param[in]  a_data1 First data holder.
    @param[in]  a_data2 Other data holder.
    @note The input data data1 and data2 must have the same number of input components. Summation is done over all components (does not have to be SpaceDim).
  */
  static void
  dotProduct(LevelData<MFCellFAB>& a_result, const LevelData<MFCellFAB>& a_data1, const LevelData<MFCellFAB>& a_data2);

  /*!
    @brief Compote the cell-wise dot product between two data holders. 
    @param[out] a_result Result. Holds the dot product in each cell.
    @param[in]  a_data1 First data holder.
    @param[in]  a_data2 Other data holder.
    @note The input data data1 and data2 must have the same number of input components. Summation is done over all components (does not have to be SpaceDim).
  */
  static void
  dotProduct(EBAMRCellData& a_result, const EBAMRCellData& a_data1, const EBAMRCellData& a_data2);

  /*!
    @brief Compote the cell-wise dot product between two data holders. 
    @param[out] a_result Result. Holds the dot product in each cell.
    @param[in]  a_data1 First data holder.
    @param[in]  a_data2 Other data holder.
    @note The input data data1 and data2 must have the same number of input components. Summation is done over all components (does not have to be SpaceDim).
  */
  static void
  dotProduct(LevelData<EBCellFAB>& a_result, const LevelData<EBCellFAB>& a_data1, const LevelData<EBCellFAB>& a_data2);

  /*!
    @brief Compote the cell-wise dot product between two data holders. 
    @param[out] a_result Result. Holds the dot product in each cell.
    @param[in]  a_data1 First data holder.
    @param[in]  a_data2 Other data holder.
    @note The input data data1 and data2 must have the same number of input components. Summation is done over all components (does not have to be SpaceDim).
  */
  static void
  dotProduct(EBCellFAB& a_result, const EBCellFAB& a_data1, const EBCellFAB& a_data2, const Box& a_box);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(MFAMRCellData& a_lhs, const MFAMRCellData& a_rhs, const Real a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<MFCellFAB>& a_lhs, const LevelData<MFCellFAB>& a_rhs, const Real a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const Real& a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real& a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRFluxData& a_lhs, const EBAMRFluxData& a_rhs, const Real& a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<EBFluxFAB>& a_lhs, const LevelData<EBFluxFAB>& a_rhs, const Real& a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBFluxFAB& a_lhs, const EBFluxFAB& a_rhs, const Real& a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs, const Real& a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs, const Real& a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRIFData& a_lhs, const EBAMRIFData& a_rhs, const Real& a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<DomainFluxIFFAB>& a_lhs, const LevelData<DomainFluxIFFAB>& a_rhs, const Real& a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @details This versions only adds irregular data from an EBAMRIVFAB. It ignores regular cells. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRCellData& a_lhs, const EBAMRIVData& a_rhs, const Real a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @details This versions only adds irregular data from an EBAMRIVFAB. It ignores regular cells. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<EBCellFAB>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs, const Real a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @details This versions only adds irregular data form an EBCellFAB to an BaseIVFAB. It ignores regular cells. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(EBAMRIVData& a_lhs, const EBAMRCellData& a_rhs, const Real a_scale);

  /*!
    @brief Function which increments data in the form a_lhs = a_lhs + a_rhs*a_scale for all components. 
    @details This versions only adds irregular data form an EBCellFAB to an BaseIVFAB. It ignores regular cells. 
    @param[inout] a_lhs   Data to be incremented.
    @param[in]    a_rhs   Incrementation.
    @param[in]    a_scale Scale factor.
  */
  static void
  incr(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const Real a_scale);

  /*!
    @brief Copy data from one data holder to another.
    @param[out] a_dst Destination data holder
    @param[in ] a_src Source data holder
  */
  static void
  copy(MFAMRCellData& a_dst, const MFAMRCellData& a_src);

  /*!
    @brief Copy data from one data holder to another.
    @param[out] a_dst Destination data holder
    @param[in ] a_src Source data holder
  */
  static void
  copy(EBAMRCellData& a_dst, const EBAMRCellData& a_src);

  /*!
    @brief Copy data from one data holder to another.
    @param[out] a_dst Destination data holder
    @param[in ] a_src Source data holder
  */
  static void
  copy(EBAMRIVData& a_dst, const EBAMRIVData& a_src);

  /*!
    @brief Divide one data holder by another. Does specified components.
    @param[inout] a_lhs Data to be divided
    @param[in]    a_rhs Divisor.
    @param[in]    a_lhsComp Destination data component
    @param[in]    a_rhsComp Source data component
  */
  static void
  divide(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_lhsComp, const int a_rhsComp);

  /*!
    @brief Divide one data holder by another. Does specified components.
    @param[inout] a_lhs Data to be divided
    @param[in]    a_rhs Divisor.
    @param[in]    a_lhsComp Destination data component
    @param[in]    a_rhsComp Source data component
  */
  static void
  divide(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs, const int a_lhsComp, const int a_rhsComp);

  /*!
    @brief Divide all components in a data holder by a scalar value
    @param[inout] a_lhs Data to be divided
    @param[in]    a_rhs Divisor. Must have exactly one component. 
  */
  static void
  divideByScalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  /*!
    @brief Divide all components in a data holder by a scalar value
    @param[inout] a_lhs Data to be divided
    @param[in]    a_rhs Divisor. Must have exactly one component. 
  */
  static void
  divideByScalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  /*!
    @brief Divide data. If the denominator is zero, set the value to a fallback option
    @param[inout] a_numerator   Dividend
    @param[in]    a_denominator Divisor
    @param[in]    a_fallback    Fallback value in case denominator is zero.
  */
  static void
  divideFallback(EBAMRCellData& a_numerator, const EBAMRCellData& a_denominator, const Real a_fallback);

  /*!
    @brief Divide data. If the denominator is zero, set the value to a fallback option
    @param[inout] a_numerator   Dividend
    @param[in]    a_denominator Divisor
    @param[in]    a_fallback    Fallback value in case denominator is zero.
    @note Does all components. 
  */
  static void
  divideFallback(LevelData<EBCellFAB>& a_numerator, const LevelData<EBCellFAB>& a_denominator, const Real a_fallback);

  /*!
    @brief Divide data. If the denominator is zero, set the value to a fallback option
    @param[inout] a_numerator   Dividend
    @param[in]    a_denominator Divisor
    @param[in]    a_fallback    Fallback data
  */
  static void
  divideFallback(EBAMRCellData& a_numerator, const EBAMRCellData& a_denominator, const EBAMRCellData& a_fallback);

  /*!
    @brief Divide data. If the denominator is zero, set the value to a fallback option
    @param[inout] a_numerator   Dividend
    @param[in]    a_denominator Divisor
    @param[in]    a_fallback    Fallback data
    @note Does all components. 
  */
  static void
  divideFallback(LevelData<EBCellFAB>&       a_numerator,
                 const LevelData<EBCellFAB>& a_denominator,
                 const LevelData<EBCellFAB>& a_fallback);

  /*!
    @brief Floor values in data holder. This sets all values below a_value to a_value
    @param[inout] a_lhs Data to be floored.
    @param[in]    a_value Floor value
    @note Also floors ghost cells. 
  */
  static void
  floor(EBAMRCellData& a_lhs, const Real a_value);

  /*!
    @brief Floor values in data holder. This sets all values below a_value to a_value
    @param[inout] a_lhs Data to be floored.
    @param[in]    a_value Floor value
    @note Also floors ghost cells. 
  */
  static void
  floor(LevelData<EBCellFAB>& a_lhs, const Real a_value);

  /*!
    @brief Floor values in data holder. This sets all values below a_value to a_value
    @param[inout] a_lhs Data to be floored.
    @param[in]    a_value Floor value
    @note Also floors ghost cells. 
  */
  static void
  floor(EBAMRIVData& a_lhs, const Real a_value);

  /*!
    @brief Floor values in data holder. This sets all values below a_value to a_value
    @param[inout] a_lhs Data to be floored.
    @param[in]    a_value Floor value
    @note Also floors ghost cells. 
  */
  static void
  floor(LevelData<BaseIVFAB<Real>>& a_lhs, const Real a_value);

  /*!
    @brief Get maximum and minimum value of specified component
    @param[out] a_max  Maximum value
    @param[out] a_min  Minium value
    @param[in]  a_data Cell-centered data
    @param[in]  a_comp Component
    @note This does the calculation over all levels, including grids that is covered by other grids. 
  */
  static void
  getMaxMin(Real& max, Real& min, EBAMRCellData& a_data, const int a_comp);

  /*!
    @brief Get maximum and minimum value of specified component
    @param[out] a_max  Maximum value
    @param[out] a_min  Minium value
    @param[in]  a_data Cell-centered data
    @param[in]  a_comp Component
  */
  static void
  getMaxMin(Real& max, Real& min, LevelData<EBCellFAB>& a_data, const int a_comp);

  /*!
    @brief Get maximum and minimum value. Assumes single-component data.
    @param[out] a_max  Maximum value
    @param[out] a_min  Minium value
    @param[in]  a_data Cell-centered data
    @note This does the calculation over all levels, including grids that is covered by other grids. 
  */
  static void
  getMaxMin(Vector<Real>& a_max, Vector<Real>& a_min, Vector<EBAMRCellData>& a_data);

  /*!
    @brief Get maximum and minimum value of normed data.
    @param[out] a_max  Maximum value
    @param[out] a_min  Minium value
    @param[in]  a_data Cell-centered data
    @note This does the calculation over all levels, including grids that is covered by other grids. 
  */
  static void
  getMaxMinNorm(Real& a_max, Real& a_min, EBAMRCellData& data);

  /*!
    @brief Get maximum and minimum value of normed data.
    @param[out] a_max  Maximum value
    @param[out] a_min  Minium value
    @param[in]  a_data Cell-centered data
  */
  static void
  getMaxMinNorm(Real& a_max, Real& a_min, LevelData<EBCellFAB>& data);

  /*!
    @brief Get maximum and minimum value of normed data.
    @param[out] a_max  Maximum value
    @param[out] a_min  Minium value
    @param[in]  a_data Cell-centered data
    @note This does the calculation over all levels, including grids that is covered by other grids. 
  */
  static void
  getMaxMinNorm(Real& a_max, Real& a_min, EBAMRIVData& data);

  /*!
    @brief Get maximum and minimum value of normed data.
    @param[out] a_max  Maximum value
    @param[out] a_min  Minium value
    @param[in]  a_data Cell-centered data
  */
  static void
  getMaxMinNorm(Real& a_max, Real& a_min, LevelData<BaseIVFAB<Real>>& data);

  /*!
    @brief Invert data
    @param[inout] Data to be inverted. On output we have a_data = 1/a_data.
    @note Does all components and does not guard against division by zero.
  */
  static void
  invert(EBAMRFluxData& a_data);

  /*!
    @brief Invert data
    @param[inout] Data to be inverted. On output we have a_data = 1/a_data.
    @note Does all components and does not guard against division by zero.
  */
  static void
  invert(LevelData<EBFluxFAB>& a_data);

  /*!
    @brief Invert data
    @param[inout] Data to be inverted. On output we have a_data = 1/a_data.
    @note Does all components and does not guard against division by zero.
  */
  static void
  invert(EBFluxFAB& a_data);

  /*!
    @brief Compute the kappa-weighted sum of the input data. Assumes one component
    @param[out] a_mass Total mass, computed as kappa * lhs
    @param[in]  a_lhs  Cell data
  */
  static void
  kappaSum(Real& a_mass, const LevelData<EBCellFAB>& a_lhs);

  /*!
    @brief Scale data by volume fraction
    @param[inout] a_data Data to be scaled. On output we have a_data = kappa * a_data
  */
  static void
  kappaScale(EBAMRCellData& a_data);

  /*!
    @brief Scale data by volume fraction
    @param[inout] a_data Data to be scaled. On output we have a_data = kappa * a_data
  */
  static void
  kappaScale(LevelData<EBCellFAB>& a_data);

  /*!
    @brief Scale data by volume fraction
    @param[inout] a_data Data to be scaled. On output we have a_data = kappa * a_data
  */
  static void
  kappaScale(MFAMRCellData& a_data);

  /*!
    @brief Scale data by volume fraction
    @param[inout] a_data Data to be scaled. On output we have a_data = kappa * a_data
  */
  static void
  kappaScale(LevelData<MFCellFAB>& a_data);

  /*!
    @brief Scale data by dx^SpaceDim
    @param[inout] a_data Data to be scaled. On output we have a_data = dx^SpaceDim * a_data
  */
  static void
  volumeScale(EBAMRCellData& a_data, const Vector<Real>& a_dx);

  /*!
    @brief Multiply data holder by another data holder
    @param[inout] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  /*!
    @brief Multiply data holder by another data holder
    @param[inout] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  /*!
    @brief Multiply data holder by another data holder
    @param[inout] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(EBAMRFluxData& a_lhs, const EBAMRFluxData& a_rhs);

  /*!
    @brief Multiply data holder by another data holder
    @param[inout] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(LevelData<EBFluxFAB>& a_lhs, const LevelData<EBFluxFAB>& a_rhs);

  /*!
    @brief Multiply data holder by another data holder
    @param[inout] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs);

  /*!
    @brief Multiply data holder by another data holder
    @param[inout] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note Does all components
  */
  static void
  multiply(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs);

  /*!
    @brief Multiply data holder by another data holder
    @param[inout] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note This multiplies all components in a_lhs with the single component in a_rhs
  */
  static void
  multiplyScalar(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  /*!
    @brief Multiply data holder by another data holder
    @param[inout] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note This multiplies all components in a_lhs with the single component in a_rhs
  */
  static void
  multiplyScalar(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  /*!
    @brief Multiply data holder by another data holder
    @param[inout] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note This multiplies all components in a_lhs with the single component in a_rhs
  */
  static void
  multiplyScalar(EBAMRIVData& a_lhs, const EBAMRIVData& a_rhs);

  /*!
    @brief Multiply data holder by another data holder
    @param[inout] a_lhs Data to be multiplied
    @param[in]    a_rhs Multiplication factor
    @note This multiplies all components in a_lhs with the single component in a_rhs
  */
  static void
  multiplyScalar(LevelData<BaseIVFAB<Real>>& a_lhs, const LevelData<BaseIVFAB<Real>>& a_rhs);

  /*!
    @brief Compute specified norm of data. Calls EBLevelDataOps
    @param[out] a_norm   Norm
    @param[in]  a_data   Data to compute the norm of
    @param[in]  a_domain Domain
    @param[in]  a_data   Norm factor, e.g. a_p = 2 => squared norm. 0 => inf norm
  */
  static void
  norm(Real& a_norm, const LevelData<EBCellFAB>& a_data, const ProblemDomain& a_domain, const int a_p);

  /*!
    @brief Compute specified norm of data. 
    @param[in]  a_data    Data to compute the norm of
    @param[in]  a_domain  Domain
    @param[in]  a_data    Norm factor, e.g. a_p = 2 => squared norm. 0 => inf norm
    @param[in]  a_regOnly If true, only include regular cells in norm computation. 
  */
  static Real
  norm(const LevelData<EBCellFAB>& a_data, const ProblemDomain& a_domain, const int a_p, const bool a_regOnly);

  /*!
    @brief General addition operator for adding together data. The user can choose which components to add. 
    @param[inout] a_lhs     Output data component
    @param[inout] a_rhs     Data to be added to a_lhs
    @param[in]    a_srcComp Component in a_rhs to add from.
    @param[in]    a_dstComp Component in a_lhs to add into.
    @param[in]    a_numComp Number of components to add. 
    @note This will add a_lhs(a_dstComp + icomp) += a_rhs(a_srcComp + icomp) for icomp = 0,1,2,..numComp-1.
  */
  static void
  plus(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs, const int a_srcComp, const int a_dstComp, const int a_numComp);

  /*!
    @brief General addition operator for adding together data. The user can choose which components to add. 
    @param[inout] a_lhs     Output data component
    @param[inout] a_rhs     Data to be added to a_lhs
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

  /*!
    @brief Scale data by factor
    @param[inout] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(MFAMRCellData& a_lhs, const Real& a_scale);

  /*!
    @brief Scale data by factor
    @param[inout] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(LevelData<MFCellFAB>& a_lhs, const Real& a_scale);

  /*!
    @brief Scale data by factor
    @param[inout] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(MFAMRFluxData& a_lhs, const Real& a_scale);

  /*!
    @brief Scale data by factor
    @param[inout] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(LevelData<MFFluxFAB>& a_lhs, const Real& a_scale);

  /*!
    @brief Scale data by factor
    @param[inout] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(EBAMRIVData& a_lhs, const Real& a_scale);

  /*!
    @brief Scale data by factor
    @param[inout] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(EBAMRCellData& a_lhs, const Real a_scale);

  /*!
    @brief Scale data by factor
    @param[inout] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(LevelData<EBCellFAB>& a_lhs, const Real a_scale);

  /*!
    @brief Scale data by factor
    @param[inout] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(EBAMRFluxData& a_lhs, const Real a_scale);

  /*!
    @brief Scale data by factor
    @param[inout] a_lhs   Data to be scaled
    @param[in]    a_scale Scaling factor
  */
  static void
  scale(LevelData<BaseIVFAB<Real>>& a_lhs, const Real& a_scale);

  /*!
    @brief Set value in covered cells. Does specified component
    @param[inout] a_lhs   Input/output data. 
    @param[ino]   a_comp  Component to set
    @param[ino]   a_value Value to set
  */
  static void
  setCoveredValue(EBAMRCellData& a_lhs, const int a_comp, const Real a_value);

  /*!
    @brief Set value in covered cells. Does specified component
    @param[inout] a_lhs   Input/output data. 
    @param[ino]   a_comp  Component to set
    @param[ino]   a_value Value to set
  */
  static void
  setCoveredValue(LevelData<EBCellFAB>& a_lhs, const int a_comp, const Real a_value);

  /*!
    @brief Version of setCoveredValue which does all component
    @param[inout] a_lhs   Data 
    @param[in]    a_value Value in covered grid cells. 
  */
  static void
  setCoveredValue(EBAMRCellData& a_lhs, const Real a_value);

  /*!
    @brief Version of setCoveredValue which does all component
    @param[inout] a_lhs   Data 
    @param[in]    a_value Value in covered grid cells. 
  */
  static void
  setCoveredValue(LevelData<EBCellFAB>& a_lhs, const Real a_value);

  /*!
    @brief Set value in cells that are covered by a finer grid.
    @param[inout] a_lhs    AMR data
    @param[in]    a_refRat Refinement ratios between levels.
    @param[in]    a_value  Value to set in cells covered by a finer grid. 
  */
  static void
  setInvalidValue(EBAMRCellData& a_lhs, const Vector<int>& a_refRat, const Real a_value);

  /*!
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
           const RealVect                             a_probLo,
           const Vector<Real>&                        a_dx,
           const int                                  a_comp);

  /*!
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

  /*!
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set
    @param[in]  a_function Function to use for setting the value. 
    @param[in]  a_probLo   Lower-left corner of physical domain. 
    @param[in]  a_dx       Grid resolutions
    @param[in]  a_comp     Component to set
    @note Uses a VofIterator everywhere to this might be slow!
  */
  static void
  setValue(EBAMRCellData&                             a_lhs,
           const std::function<Real(const RealVect)>& a_function,
           const RealVect                             a_probLo,
           const Vector<Real>&                        a_dx,
           const int                                  a_comp);

  /*!
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set
    @param[in]  a_function Function to use for setting the value. 
    @param[in]  a_probLo   Lower-left corner of physical domain. 
    @param[in]  a_dx       Grid resolutions
    @param[in]  a_comp     Component to set
    @note Uses a VofIterator everywhere to this might be slow!
  */
  static void
  setValue(LevelData<EBCellFAB>&                      a_lhs,
           const std::function<Real(const RealVect)>& a_function,
           const RealVect                             a_probLo,
           const Real                                 a_dx,
           const int                                  a_comp);

  /*!
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
           const RealVect                             a_probLo,
           const Vector<Real>&                        a_dx,
           const int                                  a_comp);

  /*!
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

  /*!
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set
    @param[in]  a_function Function to use for setting the value. 
    @param[in]  a_probLo   Lower-left corner of physical domain. 
    @param[in]  a_dx       Grid resolutions
    @param[in]  a_comp     Component to set
    @note Uses a VofIterator everywhere to this might be slow!
  */
  static void
  setValue(EBAMRIVData&                               a_lhs,
           const std::function<Real(const RealVect)>& a_function,
           const RealVect                             a_probLo,
           const Vector<Real>&                        a_dx,
           const int                                  a_comp);

  /*!
    @brief Polymorphic set value function. Takes a spatially varying function and sets the value in the specified component from that function.
    @param[out] a_lhs      Data to set
    @param[in]  a_function Function to use for setting the value. 
    @param[in]  a_probLo   Lower-left corner of physical domain. 
    @param[in]  a_dx       Grid resolutions
    @param[in]  a_comp     Component to set
    @note Uses a VofIterator everywhere to this might be slow!
  */
  static void
  setValue(LevelData<BaseIVFAB<Real>>&                a_lhs,
           const std::function<Real(const RealVect)>& a_function,
           const RealVect                             a_probLo,
           const Real                                 a_dx,
           const int                                  a_comp);

  /*!
    @brief Polymorphic set value function. Assumes that a_lhs has SpaceDim components and sets all those components from the input function.
    @param[out] a_lhs      Data to set
    @param[in]  a_function Function to use for setting the value. 
    @param[in]  a_probLo   Lower-left corner of physical domain. 
    @param[in]  a_dx       Grid resolutions
    @param[in]  a_comp     Component to set
    @note Uses a VofIterator everywhere to this might be slow!
  */
  static void
  setValue(EBAMRCellData&                                 a_lhs,
           const std::function<RealVect(const RealVect)>& a_function,
           const RealVect                                 a_probLo,
           const Vector<Real>&                            a_dx);

  /*!
    @brief Polymorphic set value function. Assumes that a_lhs has SpaceDim components and sets all those components from the input function.
    @param[out] a_lhs      Data to set
    @param[in]  a_function Function to use for setting the value. 
    @param[in]  a_probLo   Lower-left corner of physical domain. 
    @param[in]  a_dx       Grid resolutions
    @param[in]  a_comp     Component to set
    @note Uses a VofIterator everywhere to this might be slow!
  */
  static void
  setValue(LevelData<EBCellFAB>&                          a_lhs,
           const std::function<RealVect(const RealVect)>& a_function,
           const RealVect                                 a_probLo,
           const Real                                     a_dx);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(EBAMRCellData& a_lhs, const Real& a_value);

  /*!
    @brief Set value in data with specified component. 
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(EBAMRCellData& a_lhs, const Real a_value, const int a_comp);

  /*!
    @brief Set value in data with specified component. 
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<EBCellFAB>& a_lhs, const Real a_value, const int a_comp);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<EBCellFAB>& a_lhs, const Real a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<EBFluxFAB>& a_lhs, const Real a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<BaseIVFAB<Real>>& a_lhs, const Real a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(EBAMRFluxData& a_lhs, const Real& a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(EBAMRIVData& a_lhs, const Real& a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(MFAMRCellData& a_lhs, const Real& a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<MFCellFAB>& a_lhs, const Real& a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(MFAMRFluxData& a_lhs, const Real& a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<MFFluxFAB>& a_lhs, const Real& a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(MFAMRIVData& a_lhs, const Real& a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<MFBaseIVFAB>& a_lhs, const Real& a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(EBAMRIFData& a_lhs, const Real a_value);

  /*!
    @brief Set value in data
    @param[out] a_lhs   Input/output data
    @param[in]  a_value Value to set
  */
  static void
  setValue(LevelData<DomainFluxIFFAB>& a_lhs, const Real a_value);

  /*!
    @brief Compute the sum of the input value.
    @param[inout] a_value On output, contains the sum of a_value over all MPI processes.
  */
  static void
  sum(Real& a_value);

  /*!
    @brief Compute the square root of the input data
    @param[inout] a_lhs On output, contains sqrt(a_lhs)
    @note Does not do ghosted faces. 
  */
  static void
  squareRoot(EBAMRFluxData& a_lhs);

  /*!
    @brief Compute the square root of the input data
    @param[inout] a_lhs On output, contains sqrt(a_lhs)
    @note Does not do ghosted faces. 
  */
  static void
  squareRoot(LevelData<EBFluxFAB>& a_lhs);

  /*!
    @brief Compute the vector length of a data holder. Sets a_lhs = |a_rhs| where a_rhs contains SpaceDim components. 
    @param[out] a_lhs Vector length
    @param[in]  a_rhs Input vector
    @note a_lhs must have one component and a_rhs must have SpaceDim components. 
  */
  static void
  vectorLength(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  /*!
    @brief Compute the vector length of a data holder. Sets a_lhs = |a_rhs| where a_rhs contains SpaceDim components. 
    @param[out] a_lhs Vector length
    @param[in]  a_rhs Input vector
    @note a_lhs must have one component and a_rhs must have SpaceDim components. 
  */
  static void
  vectorLength(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  /*!
    @brief Compute the vector length of a data holder. Sets a_lhs = |a_rhs| where a_rhs contains SpaceDim components. 
    @param[out] a_lhs Vector length.
    @param[in]  a_rhs Input vector.
    @param[in]  a_box Computation box.
    @note a_lhs must have one component and a_rhs must have SpaceDim components. 
  */
  static void
  vectorLength(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& box);

  /*!
    @brief Compute the squared vector length of a data holder. Sets a_lhs = |a_rhs|^2 where a_rhs contains SpaceDim components. 
    @param[out] a_lhs Vector length
    @param[in]  a_rhs Input vector
    @note a_lhs must have one component and a_rhs must have SpaceDim components. 
  */
  static void
  vectorLength2(EBAMRCellData& a_lhs, const EBAMRCellData& a_rhs);

  /*!
    @brief Compute the squared vector length of a data holder. Sets a_lhs = |a_rhs|^2 where a_rhs contains SpaceDim components. 
    @param[out] a_lhs Vector length
    @param[in]  a_rhs Input vector
    @note a_lhs must have one component and a_rhs must have SpaceDim components. 
  */
  static void
  vectorLength2(LevelData<EBCellFAB>& a_lhs, const LevelData<EBCellFAB>& a_rhs);

  /*!
    @brief Compute the squared vector length of a data holder. Sets a_lhs = |a_rhs|^2 where a_rhs contains SpaceDim components. 
    @param[out] a_lhs Vector length
    @param[in]  a_rhs Input vector
    @note a_lhs must have one component and a_rhs must have SpaceDim components. 
  */
  static void
  vectorLength2(EBCellFAB& a_lhs, const EBCellFAB& a_rhs, const Box& box);

  /*!
    @brief Compute the tightest possible valid box around a cut-cell volume
    @param[out] a_lo       Low corner
    @param[out] a_hi       High corner
    @param[in]  a_normal   EB face normal
    @param[in]  a_centroid EB centroid
  */
  static void
  computeMinValidBox(RealVect& a_lo, RealVect& a_hi, const RealVect a_normal, const RealVect a_centroid);

  /*!
    @brief Check if all corners/points are inside an EB
    @param[out] a_corner   Corners
    @param[in]  a_normal   EB face normal
    @param[in]  a_centroid EB centroid
    @return Returns true if all points are inside the EB and false otherwise. 
  */
  static bool
  allCornersInsideEb(const Vector<RealVect>& a_corners, const RealVect a_normal, const RealVect a_centroid);

  /*!
    @brief Shift pointes by a specified distance
    @details On output, every point in a_corners is shifted by a_distance. 
    @param[inout] a_corners  Points to be shifted
    @param[in]    a_distance Shift vector
  */
  static void
  shiftCorners(Vector<RealVect>& a_corners, const RealVect& a_distance);
};

#include <CD_NamespaceFooter.H>

#include <CD_DataOpsImplem.H>

#endif
