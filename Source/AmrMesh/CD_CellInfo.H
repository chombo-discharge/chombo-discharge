/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_CellInfo.H
  @brief  Simple class for holding some quantities relevant in a grid cell. 
  @author Robert Marskar
*/

#ifndef CD_CELLINFO_H
#define CD_CELLINFO_H

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/**
  @brief Class for the cell-information that is often queried when merging particles inside a cell. 
*/
class CellInfo
{
public:
  /**
    @brief Default constructor - creates a regular cell
  */
  CellInfo() = delete;

  /**
    @brief Constructor which creates a regular cell 
    @param[in] a_gridIndex     Grid index
    @param[in] a_dx            Grid resolution
  */
  CellInfo(const IntVect a_gridIndex, const Real a_dx) noexcept;

  /**
    @brief Constructor which automatically constructs an approximation to the minimum box enclosing the valid region of the cell
    @param[in] a_gridIndex     Grid index
    @param[in] a_dx            Grid resolution
    @param[in] a_volFrac       Volume fraction. 
    @param[in] a_bndryCentroid Boundary centroid
    @param[in] a_bndryNormal   Boundary normal
  */
  CellInfo(const IntVect   a_gridIndex,
           const Real      a_dx,
           const Real      a_volFrac,
           const RealVect& a_bndryCentroid,
           const RealVect& a_bndryNormal) noexcept;

  /**
    @brief Destructor
  */
  virtual ~CellInfo() noexcept;

  /**
    @brief Get the grid index
    @return Grid index
  */
  IntVect&
  getGridIndex() noexcept;

  /**
    @brief Get the grid index
    @return Grid index
  */
  const IntVect&
  getGridIndex() const noexcept;

  /**
    @brief Get the grid spacing
    @return Dx
  */
  Real&
  getDx() noexcept;

  /**
    @brief Get the grid spacing
    @return Dx
  */
  const Real&
  getDx() const noexcept;

  /**
    @brief Get the volume fraction
    @return Vol frac
  */
  Real&
  getVolFrac() noexcept;

  /**
    @brief Get the volume fraction
    @return Vol frac
  */
  const Real&
  getVolFrac() const noexcept;

  /**
    @brief Get the boundary centroid
    @return Bndry centroid
  */
  RealVect&
  getBndryCentroid() noexcept;

  /**
    @brief Get the boundary centroid
    @return Bndry centroid
  */
  const RealVect&
  getBndryCentroid() const noexcept;

  /**
    @brief Get the boundary normal
    @return Bndry normal
  */
  RealVect&
  getBndryNormal() noexcept;

  /**
    @brief Get the boundary normal
    @return Bndry normal
  */
  const RealVect&
  getBndryNormal() const noexcept;

  /**
    @brief Get the lower valid region
    @return Valid lo
  */
  RealVect&
  getValidLo() noexcept;

  /**
    @brief Get the lower valid region
    @return Valid lo
  */
  const RealVect&
  getValidLo() const noexcept;

  /**
    @brief Get the upper valid region
    @return Valid hi
  */
  RealVect&
  getValidHi() noexcept;

  /**
    @brief Get the upper valid region
    @return Valid hi
  */
  const RealVect&
  getValidHi() const noexcept;

protected:
  /**
    @brief Grid index
  */
  IntVect m_gridIndex;

  /**
    @brief Grid spacing
  */
  Real m_dx;

  /**
    @brief Volume fraction
  */
  Real m_volFrac;

  /**
    @brief Boundary centroid
  */
  RealVect m_bndryCentroid;

  /**
    @brief EB boundary normal
  */
  RealVect m_bndryNormal;

  /**
    @brief Valid region in each direction (lower limit)
  */
  RealVect m_validLo;

  /**
    @brief Valid region in each direction (upper limit)
  */
  RealVect m_validHi;
};

#include <CD_NamespaceFooter.H>

#endif
