/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_Tile.cpp
  @brief  Simple class for representing a tile in arbitrary-dimensional space.
  @author Robert Marskar
*/

#ifndef CD_Tile_H
#define CD_Tile_H

// Std includes
#include <array>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for representing a tile (used in, e.g., TiledMeshRefine)
*/
template <class T, size_t DIM>
class TileI
{
public:
  /*!
    @brief Default constructor. Sets m_indices = 0
  */
  TileI() noexcept;

  /*!
    @brief Full constructor. Sets indices. Use as TileI(0,1,2,...)
  */
  template <typename... Args, typename std::enable_if<sizeof...(Args) == DIM, int>::type = 0>
  explicit TileI(Args... args) noexcept : m_indices{static_cast<T>(args)...}
  {}

  /*!
    @brief Destructor. 
  */
  virtual ~TileI() noexcept;

  /*!
    @brief Get the index in specified coordinate direction
    @param[in] a_dir Coordinate direction
  */
  inline T&
  operator[](const size_t a_dir) noexcept;

  /*!
    @brief Get the index in specified coordinate direction
    @param[in] a_dir Coordinate direction
  */
  inline const T&
  operator[](const size_t a_dir) const noexcept;

  /*!
    @brief Get the indices
    @return m_indices
  */
  inline std::array<T, DIM>&
  getIndices() noexcept;

  /*!
    @brief Get the indices
    @return m_indices
  */
  inline const std::array<T, DIM>&
  getIndices() const noexcept;

  /*!
    @brief Comparison operator. Uses lexicographical comparison.
    @param[in] a_otherTile Other tile.
  */
  inline constexpr bool
  operator<(const TileI& a_otherTile) const noexcept;

  /*!
    @brief Comparison operator. Uses lexicographical comparison.
    @param[in] a_otherTile Other tile.
  */
  inline constexpr bool
  operator>(const TileI& a_otherTile) const noexcept;

  /*!
    @brief Comparison operator. Uses lexicographical comparison.
    @param[in] a_otherTile Other tile.
  */
  inline constexpr bool
  operator==(const TileI& a_otherTile) const noexcept;

protected:
  /*!
    @brief Index in N-dimensional space
  */
  std::array<T, DIM> m_indices;
};

#include <CD_NamespaceFooter.H>

#include <CD_TileImplem.H>

#endif
