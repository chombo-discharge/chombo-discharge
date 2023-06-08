/* chombo-discharge
 * Copyright © 2023 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_TileImplem.H
  @brief  Implementation of CD_Tile.H
  @author Robert Marskar
*/

#ifndef CD_TileImplem_H
#define CD_TileImplem_H

// Std includes
#include <utility>

// Our includes
#include <CD_Tile.H>
#include <CD_NamespaceHeader.H>

// Utility functions, but std::tuple_cat should work.
// template <class T, size_t DIM, size_t... Is>
// auto
// unpack_impl(const std::array<T, DIM>& arr, std::index_sequence<Is...>) -> decltype(std::make_tuple(arr[Is]...))
// {
//   return std::make_tuple(arr[Is]...);
// }

// template <class T, size_t DIM>
// auto
// unpack(const std::array<T, DIM>& arr) -> decltype(unpack_impl(arr, std::make_index_sequence<DIM>{}))
// {
//   return unpack_impl(arr, std::make_index_sequence<DIM>{});
// }

template <class T, size_t DIM>
TileI<T, DIM>::TileI() noexcept
{
  for (auto& index : m_indices) {
    index = 0;
  }
}

// template <class T, size_t DIM>
// template <typename... Args>
// TileI<T, DIM>::TileI(Args... args) noexcept
// {}

template <class T, size_t DIM>
TileI<T, DIM>::~TileI() noexcept
{}

template <class T, size_t DIM>
inline T&
TileI<T, DIM>::operator[](const size_t a_dir) noexcept
{
  return m_indices[a_dir];
}

template <class T, size_t DIM>
inline const T&
TileI<T, DIM>::operator[](const size_t a_dir) const noexcept
{
  return m_indices[a_dir];
}

template <class T, size_t DIM>
inline std::array<T, DIM>&
TileI<T, DIM>::getIndices() noexcept
{
  return m_indices;
}

template <class T, size_t DIM>
inline const std::array<T, DIM>&
TileI<T, DIM>::getIndices() const noexcept
{
  return m_indices;
}

template <class T, size_t DIM>
inline constexpr bool
TileI<T, DIM>::operator<(const TileI& a_otherTile) const noexcept
{
  // const auto& myIndices    = unpack<T, DIM>(this->getIndices());
  // const auto& otherIndices = unpack<T, DIM>(a_otherTile.getIndices());

  const auto& myIndices    = std::tuple_cat(m_indices);
  const auto& otherIndices = std::tuple_cat(a_otherTile.getIndices());

  return std::tie(myIndices) < std::tie(otherIndices);
}

template <class T, size_t DIM>
inline constexpr bool
TileI<T, DIM>::operator>(const TileI& a_otherTile) const noexcept
{
  //  const auto& myIndices    = unpack<T, DIM>(this->getIndices());
  //  const auto& otherIndices = unpack<T, DIM>(a_otherTile.getIndices());
  const auto& myIndices    = std::tuple_cat(m_indices);
  const auto& otherIndices = std::tuple_cat(a_otherTile.getIndices());

  return std::tie(myIndices) > std::tie(otherIndices);
}

template <class T, size_t DIM>
inline constexpr bool
TileI<T, DIM>::operator==(const TileI& a_otherTile) const noexcept
{
  // const auto& myIndices    = unpack<T, DIM>(this->getIndices());
  // const auto& otherIndices = unpack<T, DIM>(a_otherTile.getIndices());

  const auto& myIndices    = std::tuple_cat(m_indices);
  const auto& otherIndices = std::tuple_cat(a_otherTile.getIndices());

  return std::tie(myIndices) == std::tie(otherIndices);
}

#include <CD_NamespaceFooter.H>

#endif
