/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_ScratchPoolImplem.H
 * @brief  Implementation of CD_ScratchPool.H
 * @author Robert Marskar
 */

#ifndef CD_SCRATCHPOOLIMPLEM_H
#define CD_SCRATCHPOOLIMPLEM_H

// Std includes
#include <type_traits>
#include <utility>

// Chombo includes
#include <CH_Timer.H>
#include <MayDay.H>

// Our includes
#include <CD_ScratchPool.H>
#include <CD_NamespaceHeader.H>

template <typename T>
Scratch<T>::Scratch(ScratchPool& a_pool, const int a_entry) noexcept : m_pool(&a_pool), m_entry(a_entry)
{}

template <typename T>
Scratch<T>::Scratch(Scratch&& a_other) noexcept : m_pool(a_other.m_pool), m_entry(a_other.m_entry)
{
  a_other.m_pool  = nullptr;
  a_other.m_entry = -1;
}

template <typename T>
void
Scratch<T>::assertOwner(const ScratchPool* const a_pool) noexcept
{
  if (a_pool == nullptr) {
    MayDay::Abort("Scratch -- using a handle that was moved from");
  }
}

template <typename T>
Scratch<T>::~Scratch() noexcept
{
  if (m_pool != nullptr) {
    m_pool->template giveBack<T>(m_entry);
  }
}

template <typename T>
Scratch<T>::operator T&() & noexcept
{
  Scratch<T>::assertOwner(m_pool);

  return m_pool->template get<T>(m_entry);
}

template <typename T>
Scratch<T>::operator const T&() const& noexcept
{
  Scratch<T>::assertOwner(m_pool);

  return m_pool->template get<T>(m_entry);
}

template <typename T>
T&
Scratch<T>::operator*() & noexcept
{
  Scratch<T>::assertOwner(m_pool);

  return m_pool->template get<T>(m_entry);
}

template <typename T>
const T&
Scratch<T>::operator*() const& noexcept
{
  Scratch<T>::assertOwner(m_pool);

  return m_pool->template get<T>(m_entry);
}

template <typename T>
T*
Scratch<T>::operator->() & noexcept
{
  Scratch<T>::assertOwner(m_pool);

  return &(m_pool->template get<T>(m_entry));
}

template <typename T>
const T*
Scratch<T>::operator->() const& noexcept
{
  Scratch<T>::assertOwner(m_pool);

  return &(m_pool->template get<T>(m_entry));
}

template <typename T, typename Build>
Scratch<T>
ScratchPool::acquire(const phase::which_phase a_phase,
                     const int                a_level,
                     const int                a_nComp,
                     const int                a_ghost,
                     const Build&             a_build)
{
  CH_TIME("ScratchPool::acquire");

  ScratchPool::assertSerial("ScratchPool::acquire");

  std::vector<Entry<T>>& list = this->entries<T>();

  const int numEntries = static_cast<int>(list.size());

  // A free entry with exactly the same key will do. The number of ghost cells is matched exactly: the irregular-cell
  // sets of EB data grow with it, so a wider temporary is not interchangeable with a narrower one.
  for (int i = 0; i < numEntries; i++) {
    Entry<T>& entry = list[i];

    if (!entry.m_inUse && (entry.m_phase == a_phase) && (entry.m_level == a_level) && (entry.m_nComp == a_nComp) &&
        (entry.m_ghost == a_ghost)) {
      entry.m_inUse = true;

      return Scratch<T>(*this, i);
    }
  }

  Entry<T> entry;

  entry.m_data  = a_build();
  entry.m_phase = a_phase;
  entry.m_level = a_level;
  entry.m_nComp = a_nComp;
  entry.m_ghost = a_ghost;
  entry.m_inUse = true;

  if (!entry.m_data) {
    MayDay::Abort("ScratchPool::acquire -- the build function returned no data");
  }

  ScratchPool::forEachLevel(*entry.m_data,
                            [&entry](const int /*a_lvl*/,
                                     const void*              a_address,
                                     const DisjointBoxLayout& a_grids,
                                     const int                a_levelComps,
                                     const IntVect&           a_levelGhost) {
                              entry.m_layout.emplace_back(a_address, a_grids, a_levelComps, a_levelGhost);
                            });

  // The temporary must have the key it is stored under, or it would be handed out later as something else.
  for (const auto& level : entry.m_layout) {
    if (std::get<0>(level) != nullptr &&
        (std::get<2>(level) != a_nComp || std::get<3>(level) != a_ghost * IntVect::Unit)) {
      MayDay::Abort("ScratchPool::acquire -- the build function returned data that does not match the key");
    }
  }

  list.emplace_back(std::move(entry));

  return Scratch<T>(*this, numEntries);
}

template <typename T>
std::vector<ScratchPool::Entry<T>>&
ScratchPool::entries() noexcept
{
  return std::get<std::vector<Entry<T>>>(m_entries);
}

template <typename T>
const std::vector<ScratchPool::Entry<T>>&
ScratchPool::entries() const noexcept
{
  return std::get<std::vector<Entry<T>>>(m_entries);
}

template <typename T>
T&
ScratchPool::get(const int a_entry) noexcept
{
  std::vector<Entry<T>>& list = this->entries<T>();

  CH_assert(a_entry >= 0 && a_entry < static_cast<int>(list.size()));
  CH_assert(list[a_entry].m_inUse);

  return *(list[a_entry].m_data);
}

template <typename T>
void
ScratchPool::giveBack(const int a_entry) noexcept
{
  CH_TIME("ScratchPool::giveBack");

  ScratchPool::assertSerial("ScratchPool::giveBack");

  std::vector<Entry<T>>& list = this->entries<T>();

  if (a_entry < 0 || a_entry >= static_cast<int>(list.size()) || !(list[a_entry].m_inUse)) {
    MayDay::Abort("ScratchPool::giveBack -- returning a temporary that is not in use");
  }

  // A level referenced from outside the pool means that a copy of the EBAMRData outlived its handle. Handing the
  // entry out again would then give two owners the same data.
  if constexpr (!std::is_same<T, LevelData<EBCellFAB>>::value) {
    const T& data = *(list[a_entry].m_data);

    for (int lvl = 0; lvl < data.size(); lvl++) {
      if (data[lvl].refCount() > 1) {
        const std::string err = std::string("ScratchPool::giveBack -- a copy of a ") + ScratchPool::typeName<T>() +
                                " temporary outlived its Scratch handle";

        MayDay::Abort(err.c_str());
      }
    }
  }

  // A temporary that was reallocated, deallocated or given other levels through its handle would be handed out again
  // under a key it no longer matches.
  bool sameLayout = true;

  ScratchPool::forEachLevel(*(list[a_entry].m_data),
                            [&](const int                a_lvl,
                                const void*              a_address,
                                const DisjointBoxLayout& a_grids,
                                const int                a_levelComps,
                                const IntVect&           a_levelGhost) {
                              const auto& recorded = list[a_entry].m_layout;

                              sameLayout = sameLayout && (a_lvl < static_cast<int>(recorded.size())) &&
                                           (std::get<0>(recorded[a_lvl]) == a_address) &&
                                           (a_address == nullptr || std::get<1>(recorded[a_lvl]) == a_grids) &&
                                           (std::get<2>(recorded[a_lvl]) == a_levelComps) &&
                                           (std::get<3>(recorded[a_lvl]) == a_levelGhost);
                            });

  if (!sameLayout ||
      static_cast<int>(list[a_entry].m_layout.size()) != ScratchPool::numLevels(*(list[a_entry].m_data))) {
    const std::string err = std::string("ScratchPool::giveBack -- a ") + ScratchPool::typeName<T>() +
                            " temporary was reallocated, redefined or given other levels through its handle";

    MayDay::Abort(err.c_str());
  }

  list[a_entry].m_inUse = false;
}

template <typename T>
void
ScratchPool::reportType(const std::vector<Entry<T>>& a_list,
                        std::ostream&                a_os,
                        const std::string&           a_prefix) const noexcept
{
  for (const Entry<T>& entry : a_list) {
    a_os << a_prefix << ScratchPool::typeName<T>() << ": phase = " << static_cast<int>(entry.m_phase)
         << ", level = " << entry.m_level << ", components = " << entry.m_nComp << ", ghost cells = " << entry.m_ghost
         << (entry.m_inUse ? ", in use" : ", free") << "\n";
  }
}

template <typename T, typename Visit>
void
ScratchPool::forEachLevel(const T& a_data, const Visit& a_visit) noexcept
{
  if constexpr (std::is_same<T, LevelData<EBCellFAB>>::value) {
    a_visit(0, static_cast<const void*>(&a_data), a_data.disjointBoxLayout(), a_data.nComp(), a_data.ghostVect());
  }
  else {
    for (int lvl = 0; lvl < a_data.size(); lvl++) {
      if (a_data[lvl].isNull()) {
        a_visit(lvl, nullptr, DisjointBoxLayout(), -1, IntVect::Zero);
      }
      else {
        a_visit(lvl,
                static_cast<const void*>(&(*a_data[lvl])),
                a_data[lvl]->disjointBoxLayout(),
                a_data[lvl]->nComp(),
                a_data[lvl]->ghostVect());
      }
    }
  }
}

template <typename T>
int
ScratchPool::numLevels(const T& a_data) noexcept
{
  if constexpr (std::is_same<T, LevelData<EBCellFAB>>::value) {
    return 1;
  }
  else {
    return a_data.size();
  }
}

template <>
inline const char*
ScratchPool::typeName<EBAMRCellData>() noexcept
{
  return "EBAMRCellData";
}

template <>
inline const char*
ScratchPool::typeName<EBAMRFluxData>() noexcept
{
  return "EBAMRFluxData";
}

template <>
inline const char*
ScratchPool::typeName<EBAMRIVData>() noexcept
{
  return "EBAMRIVData";
}

template <>
inline const char*
ScratchPool::typeName<EBAMRIFData>() noexcept
{
  return "EBAMRIFData";
}

template <>
inline const char*
ScratchPool::typeName<EBAMRFAB>() noexcept
{
  return "EBAMRFAB";
}

template <>
inline const char*
ScratchPool::typeName<MFAMRCellData>() noexcept
{
  return "MFAMRCellData";
}

template <>
inline const char*
ScratchPool::typeName<MFAMRFluxData>() noexcept
{
  return "MFAMRFluxData";
}

template <>
inline const char*
ScratchPool::typeName<MFAMRIVData>() noexcept
{
  return "MFAMRIVData";
}

template <>
inline const char*
ScratchPool::typeName<LevelData<EBCellFAB>>() noexcept
{
  return "LevelData<EBCellFAB>";
}

#include <CD_NamespaceFooter.H>

#endif
