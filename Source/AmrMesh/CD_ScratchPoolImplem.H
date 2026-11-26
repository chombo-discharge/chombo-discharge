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
Scratch<T>::operator T&() noexcept
{
  Scratch<T>::assertOwner(m_pool);

  return m_pool->template get<T>(m_entry);
}

template <typename T>
Scratch<T>::operator const T&() const noexcept
{
  Scratch<T>::assertOwner(m_pool);

  return m_pool->template get<T>(m_entry);
}

template <typename T>
T&
Scratch<T>::operator*() noexcept
{
  Scratch<T>::assertOwner(m_pool);

  return m_pool->template get<T>(m_entry);
}

template <typename T>
const T&
Scratch<T>::operator*() const noexcept
{
  Scratch<T>::assertOwner(m_pool);

  return m_pool->template get<T>(m_entry);
}

template <typename T>
T*
Scratch<T>::operator->() noexcept
{
  Scratch<T>::assertOwner(m_pool);

  return &(m_pool->template get<T>(m_entry));
}

template <typename T>
const T*
Scratch<T>::operator->() const noexcept
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
                     const Build&             a_build) noexcept
{
  CH_TIME("ScratchPool::acquire");

  ScratchPool::assertSerial("ScratchPool::acquire");

  std::vector<Entry<T>>& list = this->entries<T>();

  const int numEntries = static_cast<int>(list.size());

  // A free entry with the same phase, level and number of components, and at least the requested number of ghost
  // cells, will do. Prefer the fewest ghost cells, so that an exact match is used when there is one and the entries
  // with wide ghost regions stay free for the callers that need them.
  int best = -1;

  for (int i = 0; i < numEntries; i++) {
    const Entry<T>& entry = list[i];

    const bool match = !entry.inUse && (entry.phase == a_phase) && (entry.level == a_level) &&
                       (entry.nComp == a_nComp) && (entry.ghost >= a_ghost);

    if (match && (best < 0 || entry.ghost < list[best].ghost)) {
      best = i;
    }
  }

  if (best >= 0) {
    list[best].inUse = true;

    return Scratch<T>(*this, best);
  }

  Entry<T> entry;

  entry.data  = a_build();
  entry.phase = a_phase;
  entry.level = a_level;
  entry.nComp = a_nComp;
  entry.ghost = a_ghost;
  entry.inUse = true;

  if (!entry.data) {
    MayDay::Abort("ScratchPool::acquire -- the build function returned no data");
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
  CH_assert(list[a_entry].inUse);

  return *(list[a_entry].data);
}

template <typename T>
void
ScratchPool::giveBack(const int a_entry) noexcept
{
  CH_TIME("ScratchPool::giveBack");

  ScratchPool::assertSerial("ScratchPool::giveBack");

  std::vector<Entry<T>>& list = this->entries<T>();

  if (a_entry < 0 || a_entry >= static_cast<int>(list.size()) || !(list[a_entry].inUse)) {
    MayDay::Abort("ScratchPool::giveBack -- returning a temporary that is not in use");
  }

  // A level referenced from outside the pool means that a copy of the EBAMRData outlived its handle. Handing the
  // entry out again would then give two owners the same data.
  if constexpr (!std::is_same<T, LevelData<EBCellFAB>>::value) {
    const T& data = *(list[a_entry].data);

    for (int lvl = 0; lvl < data.size(); lvl++) {
      if (data[lvl].refCount() > 1) {
        const std::string err = std::string("ScratchPool::giveBack -- a copy of a ") + ScratchPool::typeName<T>() +
                                " temporary outlived its Scratch handle";

        MayDay::Abort(err.c_str());
      }
    }
  }

  list[a_entry].inUse = false;
}

template <typename T>
void
ScratchPool::reportType(std::ostream& a_os, const std::string& a_prefix) const noexcept
{
  const std::vector<Entry<T>>& list = this->entries<T>();

  for (const Entry<T>& entry : list) {
    a_os << a_prefix << ScratchPool::typeName<T>() << ": phase = " << static_cast<int>(entry.phase)
         << ", level = " << entry.level << ", components = " << entry.nComp << ", ghost cells = " << entry.ghost
         << (entry.inUse ? ", in use" : ", free") << "\n";
  }
}

template <typename T>
int
ScratchPool::numInUse() const noexcept
{
  int count = 0;

  for (const Entry<T>& entry : this->entries<T>()) {
    if (entry.inUse) {
      count++;
    }
  }

  return count;
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
