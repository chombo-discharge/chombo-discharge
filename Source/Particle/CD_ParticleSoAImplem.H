/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ParticleSoAImplem.H
  @brief  Implementation of CD_ParticleSoA.H
  @author Robert Marskar
*/

#ifndef CD_PARTICLESOAIMPLEM_H
#define CD_PARTICLESOAIMPLEM_H

// Std includes
#include <cmath>
#include <cstring>

// Chombo includes
#include <MayDay.H>

// Our includes
#include <CD_ParticleSoA.H>
#include <CD_NamespaceHeader.H>

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::reallocate(const std::size_t a_capacity)
{
  std::array<std::size_t, s_numColumns> off{};
  const std::size_t total  = this->computeOffsets(a_capacity, off, std::make_index_sequence<s_numColumns>{});
  void* const       newBuf = std::aligned_alloc(s_align, alignUp(total));
  if (newBuf == nullptr) {
    MayDay::Abort("ParticleSoA::reallocate - std::aligned_alloc returned null (out of memory)");
  }
  if (m_size > 0) {
    this->moveColumns(newBuf, off, std::make_index_sequence<s_numColumns>{});
  }
  std::free(m_buffer);
  m_buffer   = newBuf;
  m_capacity = a_capacity;
  this->setColPtrs(off, std::make_index_sequence<s_numColumns>{});
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::permute(const std::vector<std::size_t>& a_dest)
{
  std::array<std::size_t, s_numColumns> off{};
  const std::size_t total  = this->computeOffsets(m_capacity, off, std::make_index_sequence<s_numColumns>{});
  void* const       newBuf = std::aligned_alloc(s_align, alignUp(total));
  if (newBuf == nullptr) {
    MayDay::Abort("ParticleSoA::permute - std::aligned_alloc returned null (out of memory)");
  }
  this->scatterColumns(newBuf, off, a_dest, std::make_index_sequence<s_numColumns>{});
  std::free(m_buffer);
  m_buffer = newBuf;
  this->setColPtrs(off, std::make_index_sequence<s_numColumns>{});
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::reserve(const std::size_t a_capacity)
{
  if (a_capacity <= m_capacity) {
    return;
  }
  this->reallocate(a_capacity);
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::resize(const std::size_t a_count)
{
  this->reserve(a_count);
  m_size   = a_count;
  m_sorted = false;
}

template <typename P, typename Traits>
inline ParticleSoA<P, Traits>
ParticleSoA<P, Traits>::deepCopy() const
{
  ParticleSoA<P, Traits> dst;
  dst.reserve(m_capacity); // identical arena layout: reallocate() sets the capacity exactly
  if (m_size > 0) {
    std::memcpy(dst.m_buffer, m_buffer, this->byteSpan()); // one contiguous copy of the live region
  }
  dst.m_size      = m_size;
  dst.m_sorted    = m_sorted;
  dst.m_cellStart = m_cellStart;
  return dst;
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::deepCopyTo(ParticleSoA<P, Traits>& a_dst) const
{
  if (this == &a_dst) {
    return; // copying onto self is a no-op
  }

  a_dst.clear();         // drop a_dst's particles (keeps its arena) so the grow below moves nothing
  a_dst.reserve(m_size); // reuse a_dst's allocation if it is already large enough; grow only if not
  if (m_size > 0) {
    this->copyLiveColumnsTo(a_dst, std::make_index_sequence<s_numColumns>{});
  }
  a_dst.m_size      = m_size;
  a_dst.m_sorted    = m_sorted;
  a_dst.m_cellStart = m_cellStart;
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::append(const ParticleSoA<P, Traits>& a_other)
{
  const std::size_t n = a_other.m_size;
  if (n == 0) {
    return;
  }
  const std::size_t oldSize = m_size;
  this->reserve(oldSize + n); // may reallocate this (and a_other, when a_other == this); pointers refreshed below
  this->appendColumnsFrom(a_other, oldSize, std::make_index_sequence<s_numColumns>{});
  m_size   = oldSize + n;
  m_sorted = false;
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::appendParticle(const ParticleSoA<P, Traits>& a_src, const std::size_t a_index)
{
  CH_assert(a_index < a_src.m_size);

  const std::size_t at = m_size;
  this->reserve(m_size + 1); // may reallocate this (and a_src, when a_src == this); pointers refreshed below
  this->copyParticleColumns(a_src, a_index, at, std::make_index_sequence<s_numColumns>{});
  m_size   = at + 1;
  m_sorted = false;
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::shrinkToFit()
{
  if (m_capacity <= m_size) {
    return; // already compact
  }
  if (m_size == 0) {
    std::free(m_buffer);
    m_buffer   = nullptr;
    m_capacity = 0;
    m_colPtrs.fill(nullptr);
    return;
  }
  this->reallocate(m_size);
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::append(const RealVect& a_position, const double a_weight, const P& a_payload)
{
  this->growIfFull();
  for (int dir = 0; dir < SpaceDim; dir++) {
    this->positionColumn(dir)[m_size] = a_position[dir]; // write the reserved slot directly (m_size not yet bumped)
  }
  this->template columnByIndex<s_weightCol>()[m_size] = a_weight;
  this->template columnByIndex<s_idCol>()[m_size]     = s_invalidID;
  this->template columnByIndex<s_rankCol>()[m_size]   = static_cast<RankID>(-1);
  this->template columnByIndex<s_ghostCol>()[m_size]  = GhostType::Valid;
  this->writePayload(a_payload, std::make_index_sequence<s_numPayloadColumns>{});
  m_size++;
  m_sorted = false;
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::linearizeParticle(void* a_buffer, const std::size_t a_index) const noexcept
{
  CH_assert(a_index < m_size);

  auto* p = static_cast<unsigned char*>(a_buffer);
  this->linearizeImpl(p, a_index, std::make_index_sequence<s_numColumns>{});
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::delinearizeAndAppend(const void* a_buffer)
{
  this->growIfFull();
  const auto* p = static_cast<const unsigned char*>(a_buffer);
  this->delinearizeImpl(p, std::make_index_sequence<s_numColumns>{});
  m_size++;
  m_sorted = false;
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::h5LinearizeParticle(void* a_buffer, const std::size_t a_index) const noexcept
{
  CH_assert(a_index < m_size);

  auto* p = static_cast<unsigned char*>(a_buffer);
  for (int dir = 0; dir < SpaceDim; dir++) {
    detail::pushBytes(p, this->positionColumn(dir)[a_index]);
  }
  detail::pushBytes(p, this->template columnByIndex<s_weightCol>()[a_index]);
  this->h5LinearizePayload(p, a_index, H5PayloadColumns{});
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::h5DelinearizeAndAppend(const void* a_buffer)
{
  this->growIfFull();
  const auto* p = static_cast<const unsigned char*>(a_buffer);
  for (int dir = 0; dir < SpaceDim; dir++) {
    detail::pullBytes(p, this->positionColumn(dir)[m_size]);
  }
  detail::pullBytes(p, this->template columnByIndex<s_weightCol>()[m_size]);
  this->template columnByIndex<s_idCol>()[m_size]    = s_invalidID;
  this->template columnByIndex<s_rankCol>()[m_size]  = static_cast<RankID>(-1);
  this->template columnByIndex<s_ghostCol>()[m_size] = GhostType::Valid;

  // Non-checkpointed payload columns default-construct; then overwrite the checkpointed subset.
  const P def{};
  this->writePayload(def, std::make_index_sequence<s_numPayloadColumns>{});
  this->h5DelinearizePayload(p, H5PayloadColumns{});

  m_size++;
  m_sorted = false;
}

template <typename P, typename Traits>
inline void
ParticleSoA<P, Traits>::sortByCell(const Box& a_box, const RealVect& a_dx, const RealVect& a_probLo)
{
  // No CH_TIME here: this leaf container is the per-particle hot path; timing of bulk
  // operations like this one is done by the container/mesh drivers above.
  const std::size_t nCells  = a_box.numPts();
  const IntVect     loEnd   = a_box.smallEnd();
  const IntVect     boxSize = a_box.size();

  m_cellStart.assign(nCells + 1, 0);
  if (m_size == 0) {
    m_sorted = true;
    return;
  }

  // Pass 1: compute each particle's Fortran cell key and tally per-cell counts.
  std::vector<std::size_t> key(m_size);
  for (std::size_t i = 0; i < m_size; i++) {
    const RealVect x      = this->position(i);
    std::size_t    lin    = 0;
    std::size_t    stride = 1;
    for (int d = 0; d < SpaceDim; d++) {
      int c = static_cast<int>(std::floor((x[d] - a_probLo[d]) / a_dx[d])) - loEnd[d];
      if (c < 0) {
        c = 0;
      }
      if (c >= boxSize[d]) {
        c = boxSize[d] - 1;
      }
      lin += static_cast<std::size_t>(c) * stride;
      stride *= static_cast<std::size_t>(boxSize[d]);
    }
    key[i] = lin;
    m_cellStart[lin + 1]++;
  }

  // Prefix-sum the counts into CSR start offsets (m_cellStart[c] = start of cell c).
  for (std::size_t c = 1; c <= nCells; c++) {
    m_cellStart[c] += m_cellStart[c - 1];
  }

  // Pass 2: stable destination slot per particle, then physically reorder all columns.
  std::vector<std::size_t> cursor(m_cellStart.begin(), m_cellStart.begin() + nCells);
  std::vector<std::size_t> dest(m_size);
  for (std::size_t i = 0; i < m_size; i++) {
    dest[i] = cursor[key[i]]++;
  }
  this->permute(dest);

  m_sorted = true;
}

#include <CD_NamespaceFooter.H>

#endif
