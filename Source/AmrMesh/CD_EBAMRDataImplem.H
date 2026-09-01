/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
 * @file   CD_EBAMRDataImplem.H
 * @brief  Implementation of EBAMRData.H
 * @author Robert Marskar
 */

#ifndef CD_EBAMRDATAIMPLEM_H
#define CD_EBAMRDATAIMPLEM_H

// Chombo includes
#include <CH_Timer.H>

// Our includes
#include <CD_EBAMRData.H>
#include <CD_NamespaceHeader.H>

template <typename T>
EBAMRData<T>::EBAMRData() noexcept
{
  CH_TIME("EBAMRData<T>::EBAMRData()");
  m_reserved = false;
}

template <typename T>
EBAMRData<T>::EBAMRData(const int a_size) noexcept : EBAMRData<T>()
{
  CH_TIME("EBAMRData<T>::EBAMRData(int)");

  m_data.resize(a_size);
}

template <typename T>
EBAMRData<T>::~EBAMRData() noexcept
{
  CH_TIME("EBAMRData<T>::~EBAMRData()");
}

template <typename T>
RefCountedPtr<LevelData<T>>&
EBAMRData<T>::operator[](int a_level) noexcept
{
  return m_data[a_level];
}

template <typename T>
const RefCountedPtr<LevelData<T>>&
EBAMRData<T>::operator[](int a_level) const noexcept
{
  return m_data[a_level];
}

template <typename T>
void
EBAMRData<T>::resize(const int a_size) noexcept
{
  m_data.resize(a_size);
}

template <typename T>
int
EBAMRData<T>::size() const noexcept
{
  return m_data.size();
}

template <typename T>
void
EBAMRData<T>::clear() noexcept
{
  for (int lvl = 0; lvl < m_data.size(); lvl++) {
    m_data[lvl] = RefCountedPtr<LevelData<T>>(0);
  }

  m_data.resize(0);
}

template <typename T>
void
EBAMRData<T>::reserve() noexcept
{
  m_reserved = true;
}

template <typename T>
void
EBAMRData<T>::release() noexcept
{
  m_reserved = false;
}

template <typename T>
bool
EBAMRData<T>::isReserved() const noexcept
{
  return m_reserved;
}

template <typename T>
void
EBAMRData<T>::exchange() noexcept
{
  CH_TIMERS("EBAMRData<T>::exchange()");
  CH_TIMER("EBAMRData<T>::exchange()::post", t1);
  CH_TIMER("EBAMRData<T>::exchange()::local", t2);
  CH_TIMER("EBAMRData<T>::exchange()::wait", t3);

  // Pipeline the hierarchy instead of running one complete MPI round per level. Every level's messages go onto
  // the wire first, then the on-rank copies run while they are in flight, and only then does anything wait. The
  // exposed cost is the slowest single level rather than the sum over levels, and no level's messages are held
  // back by the on-rank work of the levels below it.
  //
  // This is safe because the send and receive buffers, and the MPI requests, live on the Copier -- and each
  // level's LevelData owns its own exchange Copier. The levels therefore never alias one another's buffers.
  // Message matching is unchanged as well: every rank still posts its levels in the same order, so MPI's
  // non-overtaking guarantee pairs sends with receives exactly as it does when the levels run one at a time.
  //
  // The IntVect::Zero test is the one LevelData::exchange() makes for itself. A level with no ghost cells has
  // nothing to exchange, so it must be skipped here rather than handed to exchangePost().
  const int numLevels = m_data.size();

  CH_START(t1);
  for (int lvl = 0; lvl < numLevels; lvl++) {
    if (m_data[lvl]->ghostVect() != IntVect::Zero) {
      m_data[lvl]->exchangePost(m_data[lvl]->exchangeCopier());
    }
  }
  CH_STOP(t1);

  CH_START(t2);
  for (int lvl = 0; lvl < numLevels; lvl++) {
    if (m_data[lvl]->ghostVect() != IntVect::Zero) {
      m_data[lvl]->exchangeLocal(m_data[lvl]->exchangeCopier());
    }
  }
  CH_STOP(t2);

  CH_START(t3);
  for (int lvl = 0; lvl < numLevels; lvl++) {
    if (m_data[lvl]->ghostVect() != IntVect::Zero) {
      m_data[lvl]->exchangeEnd();
    }
  }
  CH_STOP(t3);
}

template <typename T>
Vector<RefCountedPtr<LevelData<T>>>&
EBAMRData<T>::getData() noexcept
{
  return m_data;
}

template <typename T>
const Vector<RefCountedPtr<LevelData<T>>>&
EBAMRData<T>::getData() const noexcept
{
  return m_data;
}

template <typename T>
void
EBAMRData<T>::push_back(RefCountedPtr<LevelData<T>>& a_levelData) noexcept
{
  m_data.push_back(a_levelData);
}

template <typename T>
const std::string
EBAMRData<T>::getRealm() const noexcept
{
  return m_realm;
}

template <typename T>
void
EBAMRData<T>::setRealm(const std::string& a_realm) noexcept
{
  m_realm = a_realm;
}

// Explicit templates
template class EBAMRData<MFCellFAB>;
template class EBAMRData<MFFluxFAB>;
template class EBAMRData<MFBaseIVFAB>;
template class EBAMRData<EBCellFAB>; // Cell-centered single-phase data
template class EBAMRData<EBFluxFAB>; // Face-centered single-phase data
template class EBAMRData<EBFaceFAB>;
template class EBAMRData<BaseIVFAB<Real>>;
template class EBAMRData<DomainFluxIFFAB>;
template class EBAMRData<BaseFab<bool>>;
template class EBAMRData<FArrayBox>;

/// @cond DOXYGEN_SKIP
typedef EBAMRData<MFCellFAB>       MFAMRCellData; // Cell-centered multifluid data
typedef EBAMRData<MFFluxFAB>       MFAMRFluxData; // Face-centered multifluid data
typedef EBAMRData<MFBaseIVFAB>     MFAMRIVData;   // Irregular face multifluid data
typedef EBAMRData<EBCellFAB>       EBAMRCellData; // Cell-centered single-phase data
typedef EBAMRData<EBFluxFAB>       EBAMRFluxData; // Face-centered data in all coordinate direction
typedef EBAMRData<EBFaceFAB>       EBAMRFaceData; // Face-centered in a single coordinate direction
typedef EBAMRData<BaseIVFAB<Real>> EBAMRIVData;   // Data on irregular data centroids
typedef EBAMRData<DomainFluxIFFAB> EBAMRIFData;   // Data on domain phases
typedef EBAMRData<BaseFab<bool>>   EBAMRBool;     // For holding bool at every cell
typedef EBAMRData<FArrayBox>       EBAMRFAB;      // Just regular grid data.
/// @endcond

#include <CD_NamespaceFooter.H>

#endif
