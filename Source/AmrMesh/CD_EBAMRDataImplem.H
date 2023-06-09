/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRDataImplem.H
  @brief  Implementation of EBAMRData.H
  @author Robert Marskar
*/

#ifndef CD_EBAMRDataImplem_H
#define CD_EBAMRDataImplem_H

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
  for (int lvl = 0; lvl < m_data.size(); lvl++) {
    m_data[lvl]->exchange();
  }
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
EBAMRData<T>::setRealm(const std::string a_realm) noexcept
{
  m_realm = a_realm;
}

template <typename T>
void
EBAMRData<T>::copy(const EBAMRData<T>& a_other) noexcept
{
  MayDay::Abort("EBAMRData::copy -- should have Copier in here for efficiency!");
  for (int lvl = 0; lvl < m_data.size(); lvl++) {
    RefCountedPtr<LevelData<T>>&       src = m_data[lvl];
    const RefCountedPtr<LevelData<T>>& dst = a_other[lvl];
    if (!src.isNull() && !dst.isNull()) {
      if (m_realm == a_other.getRealm()) {
        dst->localCopyTo(*src);
      }
      else {
        dst->copyTo(*src);
      }
    }
  }
}

template <typename T>
void
EBAMRData<T>::copy(const Interval a_srcInterval, const EBAMRData<T>& a_other, const Interval a_dstInterval) noexcept
{
  MayDay::Abort("EBAMRData::copy -- should have Copier in here for efficiency!");
  for (int lvl = 0; lvl < m_data.size(); lvl++) {
    RefCountedPtr<LevelData<T>>&       src = m_data[lvl];
    const RefCountedPtr<LevelData<T>>& dst = a_other[lvl];

    if (!src.isNull() && !dst.isNull()) {
      if (m_realm == a_other.getRealm()) {
        dst->localCopyTo(a_srcInterval, *src, a_dstInterval);
      }
      else {
        dst->copyTo(a_srcInterval, *src, a_dstInterval);
      }
    }
  }
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

// Typedefs for simple typing.
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

#include <CD_NamespaceFooter.H>

#endif
