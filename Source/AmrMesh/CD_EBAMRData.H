/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_EBAMRData.H
  @brief  Class for holding data across EBAMR hierarchies
  @author Robert Marskar
*/

#ifndef CD_EBAMRData_H
#define CD_EBAMRData_H

// Chombo includes
#include <LevelData.H>
#include <RefCountedPtr.H>
#include <EBCellFAB.H>
#include <EBFluxFAB.H>
#include <EBFaceFAB.H>
#include <BaseIVFAB.H>
#include <BaseIFFAB.H>
#include <MFCellFAB.H>
#include <MFFluxFAB.H>
#include <BinItem.H>

// Our includes
#include <CD_DomainFluxIFFAB.H>
#include <CD_MFBaseIVFAB.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Default class for holding LevelData<T> data across an EBAMR realm
  @details This class is a base class for holding LevelData<T> over a realm. The class does not have arithmetic operators or 
  assignment operators because LevelData<T> does not have them. Manipulating the data is done directly by fetching the data 
  on a level through the indexing operator [], or by fetching all the data through getData. 
*/
template <typename T>
class EBAMRData
{
public:
  /*!
    @brief Default constructor. Does nothing.
  */
  EBAMRData();

  /*!
    @brief Destructor. Sets m_data to size a_size but does not allocate anything. 
  */
  EBAMRData(const int a_size);

  /*!
    @brief Destructor. Does nothing.
  */
  ~EBAMRData();

  /*!
    @brief Indexing operators. Returns m_data[a_level]
    @param[in] a_level Grid level
  */
  RefCountedPtr<LevelData<T>>&
  operator[](int a_level);

  /*!
    @brief Indexing operators. Returns const m_data[a_level]
    @param[in] a_level Grid level
  */
  const RefCountedPtr<LevelData<T>>&
  operator[](int a_level) const;

  /*!
    @brief Resize m_data to having size a_size
    @param[in] a_size Data size
  */
  void
  resize(const int a_size);

  /*!
    @brief Sets the realm for this object
    @param[in] a_realm Realm name
  */
  void
  setRealm(const std::string a_realm);

  /*!
    @brief Get size of m_data
  */
  int
  size() const;

  /*!
    @brief Exchange on all AMR levels
  */
  void
  exchange() noexcept;

  /*!
    @brief Get underlying data. Returns m_data. 
  */
  Vector<RefCountedPtr<LevelData<T>>>&
  getData();

  /*!
    @brief Get underlying data. Returns const m_data.
  */
  const Vector<RefCountedPtr<LevelData<T>>>&
  getData() const;

  /*!
    @brief General copying of another object. This does a local copy if the two realms are the same. Otherwise it does
    a communication copy. In Chombo-speak this is localCopyTo versus copyTo. 
    @param[in] a_other Other EBAMRData object
  */
  void
  copy(const EBAMRData<T>& a_other);

  /*!
    @brief General copying of another object. This does a local copy if the two realms are the same. Otherwise it does
    a communication copy. In Chombo-speak this is localCopyTo versus copyTo. 
    @param[in] a_srcInterval Source interval in a_other
    @param[in] a_other Other EBAMRData object
    @param[in] a_stInterval  Destination interval in *this
  */
  void
  copy(const Interval a_srcInterval, const EBAMRData<T>& a_other, const Interval a_dstInterval);

  /*!
    @brief Push a LevelData<T> object to the back of the data vector
    @param[in] a_levelData Data to append to vector
  */
  void
  push_back(RefCountedPtr<LevelData<T>>& a_levelData);

  /*!
    @brief Returns the string identifier for whatever realm this data is supposed to be allocated over. 
  */
  const std::string
  getRealm() const;

protected:
  /*!
    @brief Identifier for realm
  */
  std::string m_realm;

  /*!
    @brief Underlying data
  */
  Vector<RefCountedPtr<LevelData<T>>> m_data;
};

#include <CD_NamespaceFooter.H>

#include <CD_EBAMRDataImplem.H>

#endif
