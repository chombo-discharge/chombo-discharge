/*!
  @file   EBAMRDataI.H
  @brief  Implementation of EBAMRData.H
  @author Robert Marskar
  @data   Aug. 2020
*/

#ifndef _EBAMRDATAI_H_
#define _EBAMRDATAI_H_

#include "EBAMRData.H"

#include "CD_NamespaceHeader.H"

template <typename T>
EBAMRData<T>::EBAMRData(){
}

template <typename T>
EBAMRData<T>::EBAMRData(const int a_size){
  m_data.resize(a_size);
}

template <typename T>
EBAMRData<T>::~EBAMRData(){
}

template <typename T>
RefCountedPtr<LevelData<T> >& EBAMRData<T>::operator[](int a_level){
  return m_data[a_level];
}

template <typename T>
const RefCountedPtr<LevelData<T> >& EBAMRData<T>::operator[](int a_level) const {
  return m_data[a_level];
}

template <typename T>
void EBAMRData<T>::resize(const int a_size){
  m_data.resize(a_size);
}

template <typename T>
int EBAMRData<T>::size() {
  return m_data.size();
}

template <typename T>
const int EBAMRData<T>::size() const {
  return m_data.size();
}

template <typename T>
Vector<RefCountedPtr<LevelData<T> > >& EBAMRData<T>::get_data(){
  return m_data;
}

template <typename T>
const Vector<RefCountedPtr<LevelData<T> > >& EBAMRData<T>::get_data() const {
  return m_data;
}

template <typename T>
void EBAMRData<T>::push_back(RefCountedPtr<LevelData<T> >& a_levelData){
  m_data.push_back(a_levelData);
}

template <typename T>
const std::string EBAMRData<T>::getRealm() const{
  return m_Realm;
}

template <typename T>
void EBAMRData<T>::setRealm(const std::string a_realm){
  m_Realm = a_realm;
}

template <typename T>
void EBAMRData<T>::copy(const EBAMRData<T>& a_other){

  for (int lvl = 0; lvl < m_data.size(); lvl++){
    RefCountedPtr<LevelData<T> >& src       = m_data[lvl];
    const RefCountedPtr<LevelData<T> >& dst = a_other[lvl];
    if(!src.isNull() && !dst.isNull()){
      if(m_Realm == a_other.getRealm()){
	dst->localCopyTo(*src);
      }
      else{
	dst->copyTo(*src);
      }
    }
  }
}

template <typename T>
void EBAMRData<T>::copy(const Interval a_srcInterval, const EBAMRData<T>& a_other, const Interval a_dstInterval){

  for (int lvl = 0; lvl < m_data.size(); lvl++){
    RefCountedPtr<LevelData<T> >& src       = m_data[lvl];
    const RefCountedPtr<LevelData<T> >& dst = a_other[lvl];
    
    if(!src.isNull() && !dst.isNull()){
      if(m_Realm == a_other.getRealm()){
	dst->localCopyTo(a_srcInterval, *src, a_dstInterval);
      }
      else{
	dst->copyTo(a_srcInterval, *src, a_dstInterval);
      }
    }
  }
}

// Explicit templates
template class EBAMRData<MFCellFAB>;         
template class EBAMRData<MFFluxFAB>;         
template class EBAMRData<MFBaseIVFAB>;       
template class EBAMRData<EBCellFAB>;         // Cell-centered single-phase data
template class EBAMRData<EBFluxFAB>;         // Face-centered single-phase data
template class EBAMRData<EBFaceFAB>;           
template class EBAMRData<BaseIVFAB<Real> >;
template class EBAMRData<DomainFluxIFFAB>;
template class EBAMRData<BaseFab<bool> >;
template class EBAMRData<FArrayBox>;

// Typedefs for simple typing.
typedef EBAMRData<MFCellFAB>        MFAMRCellData;  // Cell-centered multifluid data
typedef EBAMRData<MFFluxFAB>        MFAMRFluxData;  // Face-centered multifluid data
typedef EBAMRData<MFBaseIVFAB>      MFAMRIVData;    // Irregular face multifluid data
typedef EBAMRData<EBCellFAB>        EBAMRCellData;  // Cell-centered single-phase data
typedef EBAMRData<EBFluxFAB>        EBAMRFluxData;  // Face-centered data in all coordinate direction
typedef EBAMRData<EBFaceFAB>        EBAMRFaceData;  // Face-centered in a single coordinate direction
typedef EBAMRData<BaseIVFAB<Real> > EBAMRIVData;    // Data on irregular data centroids
typedef EBAMRData<DomainFluxIFFAB>  EBAMRIFData;    // Data on domain phases
typedef EBAMRData<BaseFab<bool> >   EBAMRBool;      // For holding bool at every cell
typedef EBAMRData<FArrayBox>        EBAMRFAB;       // Just regular grid data.
#include "CD_NamespaceFooter.H"

#endif

