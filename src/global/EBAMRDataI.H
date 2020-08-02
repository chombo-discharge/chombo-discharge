/*!
  @file   EBAMRDataI.H
  @brief  Implementation of EBAMRData.H
  @author Robert Marskar
  @data   Aug. 2020
*/

#ifndef _EBAMRDATAI_H_
#define _EBAMRDATAI_H_

#include "EBAMRData.H"

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
const std::string EBAMRData<T>::get_realm() const{
  return m_realm;
}

template <typename T>
void EBAMRData<T>::set_realm(const std::string a_realm){
  m_realm = a_realm;
}

template <typename T>
void EBAMRData<T>::copy(const EBAMRData<T>& a_other){

  for (int lvl = 0; lvl < m_data.size(); lvl++){
    RefCountedPtr<LevelData<T> >& src       = m_data[lvl];
    const RefCountedPtr<LevelData<T> >& dst = a_other[lvl];
    if(!src.isNull() && !dst.isNull()){
      if(m_realm == a_other.get_realm()){
	dst->localCopyTo(*src);
      }
      else{
	dst->copyTo(*src);
      }
    }
  }
}

// Explicit templates
template class EBAMRData<MFCellFAB>;
template class EBAMRData<MFFluxFAB>;
template class EBAMRData<MFBaseIVFAB>;
template class EBAMRData<EBCellFAB>;
template class EBAMRData<EBFluxFAB>;
template class EBAMRData<EBFaceFAB>;
template class EBAMRData<BaseIVFAB<Real> >;
template class EBAMRData<DomainFluxIFFAB>;
template class EBAMRData<BaseFab<bool> >;

// Typedefs for simple typing.
typedef EBAMRData<MFCellFAB>        MFAMRCellData;
typedef EBAMRData<MFFluxFAB>        MFAMRFluxData;
typedef EBAMRData<MFBaseIVFAB>      MFAMRIVData;
typedef EBAMRData<EBCellFAB>        EBAMRCellData;
typedef EBAMRData<EBFluxFAB>        EBAMRFluxData;
typedef EBAMRData<EBFaceFAB>        EBAMRFaceData;
typedef EBAMRData<BaseIVFAB<Real> > EBAMRIVData;
typedef EBAMRData<DomainFluxIFFAB>  EBAMRIFData;
typedef EBAMRData<BaseFab<bool> >   EBAMRBool;

#endif

