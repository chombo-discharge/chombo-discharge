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
Vector<RefCountedPtr<LevelData<T> > >& EBAMRData<T>::data(){
  return m_data;
}

template <typename T>
const Vector<RefCountedPtr<LevelData<T> > >& EBAMRData<T>::data() const {
  return m_data;
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

#endif

