/*!
  @file   BoxType.cpp
  @brief  Implementation of BoxType.H
  @author Robert Marskar
*/

#include "BoxType.H"

int BoxType::preAllocatable(){
  return 2;
}

BoxType::BoxType(){
  m_which = -1;
}

BoxType::BoxType(const Box& box, int comps){
  m_which = -1;
}

void BoxType::setCovered(){
  m_which = 0;
}

void BoxType::setRegular(){
  m_which = 1;
}

void BoxType::setCutCell(){
  m_which = 2;
}

bool BoxType::isCovered() const{
  return m_which==0;
}

bool BoxType::isRegular() const{
  return m_which==1;
}

bool BoxType::isCutCell() const{
  return m_which==2;
}

void BoxType::define(const Box& box, int comps){
  m_which = -1;
}

void BoxType::copy(const Box& R, const Interval& Cd, const BoxType& source, const Interval Cs){
  m_which = source.m_which;
}

void BoxType::copy(const Box&      a_RegionFrom,
		   const Interval& a_Cdest,
		   const Box&      a_RegionTo,
		   const BoxType&  a_src,
		   const Interval& a_Csrc){

  if(a_RegionFrom.contains(a_RegionTo)){
    m_which = a_src.m_which;
  }
  else{
    m_which = -1;
  }
}

int BoxType::size(const Box& R, const Interval& comps) const{
  return sizeof(int);
}

void BoxType::linearOut(void* buf, const Box& R, const Interval& comps) const{
  int* buffer = (int*) buf;

  *buffer = m_which;
  buffer++;
}

void BoxType::linearIn(void* buf, const Box& R, const Interval& comps){
  int* buffer = (int*) buf;
  m_which = *buffer;
  buffer++;
}

void BoxType::operator=(const BoxType& a_src){
  m_which = a_src.m_which;
}

BoxTypeFactory::BoxTypeFactory(){

}

BoxTypeFactory::~BoxTypeFactory(){

}

BoxType* BoxTypeFactory::create(const Box& a_box, int a_ncomps, const DataIndex& a_dit) const{
  return new BoxType();
}
