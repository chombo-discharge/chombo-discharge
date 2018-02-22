/*!
  @file   tags.cpp
  @brief  Implementation of tags.H
  @author Robert Marskar
  @date   Feb. 2018
  @todo   Segregate implementation and header
*/

#include "tags.H"

tags::tags(){

}

tags::tags(const Box& a_box){
  m_ivs = DenseIntVectSet(a_box, false);
}

tags::~tags(){

}

DenseIntVectSet& tags::get_ivs() {
  return m_ivs;
}

const DenseIntVectSet& tags::get_ivs() const {
  return m_ivs;
}

void tags::copy(const Box&        a_region_from,
		const Interval&   a_dst_interval,
		const Box&        a_region_to,
		const tags&       a_src,
		const Interval&   a_src_interval){
  CH_assert(a_dst_interval == Interval(0,0));
  CH_assert(a_src_interval == Interval(0,0));
  CH_assert(a_region_from  == a_region_to);

  m_ivs = a_src.m_ivs;
}

int tags::preAllocatable(){
  return 2; // Dynamic object
}

int tags::size(const Box& a_box, const Interval& a_comps) const{
  return m_ivs.linearSize();
}

void tags::linearIn(const void* const a_in_buffer, const Box& R, const Interval& a_comps){
  m_ivs.linearIn(a_in_buffer);
}

void tags::linearOut(void* const a_out_buffer, const Box& R, const Interval& a_comps) const{
  m_ivs.linearOut(a_out_buffer);
}
