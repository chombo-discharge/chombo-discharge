/*!
  @file   dcel_edge_iterator.cpp
  @brief  Implementation of dcel_edge_iterator.H
  @author Robert Marskar
  @date   Apr. 2018
*/

#include "dcel_iterator.H"
#include "dcel_poly.H"
#include "dcel_edge.H"
#include "dcel_vert.H"

edge_iterator::edge_iterator(){

}

edge_iterator::edge_iterator(dcel_poly& a_poly){
  m_polymode  = true;
  m_begin     = a_poly.get_edge();
  m_current   = m_begin;
  m_full_loop = false;
}

edge_iterator::edge_iterator(dcel_vert& a_vert){
  m_polymode  = false;
  m_begin     = a_vert.get_edge();
  m_current   = m_begin;
  m_full_loop = false;
}

RefCountedPtr<dcel_edge>& edge_iterator::operator() (){
  return m_current;
}

void edge_iterator::reset(){
  m_current = m_begin;
}

void edge_iterator::operator++(){
  if(m_polymode){
    m_current = m_current->get_next();
  }
  else{
    m_current = m_current->get_prev()->get_pair();
    //    m_current = m_current->get_pair()->get_next();
  }

  // Signal a full loop around the polygon. 
  if(m_current == m_begin){ 
    m_full_loop = true;
  }
}

bool edge_iterator::ok(){
  if(m_current != NULL && !m_full_loop){
    return true;
  }
  else {
    return false;
  }
}
