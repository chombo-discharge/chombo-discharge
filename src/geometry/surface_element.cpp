/*!
  @brief  surface_element.cpp
  @brief  Implementation of surface_element.H
  @author Robert Marskar
  @date   Apr. 2018
  @todo   Segregate implementation
*/

#include "surface_element.H"

surface_element::surface_element(){

}

surface_element::~surface_element(){

}

bv_sphere surface_element::get_bv_sphere(){
  return m_bv;
}
