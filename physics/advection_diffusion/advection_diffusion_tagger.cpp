/*!
  @file   advection_diffusion_tagger.cpp
  @brief  Implementation of advection_diffusion_tagger.H
  @author Robert Marskar
  @date   Nov. 2017
  @todo   Rename to advection_diffusion_tagger. 
*/

#include "advection_diffusion_tagger.H"

#include <ParmParse.H> 

using namespace physics::advection_diffusion;

advection_diffusion_tagger::advection_diffusion_tagger(RefCountedPtr<cdr_solver>& a_solver){
  m_solver = a_solver;
  m_name = "advection_diffusion";
}

advection_diffusion_tagger::~advection_diffusion_tagger(){

}

void advection_diffusion_tagger::regrid(){

}

void advection_diffusion_tagger::parse_options(){
  ParmParse pp(m_name.c_str());
  pp.get("refine_curv", m_refi_curv);
  pp.get("refine_magn", m_refi_magn);
  
  parse_buffer(); // Derived from cell
}

bool advection_diffusion_tagger::tag_cells(EBAMRTags& a_tags){
  return true;
}
