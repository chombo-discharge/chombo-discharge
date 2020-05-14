/*!
  @file   ito_iterator.H
  @brief  Implementation of ito_iterator.H
  @author Robert Marskar
  @date   May 2020
*/

#include "ito_iterator.H"

ito_iterator::ito_iterator(){
  CH_TIME("ito_iterator::ito_iterator");
  m_defined = false;
}


ito_iterator::ito_iterator(ito_layout& a_layout, const species_iteration::which_species a_mode){
  CH_TIME("ito_iterator::ito_iterator");
    
  m_solvers  = a_layout.get_solvers();
  m_species  = a_layout.get_species();
  m_num      = m_solvers.size();
  m_mode     = a_mode;
  m_defined  = true;

  reset();
}

ito_iterator::~ito_iterator(){

}

int ito_iterator::num_solvers(){
  CH_assert(m_defined);
  return m_num;
}

int ito_iterator::get_solver() const {
  CH_assert(m_defined);
  return m_isolver;
}

void ito_iterator::reset(){
  CH_assert(m_defined);
  m_isolver = 0;
}

bool ito_iterator::ok(){
  return (m_isolver < m_num);
}

void ito_iterator::operator++(){
  CH_assert(m_defined);

  if(m_mode == species_iteration::all){
    m_isolver++;
  }
  else{
    m_isolver++;

    for (int isolver = m_isolver; this->ok(); ++isolver){
      const RefCountedPtr<ito_solver>& solver = m_solvers[isolver];
      const RefCountedPtr<ito_species>& species   = m_species[isolver];
      
      if(m_mode == species_iteration::charged){
	if(species->get_charge() != 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == species_iteration::negative){
	if(species->get_charge() < 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == species_iteration::positive){
	if(species->get_charge() > 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode == species_iteration::neutral){
	if(species->get_charge() == 0){
	  m_isolver = isolver;
	  break;
	}
      }
    }
  }
}

RefCountedPtr<ito_solver>& ito_iterator::operator() () {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_solvers[m_isolver];
}

RefCountedPtr<ito_species>& ito_iterator::get_species() {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_species[m_isolver];
}

