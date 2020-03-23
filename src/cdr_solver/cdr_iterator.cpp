/*!
  @file   cdr_iterator.H
  @brief  Implementation of cdr_iterator.H
  @author Robert Marskar
  @date   June Dec. 2017
*/

#include "cdr_iterator.H"

cdr_iterator::cdr_iterator(){
  CH_TIME("cdr_iterator::cdr_iterator");
  m_defined = false;
}


cdr_iterator::cdr_iterator(cdr_layout& a_layout, const species_iteration::which_species a_mode){
  CH_TIME("cdr_iterator::cdr_iterator");
    
  m_solvers  = a_layout.get_solvers();
  m_species  = a_layout.get_species();
  m_num      = m_solvers.size();
  m_mode     = a_mode;
  m_defined  = true;

  reset();
}

cdr_iterator::cdr_iterator(cdr_layout& a_layout, Vector<EBAMRCellData*> a_data){
  CH_TIME("cdr_iterator::cdr_iterator");
    
  m_solvers  = a_layout.get_solvers();
  m_species  = a_layout.get_species();
  m_num      = m_solvers.size();
  m_defined  = true;
  m_mode     = species_iteration::all;

  m_celldata = a_data;

  reset();
}

cdr_iterator::~cdr_iterator(){

}

int cdr_iterator::num_solvers(){
  CH_assert(m_defined);
  return m_num;
}

int cdr_iterator::get_solver() const {
  CH_assert(m_defined);
  return m_isolver;
}

void cdr_iterator::reset(){
  CH_assert(m_defined);
  m_isolver = 0;
}

bool cdr_iterator::ok(){
  return (m_isolver < m_num);
}

void cdr_iterator::operator++(){
  CH_assert(m_defined);

  if(m_mode == species_iteration::all){
    m_isolver++;
  }
  else{
    m_isolver++;

    for (int isolver = m_isolver; this->ok(); ++isolver){
      const RefCountedPtr<cdr_solver>& solver = m_solvers[isolver];
      const RefCountedPtr<cdr_species>& species   = m_species[isolver];
      
      if(m_mode = species_iteration::charged){
	if(species->get_charge() != 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode = species_iteration::negative){
	if(species->get_charge() < 0){
	  m_isolver = isolver;
	  break;
	}
      }
      else if(m_mode = species_iteration::positive){
	if(species->get_charge() > 0){
	  m_isolver = isolver;
	  break;
	}
      }
    }
  }
}

RefCountedPtr<cdr_solver>& cdr_iterator::operator() () {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_solvers[m_isolver];
}

RefCountedPtr<cdr_species>& cdr_iterator::get_species() {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return m_species[m_isolver];
}

EBAMRCellData& cdr_iterator::get_data() {
  CH_assert(m_defined == true);
  CH_assert(m_isolver < m_num);
  return *m_celldata[m_isolver];
}
