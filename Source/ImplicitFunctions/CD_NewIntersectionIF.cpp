/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_NewIntersectionIF.cpp
  @brief  Implementation of NewIntersectionIF.H
  @author Robert Marskar
*/

// Std includes
#include <limits>

// Our includes
#include <CD_NewIntersectionIF.H>
#include <CD_NamespaceHeader.H>

NewIntersectionIF::NewIntersectionIF(){
  m_numFuncs = 0;
  m_impFuncs.resize(0);
}

NewIntersectionIF::NewIntersectionIF(const Vector<BaseIF*>& a_impFuncs) {

  m_numFuncs = a_impFuncs.size();

  m_impFuncs.resize(m_numFuncs);

  // Make copies
  for (int i = 0; i < m_numFuncs; i++){
    if(a_impFuncs[i] == nullptr){
      m_impFuncs[i] = nullptr;
    }
    else{
      m_impFuncs[i] = a_impFuncs[i]->newImplicitFunction();
    }
  }
}

NewIntersectionIF::~NewIntersectionIF() {

  // Delete the implicit functions.
  for (int i = 0; i < m_numFuncs; i++){
    if(m_impFuncs[i] != nullptr){
      delete m_impFuncs[i];
    }
  }
}

Real NewIntersectionIF::value(const RealVect& a_point) const {
  
  // Returned value. 
  Real retval = -std::numeric_limits<Real>::max();

  // Find the maximum value and return it
  if (m_numFuncs > 0) {
      retval = m_impFuncs[0]->value(a_point);

      for (int ifunc = 1; ifunc < m_numFuncs; ifunc++){
	Real cur;
	
	cur = m_impFuncs[ifunc]->value(a_point);
	if (cur > retval) {
	  retval = cur;
	}
      }
  }

  return retval;
}

BaseIF* NewIntersectionIF::newImplicitFunction() const {
  return static_cast<BaseIF*> (new NewIntersectionIF(m_impFuncs));
}

#include <CD_NamespaceFooter.H>
