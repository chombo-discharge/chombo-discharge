/*!
  @file cdr_iterator.hpp
  @brief Iterator class for cdr_layout
  @author Robert Marskar
  @date June 2016
*/

#ifndef _CDR_ITERATOR_
#define _CDR_ITERATOR_

#include "cdr_layout.H"

/*!
  @brief Iterator class for cdr_layout
*/
class cdr_iterator {
public:  
  /*!
    @brief Empty constructor
  */
  cdr_iterator(){
    CH_TIME("cdr_iterator::cdr_iterator");
    m_defined = false;
  }

  /*!
    @brief Full constructor
  */
  cdr_iterator(cdr_layout& a_layout){
    CH_TIME("cdr_iterator::cdr_iterator");
    
    m_solvers  = a_layout.get_solvers();
    m_species  = a_layout.get_species();
    m_num      = m_solvers.size();
    m_defined  = true;

    reset();
  }

  /*!
    @brief Get number of solvers
  */
  virtual int num_solvers(){
    CH_assert(m_defined);
    return m_num;
  }

  /*!
    @brief Get current solver number
  */
  virtual int get_solver() const {
    CH_assert(m_defined);
    return m_isolver;
  }

  /*!
    @brief Reset iterator
  */
  virtual void reset(){
    CH_assert(m_defined);
    m_isolver = 0;
  }

  /*!
    @brief Ok or not
   */
  virtual bool ok(){
    return (m_isolver < m_num);
  }

  /*!
    @brief Increment
  */
  virtual void operator++(){
    CH_assert(m_defined);
    m_isolver++;
  }

  /*!
    @brief Get solver
  */
  virtual RefCountedPtr<cdr_solver>& operator() () {
    CH_assert(m_defined == true);
    CH_assert(m_isolver < m_num);
    return m_solvers[m_isolver];
  }

  /*!
    @brief Get current ion
  */
  virtual RefCountedPtr<species>& get_species() {
    CH_assert(m_defined == true);
    CH_assert(m_isolver < m_num);
    return m_species[m_isolver];
  }
  
protected:
  
  /*!
    @brief 
  */
  int m_isolver;

  /*!
    @brief Number of solvers
   */
  int m_num;;

  /*!
    @brief Defined or not
  */
  bool m_defined;

  /*!
    @brief Ion solvers
  */
  Vector<RefCountedPtr<cdr_solver> > m_solvers;

  /*!
    @brief Ions
  */
  Vector<RefCountedPtr<species> > m_species;

};
#endif
