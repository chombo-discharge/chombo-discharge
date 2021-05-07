/*!
  @file   cdr_iterator.H
  @brief  Iterator class for cdr_layout
  @author Robert Marskar
  @date   June Dec. 2017
*/

#ifndef _CDR_ITERATOR_
#define _CDR_ITERATOR_

#include "cdr_layout.H"

namespace ChomboDischarge {

  /*!
    @brief Iterator class for cdr_layout
  */
  template <class T>
  class cdr_iterator {
  public:

    /*!
      @brief Silly enum for iterating through selected species
    */
    enum class species_iteration{
      all,
      charged,
      negative,
      positive,
      neutral
    };
  
    /*!
      @brief Empty constructor
    */
    cdr_iterator();

    /*!
      @brief Full constructor
    */
    cdr_iterator(cdr_layout<T>& a_layout, const species_iteration a_mode = species_iteration::all);

    /*!
      @brief Full constructor, allows iteration through data
    */
    cdr_iterator(cdr_layout<T>& a_layout, Vector<EBAMRCellData*> a_data);

    /*!
      @brief Destructor
    */
    virtual ~cdr_iterator();

    /*!
      @brief Get number of solvers. 
    */
    virtual int num_solvers();

    /*!
      @brief Get current solver number
    */
    virtual int index() const;

    /*!
      @brief Reset iterator
    */
    virtual void reset();

    /*!
      @brief Ok or not
    */
    virtual bool ok();
  
    /*!
      @brief Increment. This is different for different iteration modes
    */
    virtual void operator++();

    /*!
      @brief Get solver
    */
    virtual RefCountedPtr<T>& operator() ();

    /*!
      @brief Get current ion
    */
    virtual RefCountedPtr<cdr_species>& get_species();

    /*!
      @brief Get current data
    */
    virtual EBAMRCellData& get_data();
  
  protected:

    /*!
      @brief Which species to iterate over
    */
    species_iteration m_mode;
  
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
    Vector<RefCountedPtr<T> > m_solvers;

    /*!
      @brief Ions
    */
    Vector<RefCountedPtr<cdr_species> > m_species;

    /*!
      @brief Cell-centered data
    */
    Vector<EBAMRCellData*> m_celldata;

  };
}

#include "cdr_iteratorI.H"

#endif
