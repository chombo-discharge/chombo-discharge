/*!
  @file lookup_table.H
  @brief Declaration of a table for looking up coefficients etc. 
  @author Robert Marskar
  @date Sept. 2016
*/

#ifndef _LOOKUPTABLE_
#define _LOOKUPTABLE_

#include <Vector.H>
#include <RealVect.H>
#include <ParmParse.H>


/*!
  @brief Class for looking up and interpolation (x,y) data. 
*/
class lookup_table {
public:

  /*!
    @brief Constructor
  */
  lookup_table();

  /*!
    @brief Copy constructor with scalability
   */
  lookup_table(const lookup_table* const a_table, const Real a_scale = 1.);

  /*!
    @brief Destructor
  */
  ~lookup_table();

  /*!
    @brief Scale table entries
  */
  virtual void scale(const Real& a_scale);
  
  /*!
    @brief Add entry
  */
  virtual void add_entry(const Real& a_x, const Real& a_y);

  /*!
    @brief Dump table to pout()
  */
  virtual void dump_table();

  /*!
    @brief Get entry
  */
  virtual Real get_entry(const Real a_x);

  /*!
    @brief Get entry
  */
  virtual const Real get_entry(const Real a_x) const;

protected:

  /*!
    @brief A base resolution for the electric field tabulations (we assumed that the rates were mapped out linearly)
   */
  Real m_dx;

  /*!
    @brief Number of entries in reaction table
  */
  int m_num_entries;

  /*!
    @brief Look up order, will be set automatically depending on the number of entries you put in
  */
  int m_order;
  
  /*!
    @brief Towsend field
  */
  Vector<Real> m_x;

  /*!
    @brief Rate coefficient
  */
  Vector<Real> m_y;
};

#endif
