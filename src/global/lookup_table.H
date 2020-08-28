/*!
  @file lookup_table.H
  @brief Declaration of a table for looking up coefficients etc. 
  @author Robert Marskar
  @date Sept. 2016
*/

#ifndef _LOOKUPTABLE_H_
#define _LOOKUPTABLE_H_

#include <RealVect.H>

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
  lookup_table(const lookup_table& a_table);

  /*!
    @brief Destructor
  */
  ~lookup_table();

  /*!
    @brief Make the table into a uniform table with the number of entries
  */
  inline
  void make_uniform(const int a_entries);

  /*!
    @brief Scale x values
  */
  inline
  void scale_x(const Real& a_scale);

  /*!
    @brief Scale y values
  */
  inline
  void scale_y(const Real& a_scale);
  
  /*!
    @brief Add entry
  */
  inline
  void add_entry(const Real& a_x, const Real& a_y);

  /*!
    @brief Add other table. 
  */
  inline
  void add_table(const lookup_table& a_otherTable, const Real a_scale);

  /*!
    @brief Dump table to pout()
  */
  inline
  void dump_table() const ;

  /*!
    @brief Get entry
  */
  inline
  Real get_entry(const Real a_x) const;

  /*!
    @brief Direct lookup, assumes non-uniform ordering. Performance killer. 
  */
  inline
  Real direct_lookup(const Real a_x) const;

  /*!
    @brief Swap x and y vectors
  */
  inline
  void swap_xy();

  /*!
    @brief Add a table to this one
  */
  inline
  void operator+=(const lookup_table& a_table);

  /*!
    @brief Subtract a table from this one
  */
  inline
  void operator-=(const lookup_table& a_table);

  /*!
    @brief Assignment operator
  */
  inline
  void operator=(const lookup_table& a_table);

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
    @brief Towsend field
  */
  std::vector<Real> m_x;

  /*!
    @brief Rate coefficient
  */
  std::vector<Real> m_y;
};

#include "lookup_tableI.H"

#endif
