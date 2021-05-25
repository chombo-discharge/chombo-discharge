/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_LookupTable.H
  @brief  Declaration of a table for looking up coefficients etc. 
  @author Robert Marskar
*/

#ifndef CD_LookupTable_H
#define CD_LookupTable_H

// Chombo includes
#include <RealVect.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for looking up and interpolation (x,y) data. 
*/
class LookupTable {
public:

  /*!
    @brief Constructor
  */
  LookupTable();

  /*!
    @brief Copy constructor with scalability
  */
  LookupTable(const LookupTable& a_table);

  /*!
    @brief Destructor
  */
  ~LookupTable();

  /*!
    @brief Make the table into a uniform table with the number of entries
  */
  inline
  void makeUniform(const int a_entries);

  /*!
    @brief Scale x values
  */
  inline
  void scaleX(const Real& a_scale);

  /*!
    @brief Scale y values
  */
  inline
  void scaleY(const Real& a_scale);
  
  /*!
    @brief Add entry
  */
  inline
  void addEntry(const Real& a_x, const Real& a_y);

  /*!
    @brief Add other table. 
  */
  inline
  void addTable(const LookupTable& a_otherTable, const Real a_scale);

  /*!
    @brief Dump table to pout()
  */
  inline
  void dumpTable() const ;

  /*!
    @brief Get entry
  */
  inline
  Real getEntry(const Real a_x) const;

  /*!
    @brief Direct lookup, assumes non-uniform ordering. Performance killer. 
  */
  inline
  Real directLookup(const Real a_x) const;

  /*!
    @brief Swap x and y vectors
  */
  inline
  void swapXY();

  /*!
    @brief Add a table to this one
  */
  inline
  void operator+=(const LookupTable& a_table);

  /*!
    @brief Subtract a table from this one
  */
  inline
  void operator-=(const LookupTable& a_table);

  /*!
    @brief Assignment operator
  */
  inline
  void operator=(const LookupTable& a_table);

protected:

  /*!
    @brief A base resolution for the electric field tabulations (we assumed that the rates were mapped out linearly)
  */
  Real m_dx;

  /*!
    @brief Number of entries in reaction table
  */
  int m_numEntries;
  
  /*!
    @brief Towsend field
  */
  std::vector<Real> m_x;

  /*!
    @brief Rate coefficient
  */
  std::vector<Real> m_y;
};

#include <CD_NamespaceFooter.H>

#include <CD_LookupTableImplem.H>

#endif
