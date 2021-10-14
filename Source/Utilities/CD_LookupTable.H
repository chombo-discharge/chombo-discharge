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

// Std includes
#include <iostream>
#include <vector>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for looking up and interpolation (x,y) data. 
*/
template<int N>
class LookupTable {
public:

  /*!
    @brief Default constructor. Creates a table without any entries. 
  */
  LookupTable();

  /*!
    @brief Copy constructor. 
  */
  LookupTable(const LookupTable& a_table);

  /*!
    @brief Destructor (does nothing). 
  */
  ~LookupTable();

  /*!
    @brief Clear all data
  */
  inline
  void clear();

  /*!
    @brief Make the table into a uniform table with the number of entries. 
    @details This will make the data uniform with constant spacing along the specified column data. 
    @note The spacing that is used is based on the smallest and largest value in the specified column. 
  */
  inline
  void makeUniform(const int a_numRows);

  /*!
    @brief Scale data in column by a specified factor. 
    @details The template parameter is the column to be scaled. 
    @param[in] a_scale       Factor to scale by
    @param[in] a_whichColumn Column to scale
  */
  template <int K>
  inline
  void scale(const Real a_scale);

  /*!
    @brief Swap two columns
    @param[in] a_column1 First column
    @param[in] a_column2 Second column
  */
  inline
  void swap(const int a_column1, const int a_column2);

  /*!
    @brief Restrict table range to maximum and minimum. Data outside these ranges are discarded. 
    @param[in] a_min         Minimum x-value in the table. 
    @param[in] a_max         Maximum x-value in the table. 
    @param[in] a_whichColumn Ranges based on this column. 
  */
  inline
  void setRange(const Real a_min, const Real a_max, const int a_whichColumn = 0);

  /*!
    @brief Sort the table. Useful if the input data was not already sorted. 
    @param[in] a_whichColumn Sort along this column. 
    @details On output, the data in the table is sorted -- smallest x goes first. 
  */
  inline
  void sort(const int a_whichColumn = 0);

  /*!
    @brief Add entry
  */
  template<typename... Ts>  
  inline
  void addEntry(const Ts&... x);

  /*!
    @brief Dump table to pout()
  */
  inline
  void dumpTable(std::ostream& a_outputStream = std::cout) const ;

  /*!
    @brief Get entry
    @note When calling this, make sure the data is sorted AND uniform. For performance reasons this is not checked when calling this routine. 
  */
  template <int K>
  inline
  Real getEntry(const Real a_x) const;

  /*!
    @brief Swap x and y vectors
  */
  inline
  void swapXY(){
    return this->swap(0,1);
  }

protected:

  /*!
    @brief For fetching/setting tuple stuff
  */
  enum StateStatus {
    Column  = 0,
    Sorted  = 1,
    Uniform = 2,
  };
  
  /*!
    @brief Table state. First entry is the column along which we sort and make the table uniform. The second entry is 
    whether or not the table is actually sorted, and the third entry is whether or not the table is uniform along the sorted column.
  */
  std::tuple<int, bool, bool> m_state;

  /*!
    @brief A base resolution for the electric field tabulations (we assumed that the rates were mapped out linearly)
  */
  Real m_delta;

  /*!
    @brief Data
  */
  std::vector<std::array<Real, N> > m_data;

  /*!
    @brief Direct lookup, assumes non-uniform ordering. 
    @details This is very slow, because it will look through the entire table. 
    @note The data MUST be sorted along the input column for this to work. 
  */
  inline
  std::array<Real, N> directLookup(const Real a_x) const;  
};

#include <CD_NamespaceFooter.H>

#include <CD_LookupTableImplem.H>

#endif
