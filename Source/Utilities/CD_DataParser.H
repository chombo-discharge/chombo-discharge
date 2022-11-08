/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DataParser.H
  @brief  Declaration of various routines for simple parsing of input data. 
  @author Robert Marskar
*/

#ifndef CD_DataParser_H
#define CD_DataParser_H

// Our includes
#include <CD_LookupTable1D.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Class for simple data parsing from IO files.
*/
namespace DataParser {

  /*!
    @brief Simple file parser which reads a file and puts the data into two columns (a lookup table). 
    @details This will read ASCII row/column data into a LookupTable1D (which is a simple x-y data structure). The user can specify which columns
    to use as the x and y coordinates in the LookupTable1D. If the user ask for a column that does not exist in the input file, the data will not be read
    into the LookupTable1D. 
    @param[in] a_fileName    Input file name. Must be an ASCII file organized into rows and columns. 
    @param[in] a_xColumn     Which column to use as the x-column. 
    @param[in] a_yColumn     Which column to use as the y-column. 
    @param[in] a_ignoreChars Characters indicating comments in the file. Lines starting with these characters are ignored. 
    @note xColumn = 0 is the FIRST column, xColumn=1 is the SECOND column and so on. 
  */
  LookupTable1D<2>
  simpleFileReadASCII(const std::string       a_fileName,
                      const int               a_xColumn     = 0,
                      const int               a_yColumn     = 1,
                      const std::vector<char> a_ignoreChars = {'#', '/'});

  /*!
    @brief Simple ASCII file parser which reads a file and puts the data into two columns (a lookup table)
    @details This version can read a partial file. It will ignore all lines until the a_startRead string is
    encountered. After that, it will read files as usual, skipping lines that begin with comments. Empty lines
    are NOT skipped, since an empty line typically identifies the end of a field from e.g. BOLSIG outputs. Instead
    to stop reading when an empty line is encountered, one can use a_stopRead="". 
    @param[in] a_fileName    Input file name. Must be an ASCII file organized into rows and columns. 
    @param[in] a_startRead   Identifier where we start parsing lines into the table
    @param[in] a_stopRead    Identifier where we stop parsing lines into the table
    @param[in] a_xColumn     Which column to use as the x-column. 
    @param[in] a_yColumn     Which column to use as the y-column. 
    @param[in] a_ignoreChars Characters indicating comments in the file. Lines starting with these characters are ignored. 
  */
  LookupTable1D<2>
  fractionalFileReadASCII(const std::string       a_fileName,
                          const std::string       a_startRead,
                          const std::string       a_stopRead,
                          const int               a_xColumn     = 0,
                          const int               a_yColumn     = 1,
                          const std::vector<char> a_ignoreChars = {'#', '/'});

} // namespace DataParser

#include <CD_NamespaceFooter.H>

#include <CD_DataParserImplem.H>

#endif
