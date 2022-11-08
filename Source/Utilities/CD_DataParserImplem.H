/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_DataParser.cpp
  @brief  Implementation of CD_DataParser.H
  @author Robert Marskar
*/

#ifndef CD_DataParserImplem_H
#define CD_DataParserImplem_H

// Std includes
#include <fstream>
#include <sstream>

// Chombo includes
#include <REAL.H>

// Our includes
#include <CD_DataParser.H>
#include <CD_NamespaceHeader.H>

namespace DataParser {
  LookupTable1D<2>
  simpleFileReadASCII(const std::string       a_fileName,
                      const int               a_xColumn,
                      const int               a_yColumn,
                      const std::vector<char> a_ignoreChars)
  {
    // This is the return table. It will be populated as we read the file.
    LookupTable1D<2> returnTable = LookupTable1D<2>();

    // Open an input file stream and start reading lines.
    std::ifstream inputFile(a_fileName);
    std::string   line;

    while (std::getline(inputFile, line)) {

      // Right trim line.
      line.erase(line.find_last_not_of(" \n\r\t") + 1);

      // Only do stuff for lines that are not empty. If we have an empty line we just proceed to the next one.
      if (!line.empty()) {

        // Check if we should parse the line. The user will have input a bunch of comment symbols which indicate lines that are comments.
        // If the line starts with any of those symbols we do not parse it.
        bool parseThisLine = true;
        for (const auto& ignoreChar : a_ignoreChars) {
          if (line.at(0) == ignoreChar) {
            parseThisLine = false;
          }
        }

        if (parseThisLine) {

          // Read the lines into a vector of values. We assume that each line consists of a number of columns. The
          // value in each column gets a place in the vector.
          std::istringstream iss(line);

          double              curVal;
          std::vector<double> values;
          while (iss >> curVal) {
            values.emplace_back(curVal);
          }

          // Read the data into the LookupTable1D IF we have enough data in the input file. Rows that do not have enough data WILL be ignored.
          const int numColumnsOnThisLine = values.size();
          if (a_xColumn < numColumnsOnThisLine && a_yColumn < numColumnsOnThisLine) {
            returnTable.addEntry(values[a_xColumn], values[a_yColumn]);
          }
        }
      }
    }

    return returnTable;
  }

  LookupTable1D<2>
  fractionalFileReadASCII(const std::string       a_fileName,
                          const std::string       a_startRead,
                          const std::string       a_stopRead,
                          const int               a_xColumn,
                          const int               a_yColumn,
                          const std::vector<char> a_ignoreChars)
  {
    // This is the return table. It will be populated as we read the file.
    LookupTable1D<2> returnTable = LookupTable1D<2>();

    // Open an input file stream and start reading lines.
    bool          parseLine = false;
    std::ifstream inputFile(a_fileName);
    std::string   line;

    while (std::getline(inputFile, line)) {

      // Right trim line.
      line.erase(line.find_last_not_of(" \n\r\t") + 1);

      // Check if we should parse the line. We start and
      // stop and we encounter the input strings.
      if (line == a_startRead) {
        parseLine = true;
      }
      else if (line == a_stopRead && parseLine) {
        break;
      }

      if (parseLine) {
        // Check if we should parse the line. The user will have input a bunch of comment symbols which
        // indicate lines that are comments. If the line starts with any of those symbols we do not parse it.
        bool lineIsCommented = false;
        for (const auto& ignoreChar : a_ignoreChars) {
          if (line.at(0) == ignoreChar) {
            lineIsCommented = true;
          }
        }

        if (!lineIsCommented) {
          // Read the lines into a vector of values. We assume that each line consists of a number of columns. The
          // value in each column gets a place in the vector.
          std::istringstream iss(line);

          double              curVal;
          std::vector<double> values;
          while (iss >> curVal) {
            values.emplace_back(curVal);
          }

          // Read the data into the LookupTable1D IF we have enough data in the input file. Rows that do not have enough
          // data WILL be ignored.
          const int numColumnsOnThisLine = values.size();
          if (a_xColumn < numColumnsOnThisLine && a_yColumn < numColumnsOnThisLine) {
            returnTable.addEntry(values[a_xColumn], values[a_yColumn]);
          }
        }
      }
    }

    return returnTable;
  }
} // namespace DataParser

#include <CD_NamespaceFooter.H>

#endif
