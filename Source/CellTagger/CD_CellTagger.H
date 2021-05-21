/* chombo-discharge
 * Copyright 2021 SINTEF Energy Research
 * Please refer to LICENSE in the chombo-discharge root directory
 */

/*!
  @file   CD_CellTagger.H
  @brief  Class for tagging cells across an AMR hierarchy
  @author Robert Marskar
*/

#ifndef CD_CellTagger_H
#define CD_CellTagger_H

// Our includes
#include <CD_ComputationalGeometry.H>
#include <CD_RealBox.H>
#include <CD_AmrMesh.H>
#include <CD_EBAMRTags.H>
#include <CD_NamespaceHeader.H>

/*!
  @brief Abstract class for tagging cells across an AMR hierarchy
*/
class CellTagger {
public:

  /*!
    @brief Weak constructor, must subsequently call define
  */
  CellTagger();

  /*!
    @brief Destructor
  */
  virtual ~CellTagger();

  /*!
    @brief Regrid function for cell tagger (in case it uses transient storage to do things)
  */
  virtual void regrid() = 0;

  /*!
    @brief Parse class options
  */
  virtual void parseOptions() = 0;

  /*!
    @brief Parse class options
  */
  virtual void parseRuntimeOptions();

  /*!
    @brief Tag cells
    @param[out] a_tags Tags on grid levels
  */
  virtual bool tagCells(EBAMRTags& a_tags) = 0;

  /*!
    @brief Get number of tracer fields
    @return Returns number of plot variables that CellTagger will write to plot files. 
  */
  virtual int getNumberOfPlotVariables();

  /*!
    @brief Write plot data
    @param[inout] a_output Output data holder
    @param[inout] a_plotVariableNames Plot variable names
    @param[inout] a_icomp Starting variable in a_output where we begin appending data.
  */
  virtual void writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp);

  /*!
    @brief Get tag buffer
  */
  int getBuffer();

protected:

  /*!
    @brief Class name
  */
  std::string m_name;

  /*!
    @brief Verbosity
  */
  int m_verbosity;

  /*!
    @brief Buffer. Grow tags by this amount (finest level only)
  */
  int m_buffer;

  /*!
    @brief Boxes used for tagging
  */
  Vector<RealBox> m_tagBoxes;

  /*!
    @brief Check if a point is inside a tagging box
  */
  bool insideTagBox(const RealVect a_pos);

  /*!
    @brief Parse boxes
  */
  virtual void parseBoxes();

  /*!
    @brief Parse buffer
  */
  virtual void parseBuffer();

  /*!
    @brief Set verbosity
  */
  virtual void parseVerbosity();
};

#include <CD_NamespaceFooter.H>

#endif
