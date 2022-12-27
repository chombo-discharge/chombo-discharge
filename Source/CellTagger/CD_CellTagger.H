/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_CellTagger.H
  @brief  Declaration of a base class for tagging cells across an AMR hierarchy. 
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
  @brief Base class for tagging cells across an AMR hierarchy.
  @details This is a pure class that users must implemented in order to enable refinement/coarsening of grid cells. This class will require implementing at least 
  three functions, one for parsing options, one for regridding, and one which flags/removes cell for refinement. This class can also provide data to plot files, in which
  case the user must overwrite getNumberOfPlotVariables() and writePlotData().
  A supporting utility for CellTagger is that the parent class provides functionality for restricting cell refinement to specified regions in space. This can be used
  to prune tags where it is (heuristically) known that cell refinement is unnecessary. 
*/
class CellTagger
{
public:
  /*!
    @brief Weak constructor (does nothing)
  */
  CellTagger();

  /*!
    @brief Destructor (does nothing). 
  */
  virtual ~CellTagger();

  /*!
    @brief Regrid function for cell tagger (in case it uses transient storage to do things)
    @details This function exists because implementations may require data to be allocated on a mesh. This function is
    called by Driver to make sure data is reallocated when it needs to. 
  */
  virtual void
  regrid() = 0;

  /*!
    @brief Parse class options.
    @note This function is called by Driver. 
  */
  virtual void
  parseOptions() = 0;

  /*!
    @brief Parse runtime class options
    @details This function is called by Driver at every time step; the user can use this for run-time updates of e.g. refinement criteria.
  */
  virtual void
  parseRuntimeOptions();

  /*!
    @brief Tag cells
    @param[inout] a_tags Tags on grid levels
    @details EBAMRTags is a data-type Vector<RefCountedPtr<LayoutData<DenseIntVectSet> > >. The vector indicates the grid level, the LayoutData indicates data ownership (in
    much the same way as LevelData. The DenseIntVectSet is essentially an IntVectSet restricted to the patch (i.e. one cannot add IntVects that are outside the patch). 
    This function cells for refinement or coarsening. The user's responsibility is to add (or remove) tags from a_tags. 
  */
  virtual bool
  tagCells(EBAMRTags& a_tags) = 0;

  /*!
    @brief Get number of plot variables that will be written to file (by Driver).
    @return Returns number of plot variables that Driver will write to plot files. 
  */
  virtual int
  getNumberOfPlotVariables() const;

  /*!
    @brief Write plot data.
    @param[inout] a_output            Output data holder
    @param[inout] a_plotVariableNames Plot variable names
    @param[inout] a_icomp             Starting variable in a_output where we begin appending data.
  */
  virtual void
  writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const;

  /*!
    @brief Get tag buffer. 
    @details This returns an integer > 0 which adds additional tags around the cells that were refined. This is used by Driver to grow the refined region. 
    @return Returns m_buffer
    @note This will look for a string in the form "m_name.buffer = <integer>". 
  */
  virtual int
  getBuffer() const;

protected:
  /*!
    @brief Cell tagger class name (default is CellTagger).
  */
  std::string m_name;

  /*!
    @brief Class verbosity. 
  */
  int m_verbosity;

  /*!
    @brief Tag buffer. Grow tags by this amount (finest level only)
  */
  int m_buffer;

  /*!
    @brief Boxes that can be used to restrict tagging. 
  */
  Vector<RealBox> m_tagBoxes;

  /*!
    @brief Refinement boxes for adding manual refinement.
  */
  std::vector<std::pair<RealBox, int>> m_refBoxes;

  /*!
    @brief Check if a point is inside a tagging box
    @param[in] a_pos Physical coordinates. 
    @return Returns true if a_pos is inside any of the boxes in m_tagBoxes. Otherwise it returns false. 
  */
  bool
  insideTagBox(const RealVect a_pos) const;

  /*!
    @brief Get the specified level for this position when doing manual refinement. 
  */
  int
  getManualRefinementLevel(const RealVect a_pos) const;

  /*!
    @brief Parse boxes.
    @details This will look for input parameters in the form m_name.num_boxes = <integer>. If <integer> = 0 then the function will process additional lines
    in the form "m_name.box1_lo = <float> <float> <float>" which indicates the lower corner of the "tag box". 
  */
  virtual void
  parseTagBoxes();

  /*!
    @brief Parse refinement boxes
  */
  virtual void
  parseRefinementBoxes();

  /*!
    @brief Parse refinement buffer.
    @details This will look for an input parameter in the form "m_name.buffer = <integer>". If this is called but the string does not exist, it will cause
    a run-time error. 
  */
  virtual void
  parseBuffer();

  /*!
    @brief Parse the verbosity from the input file or command line.
    @details This will look for an input parameter in the form "m_name.verbosity = <integer>". This is used to control "verbosity" for the cell tagger class. Turning
    on verbosity for this class is mostly for debugging purposes. 
  */
  virtual void
  parseVerbosity();
};

#include <CD_NamespaceFooter.H>

#endif
