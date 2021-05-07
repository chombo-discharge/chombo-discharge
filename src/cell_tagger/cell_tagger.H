/*!
  @file   cell_tagger.H
  @brief  Class for tagging cells across an AMR hierarchy
  @author Robert Marskar
  @date   Nov. 2017
  @todo   Rename to cell_tagger. 
*/

#ifndef _CELL_TAGGER_
#define _CELL_TAGGER_

#include "computational_geometry.H"
#include "amr_mesh.H"
#include "real_box.H"

#include "CD_NamespaceHeader.H"
  
using EBAMRTags = Vector<RefCountedPtr<LayoutData<DenseIntVectSet> > >;

/*!
  @brief Abstract class for tagging cells across an AMR hierarchy
*/
class cell_tagger {
public:

  /*!
    @brief Weak constructor, must subsequently call define
  */
  cell_tagger();

  /*!
    @brief Destructor
  */
  virtual ~cell_tagger();

  /*!
    @brief Regrid storage for tracer field
  */
  virtual void regrid() = 0;

  /*!
    @brief Parse class options
  */
  virtual void parse_options() = 0;

  /*!
    @brief Parse class options
  */
  virtual void parse_runtime_options();

  /*!
    @brief Tag cells
  */
  virtual bool tag_cells(EBAMRTags& a_tags) = 0;

  /*!
    @brief Get number of tracer fields
  */
  virtual int get_num_plot_vars();

  /*!
    @brief Write plot data
  */
  virtual void write_plot_data(EBAMRCellData& a_output, Vector<std::string>& a_plotvar_names, int& a_icomp);

  /*!
    @brief Get tag buffer
  */
  int get_buffer();

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
  Vector<real_box> m_tagboxes;

  /*!
    @brief Check if a point is inside a tagging box
  */
  bool inside_tag_box(const RealVect a_pos);

  /*!
    @brief Parse boxes
  */
  virtual void parse_boxes();

  /*!
    @brief Parse buffer
  */
  virtual void parse_buffer();

  /*!
    @brief Set verbosity
  */
  virtual void parse_verbosity();
};
#include "CD_NamespaceFooter.H"

#endif
