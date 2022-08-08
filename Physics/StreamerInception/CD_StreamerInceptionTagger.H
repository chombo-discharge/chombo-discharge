/* chombo-discharge
 * Copyright © 2022 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_StreamerInceptionTagger.H
  @brief  Cell tagger for StreamerInception module. 
  @author Robert Marskar
*/

#ifndef CD_StreamerInceptionTagger_H
#define CD_StreamerInceptionTagger_H

// Our includes
#include <CD_CellTagger.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace StreamerInception {

    /*!
      @brief Class for tagging cells in the StreamerInception module. 
      @details This flags cells for refinement if |grad(|E|)|dx /|E| > threshold. 
    */
    class StreamerInceptionTagger : public CellTagger
    {
    public:

      // Disallowed constructors. 
      StreamerInceptionTagger() = delete;
      StreamerInceptionTagger(const StreamerInceptionTagger&) = delete;
      StreamerInceptionTagger(const StreamerInceptionTagger&&) = delete;            

      /*!
	@brief Full constructor.
	@param[in] a_amrMesh       Amr mesh. 
	@param[in] a_electricField Pointer to electric field. 
	@param[in] a_phase         Default phase. 
      */
      StreamerInceptionTagger(const RefCountedPtr<AmrMesh>& a_amrMesh,
			      const EBAMRCellData* const a_electricField,
			      const phase::which_phase a_phase = phase::gas);

      /*!
	@brief Destructor
      */
      virtual ~StreamerInceptionTagger();

      /*!
	@brief Parse options
      */
      virtual void
      parseOptions() override final;

      /*!
	@brief Parse options
      */
      virtual void
      regrid() override final;

      /*!
	@brief Tag cells for refinement. 
	@details This tags based on |grad(E)|*dx/E. 
	@param[inout] a_tags Tags on grid levels.
      */
      virtual bool
      tagCells(EBAMRTags& a_tags) override final;      

    protected:

      /*!
	@brief AmrMesh reference. 
      */
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Pointer to electric field.
      */
      const EBAMRCellData* m_electricField;

      /*!
	@brief Tracer field. Will be populated with |grad(E)|*dx/|E|
      */
      EBAMRCellData m_tracerField;

      /*!
	@brief Phase
      */
      phase::which_phase m_phase;

      /*!
	@brief Realm
      */
      std::string m_realm;

      /*!
	@brief Curvature refinement factor
      */
      Real m_refCurv;
    };
  }
}

#include <CD_NamespaceFooter.H>

#endif
