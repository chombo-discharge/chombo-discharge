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
      @details This flags cells for refinement if alpha_eff * dx > threshold. User must supply the
      effective ionization coefficient through the constructor. 
    */
    class StreamerInceptionTagger : public CellTagger
    {
    public:
      /*!
	@brief Disallowed constructor
      */
      StreamerInceptionTagger() = delete;

      /*!
	@brief Disallowed constructor
      */
      StreamerInceptionTagger(const StreamerInceptionTagger&) = delete;

      /*!
	@brief Disallowed constructor
      */
      StreamerInceptionTagger(const StreamerInceptionTagger&&) = delete;

      /*!
	@brief Full constructor.
	@param[in] a_amrMesh       Amr mesh. 
	@param[in] a_electricField Pointer to electric field data.
	@param[in] a_alphaEff      Ionization coefficient
	@param[in] a_phase         Default phase. 
      */
      StreamerInceptionTagger(const RefCountedPtr<AmrMesh>&            a_amrMesh,
                              const EBAMRCellData* const               a_electricField,
                              const std::function<Real(const Real E)>& a_alphaEff,
                              const phase::which_phase                 a_phase = phase::gas);

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
      mutable EBAMRCellData m_tracerField;

      /*!
	@brief Ionization coefficient
      */
      std::function<Real(const Real& E)> m_alphaEff;

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
      Real m_refAlpha;

      /*!
	@brief Max voltage
      */
      Real m_maxVoltage;

      /*!
	@brief Plot tracer field or not
      */
      bool m_plot;

      /*!
	@brief Compute alpha_eff*dx at max voltage.
      */
      virtual void
      computeTracerField() const noexcept;
    };
  } // namespace StreamerInception
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
