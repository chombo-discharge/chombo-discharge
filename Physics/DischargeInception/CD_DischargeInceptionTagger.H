/*
 * SPDX-FileCopyrightText: 2022-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_DischargeInceptionTagger.H
   @brief  Declaration of the Physics::DischargeInception::DischargeInceptionTagger cell tagger
   @author Robert Marskar
*/

#ifndef CD_DISCHARGEINCEPTIONTAGGER_H
#define CD_DISCHARGEINCEPTIONTAGGER_H

// Our includes
#include <CD_CellTagger.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace DischargeInception {

    /**
       @brief CellTagger for the DischargeInception module.
       @details Flags cells for refinement where the effective ionization coefficient scaled by the
       local grid spacing exceeds a threshold: @f$ \alpha_\mathrm{eff} \cdot \Delta x > \text{threshold} @f$.
       The effective ionization coefficient @f$ \alpha_\mathrm{eff}(E, \mathbf{x}) @f$ is supplied by
       the caller and is typically the net Townsend coefficient @f$ \alpha - \eta @f$.
    */
    class DischargeInceptionTagger : public CellTagger
    {
    public:
      /**
         @brief Disallowed constructor
      */
      DischargeInceptionTagger() = delete;

      /**
         @brief Disallowed constructor
      */
      DischargeInceptionTagger(const DischargeInceptionTagger&) = delete;

      /**
         @brief Disallowed constructor
      */
      DischargeInceptionTagger(const DischargeInceptionTagger&&) = delete;

      /**
         @brief Full constructor.
         @param[in] a_amrMesh       AMR mesh.
         @param[in] a_electricField Pointer to the cell-centred electric field data (not owned).
         @param[in] a_alphaEff      Net effective ionization coefficient as a function of field magnitude E and position x.
         @param[in] a_phase         Phase (gas/solid) on which the tagger operates.
      */
      DischargeInceptionTagger(const RefCountedPtr<AmrMesh>&                              a_amrMesh,
                               const EBAMRCellData* const                                 a_electricField,
                               const std::function<Real(const Real E, const RealVect x)>& a_alphaEff,
                               const phase::which_phase                                   a_phase = phase::gas);

      /**
         @brief Destructor
      */
      virtual ~DischargeInceptionTagger();

      /**
         @brief Parse options
      */
      virtual void
      parseOptions() override final;

      /**
         @brief Regrid hook. Re-allocates internal data holders after the AMR hierarchy changes.
      */
      virtual void
      regrid() override final;

      /**
         @brief Tag cells for refinement. 
         @details This tags based on |grad(E)|*dx/E. 
         @param[in,out] a_tags Tags on grid levels.
      */
      virtual bool
      tagCells(EBAMRTags& a_tags) override final;

      /**
         @brief Get number of plot variables that will be written to file (by Driver).
         @return Returns number of plot variables that Driver will write to plot files. 
      */
      virtual int
      getNumberOfPlotVariables() const override;

      /**
         @brief Write plot data.
         @param[in,out] a_output            Output data holder
         @param[in,out] a_plotVariableNames Plot variable names
         @param[in,out] a_icomp             Starting variable in a_output where we begin appending data.
      */
      virtual void
      writePlotData(EBAMRCellData& a_output, Vector<std::string>& a_plotVariableNames, int& a_icomp) const;

    protected:
      /**
         @brief AmrMesh reference. 
      */
      RefCountedPtr<AmrMesh> m_amr;

      /**
         @brief Pointer to electric field.
      */
      const EBAMRCellData* m_electricField;

      /**
         @brief Tracer field. Will be populated with |grad(E)|*dx/|E|
      */
      mutable EBAMRCellData m_tracerField;

      /**
         @brief Net effective ionization coefficient @f$ \alpha_\mathrm{eff}(E, x) @f$ callable.
      */
      std::function<Real(const Real& E, const RealVect x)> m_alphaEff;

      /**
         @brief Phase (gas/solid) on which the tagger operates.
      */
      phase::which_phase m_phase;

      /**
         @brief Realm name where solver data lives.
      */
      std::string m_realm;

      /**
         @brief Threshold for refinement: refine if @f$ \alpha_\mathrm{eff} \cdot \Delta x > m\_refAlpha @f$.
      */
      Real m_refAlpha;

      /**
         @brief Reference voltage used to scale the electric field when computing the tracer field.
      */
      Real m_maxVoltage;

      /**
         @brief If true, write the tracer field to the plot file.
      */
      bool m_plot;

      /**
         @brief Compute alpha_eff*dx at max voltage.
      */
      virtual void
      computeTracerField() const noexcept;
    };
  } // namespace DischargeInception
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
