/* chombo-discharge
 * Copyright © 2021 SINTEF Energy Research.
 * Please refer to Copyright.txt and LICENSE in the chombo-discharge root directory.
 */

/*!
  @file   CD_ItoPlasmaTagger.H
  @brief  Parent tagger class for ItoPlasma
  @author Robert Marskar
*/

#ifndef CD_ItoPlasmaTagger_H
#define CD_ItoPlasmaTagger_H

// Our includes
#include <CD_CellTagger.H>
#include <CD_ItoPlasmaStepper.H>
#include <CD_ItoPlasmaPhysics.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoPlasma {

    /*!
      @brief Class for tagging grid cells using the Ito-KMC-Poisson module. Note that S is the time integrator
      type that is used. 
    */
    template <typename S>
    class ItoPlasmaTagger : public CellTagger
    {
    public:
      /*!
	@brief Weak constructor. User MUST subsequently call the define function. 
      */
      ItoPlasmaTagger();

      /*!
	@brief Full constructor
	@param[in] a_physics Physics implementation. 
	@param[in] a_timeStepper Time integrator. 
	@param[in] a_amr Reference to AmrMesh 
      */
      ItoPlasmaTagger(const RefCountedPtr<ItoPlasmaPhysics>& a_physics,
                      const RefCountedPtr<S>&                a_timeStepper,
                      const RefCountedPtr<AmrMesh>&          a_amr) noexcept;

      /*!
	@brief Destructor.
      */
      virtual ~ItoPlasmaTagger() noexcept;

      /*!
	@brief Define function
	@param[in] a_physics Physics implementation. 
	@param[in] a_timeStepper Time integrator. 
	@param[in] a_amr Reference to AmrMesh 
      */
      virtual void
      define(const RefCountedPtr<ItoPlasmaPhysics>& a_physics,
             const RefCountedPtr<S>&                a_timeStepper,
             const RefCountedPtr<AmrMesh>&          a_amr) noexcept;

      /*!
	@brief Compute tagging fields.
      */
      virtual void
      computeTagFields() const noexcept = 0;

      /*!
	@brief Parse class options
      */
      virtual void
      parseOptions() = 0;

      /*!
	@brief Parse run-time configurable class options
      */
      virtual void
      parseRuntimeOptions() = 0;

      /*!
	@brief Get number of plot variables for plot files.
      */
      virtual int
      getNumberOfPlotVariables() const noexcept override;

      /*!
	@brief Regrid this class. Note that there is no preRegrid method.
      */
      virtual void
      regrid() noexcept override;

      /*!
	@brief Tag cells for refinement and coarsening
	@param[inout] a_tags Cell tags that become the new grids. 
      */
      virtual bool
      tagCells(EBAMRTags& a_tags) noexcept override;

      /*!
	@brief Set the phase where we do the tagging.
	@param[in] a_phase Phase
      */
      virtual void
      setPhase(const phase::which_phase a_phase) noexcept;

      /*!
	@brief Get tag fields. This is data that live on the primal realm and specified input phase. 
      */
      virtual Vector<EBAMRCellData>&
      getTagFields() noexcept;

      /*!
	@brief Write plot data
	@param[inout] a_output            Output data holder
	@param[inout] a_plotVariableNames Plot variable names
	@param[inout] a_icomp             Starting variable in a_output where we begin appending data.
      */
      virtual void
      writePlotData(EBAMRCellData&       a_output,
                    Vector<std::string>& a_plotVariableNames,
                    int&                 a_icomp) const noexcept override;

    protected:
      /*!
	@brief Defined or not
      */
      bool m_isDefined;

      /*!
	@brief Plot data or not
      */
      bool m_plotData;

      /*!
	@brief Realm where data is allocated.
      */
      std::string m_realm;

      /*!
	@brief Phase where data is allocated.
      */
      phase::which_phase m_phase;

      /*!
	@brief Physics implementation. 
      */
      RefCountedPtr<ItoPlasmaPhysics> m_physics;

      /*!
	@brief Time stepper
      */
      RefCountedPtr<S> m_timeStepper;

      /*!
	@brief AmrMesh. Needed for allocating memory and things like that. 
      */
      RefCountedPtr<AmrMesh> m_amr;

      /*!
	@brief Number of tagging fields. 
      */
      int m_numTagFields;

      /*!
	@brief Tagging fields. 
      */
      mutable Vector<EBAMRCellData> m_tagFields;

      /*!
	@brief Gradient of tagging fields
      */
      mutable Vector<EBAMRCellData> m_gradTagFields;

      /*!
	@brief Per-box refinement tagging routine. 
	@param[out] a_refinedCells Cells to be refined
	@param[out] a_coarsenedCells Cells to be coarsened
	@param[in] a_tagFields Tagging fields
	@param[in] a_gradTagFields Gradient of tagging fields
	@param[in] a_lvl Grid level
	@param[in] a_dit Grid index
	@param[in] a_box Cell-centered box
	@param[in] a_ebisbox EBIS box
	@param[in] a_time Current time
	@param[in] a_dx Grid resolution
	@param[in] a_probLo Lower-left corner of physical domain
      */
      virtual void
      tagCellsBox(DenseIntVectSet&          a_refinedCells,
                  DenseIntVectSet&          a_coarsenedCells,
                  const Vector<EBCellFAB*>& a_tagFields,
                  const Vector<EBCellFAB*>& a_gradTagFields,
                  const int                 a_lvl,
                  const DataIndex           a_dit,
                  const Box                 a_box,
                  const EBISBox&            a_ebisbox,
                  const Real                a_time,
                  const Real                a_dx,
                  const RealVect            a_probLo) const noexcept;

      /*!
	@brief Determine if a particular cell should be coarsened or not. 
	@param[in] a_pos Physical position. 
	@param[in] a_time Current time
	@param[in] a_dx Grid resolution
	@param[in] a_lvl Grid level 
	@param[in] a_tagFields Cell tagging fields
	@param[in] a_gradTagFields Gradient of cell tagging fields
      */
      virtual bool
      coarsenCell(const RealVect         a_pos,
                  const Real             a_time,
                  const Real             a_dx,
                  const int              a_lvl,
                  const Vector<Real>     a_tagFields,
                  const Vector<RealVect> a_gradTagFields) const noexcept = 0;

      /*!
	@brief Determine if a particular cell should be refined or not. 
	@param[in] a_pos Physical position. 
	@param[in] a_time Current time
	@param[in] a_dx Grid resolution
	@param[in] a_lvl Grid level 
	@param[in] a_tagFields Cell tagging fields
	@param[in] a_gradTagFields Gradient of cell tagging fields
      */
      virtual bool
      refineCell(const RealVect         a_pos,
                 const Real             a_time,
                 const Real             a_dx,
                 const int              a_lvl,
                 const Vector<Real>     a_tagFields,
                 const Vector<RealVect> a_gradTagFields) const noexcept = 0;
    };
  } // namespace ItoPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoPlasmaTaggerImplem.H>

#endif
