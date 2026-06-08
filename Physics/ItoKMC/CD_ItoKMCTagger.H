/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
  @file   CD_ItoKMCTagger.H
  @brief  Declaration of the Physics::ItoKMC::ItoKMCTagger abstract CellTagger
  @author Robert Marskar
*/

#ifndef CD_ITOKMCTAGGER_H
#define CD_ITOKMCTAGGER_H

// Our includes
#include <CD_CellTagger.H>
#include <CD_ItoKMCStepper.H>
#include <CD_ItoKMCPhysics.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace ItoKMC {

    /**
      @brief Abstract CellTagger for the Ito-KMC-Poisson physics module.
      @details The template parameter @p S is the concrete ItoKMCStepper type. Subclasses implement
      `computeTagFields()`, `refineCell()`, and `coarsenCell()` to define the refinement criterion.
      Tag fields are computed on the Primal realm and live on the phase set by `setPhase()`.
    */
    template <typename S>
    class ItoKMCTagger : public CellTagger
    {
    public:
      /**
	@brief Weak constructor. User MUST subsequently call the define function. 
      */
      ItoKMCTagger();

      /**
	@brief Full constructor
	@param[in] a_physics Physics implementation. 
	@param[in] a_timeStepper Time integrator. 
	@param[in] a_amr Reference to AmrMesh 
      */
      ItoKMCTagger(const RefCountedPtr<ItoKMCPhysics>& a_physics,
                   const RefCountedPtr<S>&             a_timeStepper,
                   const RefCountedPtr<AmrMesh>&       a_amr) noexcept;

      /**
	@brief Destructor.
      */
      virtual ~ItoKMCTagger() noexcept;

      /**
	@brief Define function
	@param[in] a_physics Physics implementation. 
	@param[in] a_timeStepper Time integrator. 
	@param[in] a_amr Reference to AmrMesh 
      */
      virtual void
      define(const RefCountedPtr<ItoKMCPhysics>& a_physics,
             const RefCountedPtr<S>&             a_timeStepper,
             const RefCountedPtr<AmrMesh>&       a_amr) noexcept;

      /**
	@brief Get number of plot variables that will be written to file (by Driver).
	@return Returns number of plot variables that Driver will write to plot files. 
      */
      virtual int
      getNumberOfPlotVariables() const noexcept override;

      /**
	@brief Get plot variable names.
	@return Vector of variable name strings for the tag fields.
      */
      virtual Vector<std::string>
      getPlotVariableNames() const noexcept override;

      /**
	@brief Write plot data.
	@param[in,out] a_output Output data holder
	@param[in,out] a_icomp Starting variable in a_output where we begin appending data.
	@param[in] a_level Grid level
      */
      virtual void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string     a_outputRealm,
                    const int             a_level) const noexcept override;

      /**
	@brief Compute tagging fields.
      */
      virtual void
      computeTagFields() const noexcept = 0;

      /**
	@brief Parse class options
      */
      virtual void
      parseOptions() override = 0;

      /**
	@brief Parse run-time configurable class options
      */
      virtual void
      parseRuntimeOptions() override = 0;

      /**
	@brief Regrid this class. Note that there is no preRegrid method.
      */
      virtual void
      regrid() noexcept override;

      /**
	@brief Tag cells for refinement and coarsening
	@param[in,out] a_tags Cell tags that become the new grids. 
      */
      virtual bool
      tagCells(EBAMRTags& a_tags) noexcept override;

      /**
	@brief Set the phase where we do the tagging.
	@param[in] a_phase Phase
      */
      virtual void
      setPhase(const phase::which_phase a_phase) noexcept;

      /**
	@brief Get the tagging fields (mutable reference).
	@return Reference to the vector of cell-centred tagging fields on the Primal realm.
      */
      virtual Vector<EBAMRCellData>&
      getTagFields() noexcept;

    protected:
      /**
	@brief Defined or not
      */
      bool m_isDefined;

      /**
	@brief Plot data or not
      */
      bool m_plotData;

      /**
	@brief Realm where data is allocated.
      */
      std::string m_realm;

      /**
	@brief Phase where data is allocated.
      */
      phase::which_phase m_phase;

      /**
	@brief Physics implementation. 
      */
      RefCountedPtr<ItoKMCPhysics> m_physics;

      /**
	@brief Time stepper
      */
      RefCountedPtr<S> m_timeStepper;

      /**
	@brief AmrMesh. Needed for allocating memory and things like that. 
      */
      RefCountedPtr<AmrMesh> m_amr;

      /**
	@brief Number of tagging fields. 
      */
      int m_numTagFields;

      /**
	@brief Tagging fields. 
      */
      mutable Vector<EBAMRCellData> m_tagFields;

      /**
	@brief Gradient of tagging fields
      */
      mutable Vector<EBAMRCellData> m_gradTagFields;

      /**
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

      /**
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

      /**
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
  } // namespace ItoKMC
} // namespace Physics

#include <CD_NamespaceFooter.H>

#include <CD_ItoKMCTaggerImplem.H>

#endif
