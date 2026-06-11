/*
 * SPDX-FileCopyrightText: 2021-2026 SINTEF Energy Research
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

/**
   @file   CD_CdrPlasmaTagger.H
   @brief  Declaration of the Physics::CdrPlasma::CdrPlasmaTagger abstract class
   @author Robert Marskar
*/

#ifndef CD_CDRPLASMATAGGER_H
#define CD_CDRPLASMATAGGER_H

// Our includes
#include <CD_CellTagger.H>
#include <CD_CdrPlasmaStepper.H>
#include <CD_NamespaceHeader.H>

namespace Physics {
  namespace CdrPlasma {

    /**
       @brief Abstract CellTagger base class for the CDR plasma physics module.
       @details This class provides subclasses with access to solver data and plasma physics
       through references to CdrPlasmaStepper, AmrMesh, CdrPlasmaPhysics, and
       ComputationalGeometry. It implements the CellTagger interface using a tracer-field
       approach: subclasses compute one or more scalar tracer fields via computeTracers(),
       and the base class then calls the per-cell refineCell() and coarsenCell() methods to
       determine which cells should be flagged.

       Subclasses must implement computeTracers(), parseOptions(), refineCell(), and coarsenCell().
       The number of tracer fields is controlled by setting m_numTracers in the subclass constructor.
    */
    class CdrPlasmaTagger : public CellTagger
    {
    public:
      /**
         @brief Default constructor. Sets default realm and phase.
      */
      CdrPlasmaTagger();

      /**
         @brief Full constructor. Leaves the object in a usable state.
         @param[in] a_physics               Plasma kinetics interface.
         @param[in] a_timeStepper           Time stepper providing solver data.
         @param[in] a_amr                   AMR mesh for grid iteration and operators.
         @param[in] a_computationalGeometry Computational geometry (electrodes and dielectrics).
      */
      CdrPlasmaTagger(const RefCountedPtr<CdrPlasmaPhysics>&      a_physics,
                      const RefCountedPtr<CdrPlasmaStepper>&      a_timeStepper,
                      const RefCountedPtr<AmrMesh>&               a_amr,
                      const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

      /**
         @brief Destructor.
      */
      virtual ~CdrPlasmaTagger();

      /**
         @brief Define function. Leaves the object in a usable state.
         @param[in] a_physics               Plasma kinetics interface.
         @param[in] a_timeStepper           Time stepper providing solver data.
         @param[in] a_amr                   AMR mesh for grid iteration and operators.
         @param[in] a_computationalGeometry Computational geometry.
      */
      void
      define(const RefCountedPtr<CdrPlasmaPhysics>&      a_physics,
             const RefCountedPtr<CdrPlasmaStepper>&      a_timeStepper,
             const RefCountedPtr<AmrMesh>&               a_amr,
             const RefCountedPtr<ComputationalGeometry>& a_computationalGeometry);

      /**
         @brief Compute tracer fields on the AMR hierarchy. Must be implemented by subclasses.
         @details Results are stored in m_tracers. The number of tracer fields must equal m_numTracers.
      */
      virtual void
      computeTracers() const = 0;

      /**
         @brief Parse class options from ParmParse. Must be implemented by subclasses.
      */
      virtual void
      parseOptions() override = 0;

      /**
         @brief Pre-plot hook. Computes tracer fields so they appear in plot output.
      */
      virtual void
      prePlot() const noexcept override;

      /**
         @brief Pre-regrid hook. Clears tracer field storage before the grid changes.
      */
      virtual void
      preRegrid() noexcept override;

      /**
         @brief Regrid hook. Reallocates tracer field and gradient storage on the new mesh.
      */
      virtual void
      regrid() override final;

      /**
         @brief Flag cells for refinement and coarsening using the tracer-field approach.
         @details Calls computeTracers(), then iterates over grid cells calling refineCell()
         and coarsenCell() on each. Returns true if the tag set changed on any MPI rank.
         @param[in,out] a_tags Cell tags across the AMR hierarchy.
         @return True if any new tags were added or removed.
      */
      virtual bool
      tagCells(EBAMRTags& a_tags) override final;

      /**
         @brief Return the number of plot variables (one per tracer field).
         @return Number of tracer fields (m_numTracers).
      */
      virtual int
      getNumberOfPlotVariables() const override final;

      /**
         @brief Return the names of the tracer field plot variables.
         @return Vector of names of the form "Tracer field-i".
      */
      virtual Vector<std::string>
      getPlotVariableNames() const override final;

      /**
         @brief Write tracer field data into the output holder.
         @param[in,out] a_output      Output data holder.
         @param[in,out] a_icomp       Starting component index; incremented on return.
         @param[in]    a_outputRealm Realm that owns a_output.
         @param[in]    a_level       AMR level to write.
      */
      virtual void
      writePlotData(LevelData<EBCellFAB>& a_output,
                    int&                  a_icomp,
                    const std::string&    a_outputRealm,
                    const int             a_level) const override final;

    protected:
      /**
         @brief Realm on which tracer field storage is allocated.
      */
      std::string m_realm;

      /**
         @brief Phase on which tracer field storage is allocated.
      */
      phase::which_phase m_phase;

      /**
         @brief Plasma kinetics interface.
      */
      RefCountedPtr<CdrPlasmaPhysics> m_physics;

      /**
         @brief Time stepper providing access to solver data and electric field.
      */
      RefCountedPtr<CdrPlasmaStepper> m_timeStepper;

      /**
         @brief AMR mesh for grid iteration, gradient computation, and interpolation.
      */
      RefCountedPtr<AmrMesh> m_amr;

      /**
         @brief Computational geometry providing electrode and dielectric information.
      */
      RefCountedPtr<ComputationalGeometry> m_computationalGeometry;

      /**
         @brief Number of scalar tracer fields used for tagging decisions.
      */
      int m_numTracers;

      /**
         @brief Tracer fields on the AMR hierarchy. Allocated and filled by computeTracers().
      */
      mutable Vector<EBAMRCellData> m_tracers;

      /**
         @brief Gradients of the tracer fields on the AMR hierarchy.
      */
      mutable Vector<EBAMRCellData> m_gradTracers;

      /**
         @brief Flag cells for refinement within a single grid patch.
         @param[out] a_refinedCells Cells flagged for refinement on this patch.
         @param[in]  a_tracers      Tracer field data on this patch.
         @param[in]  a_gradTracers  Gradient of tracer fields on this patch.
         @param[in]  a_lvl          AMR level.
         @param[in]  a_dit          Grid patch index.
         @param[in]  a_box          Cell-centered box for this patch.
         @param[in]  a_ebisbox      EB geometry information for this patch.
         @param[in]  a_time         Current simulation time.
         @param[in]  a_dx           Grid cell spacing.
         @param[in]  a_probLo       Lower-left corner of the simulation domain.
      */
      virtual void
      refineCellsBox(DenseIntVectSet&          a_refinedCells,
                     const Vector<EBCellFAB*>& a_tracers,
                     const Vector<EBCellFAB*>& a_gradTracers,
                     const int                 a_lvl,
                     const DataIndex&          a_dit,
                     const Box&                a_box,
                     const EBISBox&            a_ebisbox,
                     const Real                a_time,
                     const Real                a_dx,
                     const RealVect&           a_probLo);

      /**
         @brief Flag cells for coarsening within a single grid patch.
         @param[out] a_coarsenedCells Cells flagged for coarsening on this patch.
         @param[in]  a_tracers        Tracer field data on this patch.
         @param[in]  a_gradTracers    Gradient of tracer fields on this patch.
         @param[in]  a_lvl            AMR level.
         @param[in]  a_dit            Grid patch index.
         @param[in]  a_box            Cell-centered box for this patch.
         @param[in]  a_ebisbox        EB geometry information for this patch.
         @param[in]  a_time           Current simulation time.
         @param[in]  a_dx             Grid cell spacing.
         @param[in]  a_probLo         Lower-left corner of the simulation domain.
      */
      virtual void
      coarsenCellsBox(DenseIntVectSet&          a_coarsenedCells,
                      const Vector<EBCellFAB*>& a_tracers,
                      const Vector<EBCellFAB*>& a_gradTracers,
                      const int                 a_lvl,
                      const DataIndex&          a_dit,
                      const Box&                a_box,
                      const EBISBox&            a_ebisbox,
                      const Real                a_time,
                      const Real                a_dx,
                      const RealVect&           a_probLo);

      /**
         @brief Decide whether a single cell should be coarsened. Must be implemented by subclasses.
         @param[in] a_pos         Cell centroid position.
         @param[in] a_time        Current simulation time.
         @param[in] a_dx          Grid cell spacing.
         @param[in] a_lvl         AMR level of this cell.
         @param[in] a_tracers     Tracer field values at this cell.
         @param[in] a_gradTracers Gradient of tracer fields at this cell.
         @return True if the cell should be coarsened.
      */
      virtual bool
      coarsenCell(const RealVect&         a_pos,
                  const Real              a_time,
                  const Real              a_dx,
                  const int               a_lvl,
                  const Vector<Real>&     a_tracers,
                  const Vector<RealVect>& a_gradTracers) const = 0;

      /**
         @brief Decide whether a single cell should be refined. Must be implemented by subclasses.
         @param[in] a_pos         Cell centroid position.
         @param[in] a_time        Current simulation time.
         @param[in] a_dx          Grid cell spacing.
         @param[in] a_lvl         AMR level of this cell.
         @param[in] a_tracers     Tracer field values at this cell.
         @param[in] a_gradTracers Gradient of tracer fields at this cell.
         @return True if the cell should be refined.
      */
      virtual bool
      refineCell(const RealVect&         a_pos,
                 const Real              a_time,
                 const Real              a_dx,
                 const int               a_lvl,
                 const Vector<Real>&     a_tracers,
                 const Vector<RealVect>& a_gradTracers) const = 0;
    };
  } // namespace CdrPlasma
} // namespace Physics

#include <CD_NamespaceFooter.H>

#endif
